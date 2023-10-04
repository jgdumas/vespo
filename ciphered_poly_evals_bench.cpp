// ==========================================================================
// Evaluation of LHE enciphered Polynomial on a public oint
// Authors: J-G Dumas
// Time-stamp: <04 Oct 23 10:18:54 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/****************************************************************
 * Uses libvespo https://github.com/jgdumas/vespo
 * Refzrence: "VESPo: Verified Evaluation of Secret Polynomials"
        PETS 2023: 23rd Privacy Enhancing Technologies Symposium,
        J-G. Dumas, A. Maignan, C. Pernet, D. S. Roche,
        https://petsymposium.org/popets/2023/popets-2023-0085.phpV
 ****************************************************************/

#include "vespo_library.h"
#include "vespo_library.inl"

//========================================================
// Main
//========================================================
int main(int argc, char * argv[]) {
        //argv[1]: Degree>1 of the polynomial
        //argv[2]: modulus size in bits for Paillier
        //argv[3]: number of experiments
        //argv[4]: number of parallel tasks (default is number of threads)

    if (argc < 1) {
        std::cerr << "Usage: " << argv[0] << " <Polynomial degree> <Paillier bitsize> <Verification iterations> [tasks (default=threads)]" << std::endl;
        return 1;
    }

    std::srand(std::time(nullptr));
    int64_t numthreads(1);

#if defined(_OPENMP)
#pragma omp parallel
#pragma omp single
    {
        numthreads = omp_get_num_threads();
    }
#endif

		// Polynomial degree
    const int64_t degree = std::max(0,atoi(argv[1]));
    const int64_t psize(degree+1);
        // Paillier bitsize
    const uint64_t pailliersize = (argc>2?atoi(argv[2]):2048);
		// Number of evaluation points
    const uint64_t nb_points = std::max(1,(argc>3?atoi(argv[3]):3));
        // tasks
    const int nb_tasks = (argc>4?atoi(argv[4]):std::min(numthreads,psize));

        /********************************************************************
         * RELIC SETUP
         *******************************************************************/
    if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}
	conf_print();

	if (pc_param_set_any() != RLC_OK) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 1;
	}
	pc_param_print();
	std::clog << "  security level: " << pc_param_level() << std::endl;

        /********************************************************************
         * VESPo: parameters
         *******************************************************************/
	util_banner("VESPO:", 0);

    std::clog << "  MAX_ALLOC: " << VESPO_RELIC_LIMIT_MAX_ALLOC << std::endl;
    std::clog << "  SHAMIR_WD: " << RLC_MAX(8, RLC_WIDTH) << std::endl;
    std::clog << "  CLOCK TYP: " << CLOCKTYPE << std::endl;
	std::clog << "  OMP CORES: " << numthreads << std::endl;
	std::clog << "  NUM POINT: " << nb_points << std::endl;
	std::clog << "  NUM TASKS: " << nb_tasks << std::endl;

#ifdef VESPO_NOTSECURE
    std::clog << "  REUSEDPRD: " << VESPO_NOTSECURE << std::endl;
#endif

    double time_i(0.), time_r(0.), time_d(0.), time_e(0.);
    Chrono c_step; c_step.start();

        /********************************************************************
         * VESPo: group order
         *******************************************************************/
    bn_t group_mod; bn_null(group_mod); bn_new(group_mod); pc_get_ord(group_mod);
	std::clog << "  MOD BITS : " << bn_bits(group_mod) << std::endl;
	std::clog << "  PAIL SIZE: " << pailliersize << std::endl;

        /********************************************************************
         * VESPo: Benchmarking with a random polynomial
         *******************************************************************/
    Polynomial<bn_t> P(degree, group_mod); P.random();

#ifdef DEBUG
    std::clog << "[RANDOM P]: " << P << std::endl;
#endif


    c_step.start();


    paillier_pubkey_t pub;  // Paillier keys
    paillier_prvkey_t prv;  //     public/private
    paillier_newpubkey(pub);
    paillier_newprvkey(prv);

    pc_get_ord( prv->a );
    cp_shpe_gen( pub.rlc, prv, bn_bits( prv->a ), pailliersize);
    bn_sqr(pub.nsq, pub.rlc->crt->n); // store n^2 too

        // Enciphered Polynomial, cut into 'blocks' chunks
        //     last ones of size [ (degree+1)/Blocks ]
        //     first ones of that size + 1
        //     so that sum of sizes is 'degree+1'
        //     W is the vector of these blocks
    std::vector<Polynomial<paillier_ciphertext_t>> W;
    const int64_t blocks = std::max(2,nb_tasks);
    const int64_t bigblocks = setup_block_chunks(W, degree, blocks, group_mod);

    std::clog << "[Setup Paillier mod "<< pub.rlc->crt->n << ']' << std::endl;

        // parallel store W on server
    encrypt_poly(W, P, pub, prv);
        // parallel hide the polynomial

    time_i = c_step.stop();
    std::clog << "[TIME SETUP CIPHER P]: " << time_i << " (" <<
#ifndef VESPO_NOTSECURE
        psize
#else
        std::min((int64_t)VESPO_NOTSECURE,psize)
#endif
              << " randomness)" << std::endl;
    c_step.start();

		// Evaluation points
    Polynomial<bn_t> r(nb_points-1, group_mod); r.random();

        // The correct polynomial evaluations
    Polynomial<bn_t> Pr(nb_points-1, group_mod);
    for(size_t l=0; l<nb_points; ++l)
        P.eval(Pr[l], r[l]);

    time_r = c_step.stop();
    std::clog << "[TIME Ref EVALUATION]: " << time_r << " (" << nb_points << " Horner evaluations)" << std::endl;
    c_step.start();


        // Enciphered polynomial evaluation
    paillier_plaintext_t z;

    const int64_t nbblocks(W.size());
    const int64_t dblocks(nbblocks-1);

    paillier_ciphertext_t zeta;
    Polynomial<paillier_ciphertext_t> bzeta(dblocks, group_mod);

    bn_t tmpz; bn_new(tmpz);
    bn_t rpowsm; bn_new(rpowsm);
    bn_t rpows; bn_new(rpows);

    const int64_t step( std::ceil( (double)degree/(double)nb_tasks) );
    const int64_t powsr_deg(std::max(W[0].degree(),step));
    const int64_t apows_deg(std::max(powsr_deg, std::max(bigblocks,dblocks)));

#ifdef DEBUG
    std::clog << "pows_r: d°" << powsr_deg << '|' << apows_deg << std::endl;
#endif

    Polynomial<bn_t> pows_r(apows_deg,group_mod);


    size_t correct_geomprog(0), correct_evals(0);


    for(size_t l=0; l<nb_points; ++l) {

        geo_progression(pows_r, 0, r[l], powsr_deg, r[l], group_mod, time_r);

        bn_copy(rpowsm, pows_r[W[0].degree()]);	// last power of r;
        bn_mul(rpows, rpowsm, r[l]);			// next power of r;
        bn_mod(rpows, rpows, group_mod);

        Polynomial<bn_t> revpow(step,group_mod);
        for(int64_t i=0; i<(step+1); ++i) {
            bn_copy(revpow[i],pows_r[step-i]);
        }
            // compute evaluation by blocks
#pragma omp parallel for shared(bzeta,pub,W,pows_r)
        for(int64_t i=0; i<nbblocks; ++i)
                // server computation : zeta baby-steps
            paillier_hom_dp(bzeta[i], pub, W[i], pows_r, W[i].degree());

        // server computation : zeta giant-step
        geo_progression(pows_r, bigblocks, rpows, dblocks, rpowsm,
                        group_mod, time_r);

        paillier_hom_dp(zeta, pub, bzeta, pows_r, dblocks);

        time_e += c_step.stop();

            // Is the geometric progression correct ?
        bn_mxp_dig(tmpz, r[l], W[0].degree()+1, group_mod);
        if(bn_cmp(tmpz, rpows) == RLC_EQ) ++correct_geomprog;

        c_step.start();

            // Client: Deciphering Paillier to get z
        PAILLIER_DECIPHER(z, pub, prv, zeta);
        bn_mod(z.m,z.m,group_mod);


        time_d += c_step.stop();

        // Is the evaluation correct ?
        if(bn_cmp(Pr[l],z.m) == RLC_EQ) ++correct_evals;
    }


    if(correct_geomprog == nb_points) {
        std::clog << "[CLEARTEXT GEOM-PROG] \033[1;32mOK\033[0m ";
    } else {
        std::cerr << "[CLEARTEXT GEOM-PROG] \033[1;31m***** FAIL *****\033[0m ";
    }
    std::clog << correct_geomprog << '/' << nb_points << std::endl;


    if(correct_evals == nb_points) {
        std::clog << "[CIPHERED EVALUATION] \033[1;32mOK\033[0m ";
    } else {
        std::cerr << "[CIPHERED EVALUATION] \033[1;31m***** FAIL *****\033[0m ";
    }
    std::clog << correct_evals << '/' << nb_points << std::endl;


    std::clog << "[TIME CIPHERED EVALS]: " << time_e << " (" << (nb_points*psize) << " Paillier operations on ciphers)" << std::endl;
    std::clog << "[TIME DECIPHERING Z ]: " << time_d << " (" << nb_points << " Paillier decipherings)" << std::endl;


        /********************************************************************
         * VESPo: clean-up
         *******************************************************************/
    bn_free(group_mod);
    bn_free(tmpz);
    bn_free(rpowsm);
    bn_free(rpows);
    paillier_freeprvkey(this->prv);
    paillier_freepubkey(this->pub);

	core_clean();

	return 0;
}
