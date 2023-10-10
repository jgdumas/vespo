// ==========================================================================
// VESPo: Verified Evaluation of Secret Polynomials
// Reference: [ https://arxiv.org/abs/2110.02022
//              J-G. Dumas, A. Maignan, C. Pernet, D. S. Roche ]
// Authors: J-G Dumas
// Time-stamp: <04 Oct 23 09:40:05 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/****************************************************************
 * Example of parameters, see file: preset_x64-pbc-bn254_vespo.sh
 * Or, parameters RELIC via: ccmake relic-target
 * ARITH                            x64-asm-4l	(or gmp)
 * BN_PRECI                         4096
 * FP_PRIME                         254			(or 256)
 * see relic-git-src/cmake/bn.cmake: set(BN_PRECI 4096)
 * see relic-git-src/cmake/fp.cmake: set(FP_PRIME 254) or 256 ...
 ****************************************************************/

#include "vespo_library.h"
#include "vespo_library.inl"

//========================================================
// Main
//========================================================
int main(int argc, char * argv[]) {
	//argv[1]: size of the vector (polynomial degree + 1)
	//argv[2]: modulus size in bits for Paillier
	//argv[3]: number of experiments
	//argv[4]: number of parallel tasks (default is number of threads)

    if (argc < 1) {
	std::cerr << "Usage: " << argv[0] << " <Vector dimension> <Paillier bitsize> <Verification iterations> [tasks (default=threads)]" << std::endl;
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
    const int64_t degree = std::max(0,atoi(argv[1]) - 1);
	// Paillier bitsize
    const uint64_t pailliersize = (argc>2?atoi(argv[2]):2048);
		// Verification iterations
    const uint64_t nb_iter = (argc>3?atoi(argv[3]):3);
	// tasks
    const int nb_tasks = (argc>4?atoi(argv[4]):std::min(numthreads,degree+1));

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

    const uint64_t group_bits = pc_param_level()<<1;
	std::clog << "  security level: " << pc_param_level() << std::endl;

	/********************************************************************
	 * VESPo: parameters
	 *******************************************************************/
	util_banner("VESPO:", 0);

    std::clog << "  MAX_ALLOC: " << VESPO_RELIC_LIMIT_MAX_ALLOC << std::endl;
    std::clog << "  SHAMIR_WD: " << RLC_MAX(8, RLC_WIDTH) << std::endl;
    std::clog << "  CLOCK TYP: " << CLOCKTYPE << std::endl;
	std::clog << "  OMP CORES: " << numthreads << std::endl;
	std::clog << "  NUM ITERS: " << nb_iter << std::endl;
	std::clog << "  NUM TASKS: " << nb_tasks << std::endl;

#ifdef VESPO_NOTSECURE
    std::clog << "  REUSEDPRD: " << VESPO_NOTSECURE << std::endl;
#endif

	/********************************************************************
	 * VESPo: group order
	 *******************************************************************/
    bn_t group_mod; bn_null(group_mod); bn_new(group_mod); pc_get_ord(group_mod);
	std::clog << "  MOD BITS : " << bn_bits(group_mod) << std::endl;
	std::clog << "  PAIL SIZE: " << pailliersize << std::endl;

    std::vector<double> time_c(nb_iter), time_s(nb_iter), time_u((nb_iter*(nb_iter+1))>>1);
    double time_i=0., time_eval=0.;

	/********************************************************************
	 * VESPo: Benchmarking with a random polynomial
	 *******************************************************************/
    Polynomial<bn_t> P(degree, group_mod); P.random();

#ifdef DEBUG
    std::clog << "[RANDOM P]: " << P << std::endl;
#endif

	/********************************************************************
	 * VESPo: BSetup client/server
	 *******************************************************************/
    Chrono c_setup; c_setup.start();

    client_t client(degree);
	// For now 2 threads for xi1 and xi2 the rest (at least 2) for Paillier
    server_t server(degree, std::max(2,nb_tasks), group_mod);

    setup(client, server, pailliersize, P, time_eval);

    time_i = c_setup.stop();

#ifdef VESPO_TIMINGS
    std::vector<double> time_e(nb_iter), time_p(nb_iter), time_ce(nb_iter), time_cg(nb_iter), time_cp(nb_iter);
    std::clog << "[CLIENT STATE]: " << client.private_size(group_bits) << " bits" << std::endl;
#endif

	/********************************************************************
	 * VESPo: Iterate nb_iter client/server audits
	 *******************************************************************/
    size_t updnum(0);
    bn_t r; bn_new(r);
    for(uint64_t i = 0; i < nb_iter; ++i) {

	    // Perform i random updates at random monomials
	for(uint64_t j = 0; j < i; ++j) {
	    bn_rand_mod(r, group_mod);
	    update(client, server, P, r, std::rand() % (degree+1), time_u[updnum++]);
	}

	    // Evaluation point
	bn_rand_mod(r, group_mod);

#ifdef DEBUG
	printf("r: "); bn_print(r); printf("\n");
#endif
#ifdef VESPO_CHECKERS
	    // The correct polynomial evaluation
	bn_t Pr; bn_new(Pr); P.eval(Pr, r);
#endif

	paillier_plaintext_t z;
	    // Verified evaluation
	bool pass = eval(z, client, server, r, nb_tasks, time_c[i], time_s[i]
#ifdef VESPO_TIMINGS
			 , time_e[i], time_p[i]
			 , time_ce[i], time_cg[i], time_cp[i]
#endif
			 );

	if (pass) {
	    std::clog << "[ VAUDIT ] " << (i+1) << '/' << nb_iter << " \033[1;32mOK : EVERYTHING IS WORKING !\033[0m" << std::endl;
	} else {
	    std::cerr << "[ VAUDIT ] " << (i+1) << '/' << nb_iter << " \033[1;31m****** FAIL there is a PROBLEM ******\033[0m" << std::endl;
	}

#ifdef VESPO_CHECKERS
	    // Audit has passed, but is the protocol sound?
	if(bn_cmp(Pr,z.m) == RLC_EQ) {
	    std::clog << "[SOUNDNESS] \033[1;32mOK\033[0m" << std::endl;
	} else {
	    std::cerr << "[SOUNDNESS] \033[1;31m****** FAIL ******\033[0m" << std::endl;
	}
	bn_free(Pr);
#endif
    }


	/********************************************************************
	 * VESPo: clean-up
	 *******************************************************************/
    bn_free(r);
    bn_free(group_mod);
	core_clean();

	/********************************************************************
	 * VESPo: Compute medians and deviations
	 *******************************************************************/
    std::sort(time_c.begin(),time_c.end());
    std::sort(time_s.begin(),time_s.end());
    std::sort(time_u.begin(),time_u.end());

#ifdef VESPO_TIMINGS
    std::sort(time_e.begin(),time_e.end());
    std::sort(time_p.begin(),time_p.end());
    std::sort(time_ce.begin(),time_ce.end());
    std::sort(time_cg.begin(),time_cg.end());
    std::sort(time_cp.begin(),time_cp.end());
    printf("[Eval DEVI] %lu | zeta : %f%% | xi : %f%% | C-powm : %f%% | C-gsum : %f%% | C-H_dec: %f%%\n", degree+1, mediandeviation(time_e), mediandeviation(time_p), mediandeviation(time_ce), mediandeviation(time_cg), mediandeviation(time_cp));
    printf("[Eval TIME] %lu | zeta : %f | xi : %f | C-powm : %f | C-gsum : %f | C-H_dec: %f | pure-horner : %f\n", degree+1, time_e[nb_iter>>1], time_p[nb_iter>>1], time_ce[nb_iter>>1], time_cg[nb_iter>>1], time_cp[nb_iter>>1], time_eval);
#endif

    printf("[DEVIATIO] | %lu %lu %lu | audit-client : %f%% | audit-server : %f%%\n", degree+1, group_bits, pailliersize, mediandeviation(time_c), mediandeviation(time_s));

    printf("[TIMINGS ] | %lu %lu %lu | setup : %f | audit-client : %f | audit-server : %f | pure-horner : %f | update : %f  \n=== end ===\n\n", degree+1, group_bits, pailliersize, time_i, time_c[time_c.size()>>1], time_s[time_s.size()>>1], time_eval, time_u[time_u.size()>>1]);

	return 0;
}
