// ==========================================================================
// VESPo: Verified Evaluation of Secret Polynomials
// Reference: [ https://arxiv.org/abs/2110.02022
//              J-G. Dumas, A. Maignan, C. Pernet, D. S. Roche ]
// Authors: J-G Dumas
// Time-stamp: <04 Oct 23 09:54:10 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/****************************************************************
 * VESPO Library inline implementations
 ****************************************************************/

#include "vespo_library.h"

/****************************************************************
 * Median deviation
 ****************************************************************/
template<typename Vect> inline double mediandeviation(const Vect& v) {
    assert(v.size()>0);
    typename Vect::value_type median(v[v.size()/2]);
    double t1( median-v.front() );
    double t2( v.back()-median );
    return 100.*std::max(t1,t2)/median;
}


/****************************************************************
 * relic bn_t utilities
 ****************************************************************/
inline std::ostream& operator<<(std::ostream& out, const bn_t& numb) {
    size_t bits = bn_size_str(numb, 10);
    char str[RLC_BN_BITS + 2];
    bn_write_str(str, bits, numb, 10);
    return out << str;
}

inline void vespo_init_set(bn_t& a, const bn_t b) {
    bn_null(a); bn_new(a); bn_copy(a,b);
}

inline void vespo_init_set_ui(bn_t& a, const dig_t b) {
    bn_null(a); bn_new(a); bn_set_dig(a,b);
}

// Set a to a + b * c.
inline void vespo_addinmul(bn_t& a, const bn_t b, const bn_t c) {
    bn_t tmp; bn_new(tmp);
    bn_mul(tmp,b,c);
    bn_add(a,a,tmp);
    bn_free(tmp);
}


//=====================================================================
// VeSPo data structures
//       Polynomial types
//=====================================================================
template<> inline void freeer<bn_t>(bn_t& e) { bn_free(e); }
template<> inline void freeer<g1_t>(g1_t& e) { g1_free(e); }
template<> inline void freeer<g2_t>(g2_t& e) { g2_free(e); }
template<> inline void freeer<gt_t>(gt_t& e) { gt_free(e); }
template<> inline void freeer<paillier_ciphertext_t>(paillier_ciphertext_t& e) { }

template<> inline void newer<bn_t>(bn_t& e) { bn_null(e); bn_new(e); }
template<> inline void newer<g1_t>(g1_t& e) { g1_null(e); g1_new(e); }
template<> inline void newer<g2_t>(g2_t& e) { g2_null(e); g2_new(e); }
template<> inline void newer<gt_t>(gt_t& e) { gt_null(e); gt_new(e); }
template<> inline void newer<paillier_ciphertext_t>(paillier_ciphertext_t& e) { }

template<> inline std::ostream& printer<bn_t>(std::ostream& o, const bn_t& e) { return o<<e;}
template<> inline std::ostream& printer<g1_t>(std::ostream& o, const g1_t& e) { g1_print(e); return o;}
template<> inline std::ostream& printer<g2_t>(std::ostream& o, const g2_t& e) { g2_print(e); return o;}
template<> inline std::ostream& printer<gt_t>(std::ostream& o, const gt_t& e) { gt_print(e); return o;}


template<typename Element>
inline Polynomial<Element>::Polynomial(const int64_t degree,
                                       const bn_t& modulus) :
        mod(modulus), d(degree), data(new Element[degree+1]) {
    for (int64_t i = 0; i <= degree; ++i) {
        newer(this->data[i]);
    }
}


template<typename Element>
inline Polynomial<Element>::~Polynomial() {
    for (int64_t i = 0; i <= this->d; ++i) {
        freeer(this->data[i]);
    }
    delete [] data;
}


template<typename Element>
inline const Element& Polynomial<Element>::operator[](size_t i) const {
#ifdef DEBUG
    if ( (int64_t)i > this->degree()) {
        std::cerr << "****************** ERROR ******************* Polynomial access at " << i << " out of " << this->degree() << std::endl;
        throw(1);
    }
#endif
    return this->data[i];
}

template<typename Element>
inline Element& Polynomial<Element>::operator[](size_t i) {
#ifdef DEBUG
    if ( (int64_t)i > this->degree()) {
        std::cerr << "****************** ERROR ******************* Polynomial access at " << i << " out of " << this->degree() << std::endl;
        throw(1);
    }
#endif
    return this->data[i];
}

template<typename Element>
inline std::ostream& operator<<(std::ostream& out, const Polynomial<Element>& P) {
    printer(out,P.data[0]) << " + ";
    for (int64_t i = 1; i < P.d; ++i) {
        printer(out,P.data[i]) << "X^" << i << " + ";
    }
    return printer(out,P.data[P.d]) << "X^" << P.d;
}


template<> inline void Polynomial<bn_t>::random() {
    std::clog << "[Rand P] BEG: " << this->d << ", mod " << this->modulus() << std::endl;
    for(int64_t i = 0; i <= this->d; ++i) {
        bn_rand_mod(this->data[i], this->modulus());
    }
    std::clog << "[Rand P] END" << std::endl;
}

template<> inline void Polynomial<bn_t>::eval(bn_t& eval, const bn_t& s) const {
#ifdef DEBUG
    std::clog << "[Horner] BEG: " << this->d << std::endl;
#endif
    bn_copy(eval,this->data[this->d]);

#ifdef DEBUG
    printf("mod : "); bn_print(this->modulus()); printf("\n");
    printf("P_%lu: ", this->d); bn_print(eval); printf("\n");
#endif

    for(int64_t i = this->d; i>0; --i) {	// Horner
        bn_mul(eval,eval,s);
        bn_add(eval,eval,this->data[i-1]);
        bn_mod(eval,eval, this->modulus());
    }

#ifdef DEBUG
    bn_print(eval); printf(" = P("); bn_print(s); printf(")\n");
    std::clog << "[Horner] END" << std::endl;
#endif
}

/****************************************************************
 * Paillier utilities
 ****************************************************************/
inline int paillier_mul(paillier_ciphertext_t& res,
                 const paillier_ciphertext_t& a,
                 const paillier_ciphertext_t& b,
                 const paillier_pubkey_t& kpub) {
	int result = RLC_OK;

    bn_mul(res.c, a.c, b.c);
    bn_mod(res.c, res.c, kpub.nsq);

	return result;
}

/****************************************************************
 * Paillier generation
 ****************************************************************/
//=====================================================================
// WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
// if VESPO_NOTSECURE is defined, Paillier masks are reused.
// Reusing the same randomness for several Paillier encryption.
//  This is for accelerating setup, but it is NOT SECURE,
//  however does not change the client nor server execution time
//  should be used only for benchmarking the audit part of the protocol
//=====================================================================
// #define VESPO_NOTSECURE 100

inline void random_precompute(Polynomial<bn_t>& reusedrands,
                              const paillier_pubkey_t& pub,
                              const paillier_prvkey_t& prv,
                              const int64_t psize) {
    if (reusedrands.degree()<(psize-1)) {
        std::clog << "*WARNING WARNING* Paillier randomness reused\n"
                  << "*WARNING WARNING* " << (reusedrands.degree()+1)
                  << " instead of " << psize
                  << "\n*WARNING WARNING* insecure, only benchmarking"
                  << std::endl;
    }
#pragma omp parallel for shared(reusedrands,prv,pub)
    for (int64_t i = 0; i <= reusedrands.degree(); ++i) {
		bn_rand_mod(reusedrands[i], prv->a);
		bn_mxp(reusedrands[i], prv->gn, reusedrands[i], pub.nsq);
    }
}

//=====================================================================

//=====================================================================
// VeSPo: linear algebra
//=====================================================================

inline void vector::init() {
    bn_null(this->first);
    bn_null(this->second);
    bn_new(this->first);
    bn_new(this->second);
}
inline vector::vector() { init(); }
inline vector::~vector() {
    bn_free(this->first);
    bn_free(this->second);
}
inline vector::vector(const bn_t& f, const bn_t& s) : vector() {
    bn_copy(this->first, f);
    bn_copy(this->second,s);
}
inline vector::vector(const vector& v) : vector() {
    this->copy(v);
}
inline void vector::randomize(const bn_t mod) {
    bn_rand_mod(this->first, mod);
    bn_rand_mod(this->second, mod);
}
inline void vector::modin(const bn_t& r) {
    bn_mod(this->first, this->first, r);
    bn_mod(this->second, this->second, r);
}
inline void vector::copy(const vector& v) {
    bn_copy(this->first, v.first);
    bn_copy(this->second,v.second);
}
inline void vector::addin(const vector& v) {
    bn_add(this->first, this->first, v.first);
    bn_add(this->second, this->second, v.second);
}
inline bool vector::areEqual(const vector& v) const {
    return ( (bn_cmp(this->first,v.first)== RLC_EQ) &&
             (bn_cmp(this->second,v.second)== RLC_EQ) );
}

inline std::ostream& operator<<(std::ostream& out, const vector& v) {
    return out << '<' << v.first << ',' << v.second << '>';
}

inline void matrix::init() {
    bn_null(a);bn_null(b);bn_null(c);bn_null(d);
    bn_new(a); bn_new(b); bn_new(c); bn_new(d);
}
inline matrix::~matrix() { bn_free(a); bn_free(b); bn_free(c); bn_free(d); }

inline matrix::matrix() { init(); }

inline matrix::matrix(const bn_t& Aa, const bn_t& Ab, const bn_t& Ac, const bn_t& Ad) : matrix() {
    bn_copy(a, Aa);
    bn_copy(b, Ab);
    bn_copy(c, Ac);
    bn_copy(d, Ad);
}
inline matrix::matrix(const matrix& A) : matrix(A.a,A.b,A.c,A.d) {}

inline void matrix::randomize(const bn_t& mod) {
    bn_rand_mod(a, mod);
    bn_rand_mod(b, mod);
    bn_rand_mod(c, mod);
    bn_rand_mod(d, mod);
}
inline void matrix::mulin(const bn_t& r) {
    bn_mul(a, a, r);
    bn_mul(b, b, r);
    bn_mul(c, c, r);
    bn_mul(d, d, r);
}
inline void matrix::modin(const bn_t& r) {
    bn_mod(a, a, r);
    bn_mod(b, b, r);
    bn_mod(c, c, r);
    bn_mod(d, d, r);
}
inline void matrix::copy(const matrix& A) {
    bn_copy(a, A.a);
    bn_copy(b, A.b);
    bn_copy(c, A.c);
    bn_copy(d, A.d);
}

std::ostream& operator<<(std::ostream& out, const matrix& m) {
    return out << "<<" << m.a << '|' << m.b << ">,<" << m.c << '|' << m.d << ">>";
}



inline void sqr_modcharp(bn_t& s0, bn_t& s1,
                         const bn_t& p0, const bn_t& p1,
                         const bn_t& t0, const bn_t& t1,
                         const bn_t& mod) {
        // (t0 + t1 X)^2 mod p0 + p1 X + X^2
        // --> t0^2 + 2t0t1 X + t1^2 ( -p0 -p1 X)
        // --> (t0^2 -p0t1^2) + (2t0t1 -p1t1^2) X

    bn_sqr(s0, t0);						// T0^2
    bn_mod(s0, s0, mod);
    bn_mul(s1, t0, t1);					// T0T1
    bn_mod(s1, s1, mod);
    bn_add(s1, s1, s1);					// 2T0T1

    bn_t tmp, smp; bn_new(tmp); bn_new(smp);

    bn_sqr(tmp, t1);					// T1^2
    bn_mod(tmp, tmp, mod);
    bn_mul(smp, tmp, p1);				// T1^2P1
    bn_mod(smp, smp, mod);
	bn_mul(tmp, tmp, p0);				// T1^2P0
	bn_mod(tmp, tmp, mod);
    bn_sub(s0, s0, tmp);				// T0^2-T1^2P0
    bn_mod(s0, s0, mod);
    bn_sub(s1, s1, smp);				// 2T0T1-T1^2P1
    bn_mod(s1, s1, mod);

    bn_free(tmp);bn_free(smp);
}


inline void mulinx_modcharp(bn_t& s0, bn_t& s1,
                            const bn_t& p0, const bn_t& p1,
                            const bn_t& mod) {
        // (s0 + s1 X) X mod p0 + p1 X + X^2
        // --> s0 X + s1 (-p0 -p1 X)
        // --> (-p0 s1) + (s0 - s1p1) X

	bn_t tmp; bn_t smp; bn_new(tmp); bn_new(smp);

	bn_mul(tmp, s1, p1);			// S1P1
	bn_mod(tmp, tmp, mod);
	bn_sub(smp, s0, tmp);			// S0-S1P1

	bn_mul(tmp, s1, p0);			// S1P0
	bn_neg(tmp, tmp);				// -S1P0

	bn_mod(s0, tmp, mod);			// Can now be overwritten
	bn_mod(s1, smp, mod);			// Can now be overwritten

	bn_free(tmp); bn_free(smp);
}

inline void TMMP(bn_t& s0, bn_t& s1,
          const bn_t& p0, const bn_t& p1,
          const int64_t d, const bn_t& mod) {
        // Computes s0+s1Z=Z^d modulo P(X)=p0+p1Z+Z^2
#ifdef DEBUG
    std::clog << "[2-MMP] " << d << " BEG" << std::endl;
#endif
    if (d<=1) {
        if (d==1) {
            bn_set_dig(s0,0);
            bn_set_dig(s1,1);
            return;
        }
        if (d==0) {
            bn_set_dig(s0,1);
            bn_set_dig(s1,0);
            return;
        }
        if (d<0) {
            std::cerr << "**** ERROR **** negative exponent not implemented\n";
            return;
        }
    }

    bn_t t0, t1; bn_null(t0); bn_null(t1); bn_new(t0); bn_new(t1);
    TMMP(t0, t1, p0, p1, d>>1, mod);			// Recursive call

    sqr_modcharp(s0, s1, p0, p1, t0, t1, mod);	// T^2		mod P
    bn_free(t0);bn_free(t1);

    if (d & 0x1) {
        mulinx_modcharp(s0, s1, p0, p1, mod);	// X T^2	mod P
    }

#ifdef DEBUG
    std::clog << "[2MMP] " << d << " END" << std::endl;
#endif
}


inline void moddet(bn_t& det,
            const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
            const bn_t& mod) {
        // Determinant of <<a|b>,<c|d>>
    bn_t tmp; bn_null(tmp); bn_new(tmp);

    bn_mul(det, a, d);		// ad
    bn_mul(tmp, b, c);		// bc
    bn_sub(det, det, tmp);	// ad-bc
    bn_mod(det, det, mod);

    bn_free(tmp);
}

inline void modcharpoly(bn_t& p0, bn_t& p1,
                 const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
                 const bn_t& mod) {
        // Characteristic polynomial of <<a|b>,<c|d>>
        // (ad-bc) - (a+d) X + X^2

    bn_t tmp; bn_null(tmp); bn_new(tmp);

    bn_add(p1, a, d);	// a+d
    bn_neg(p1, p1);
    bn_mod(p1, p1, mod);

    bn_free(tmp);

    moddet(p0, a, b, c, d, mod);
}

inline void matvectmod(bn_t& u0, bn_t& u1,
                const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
                const bn_t& v0, const bn_t& v1,
                const bn_t& mod) {
        // <u0,u1> <-- <<a|b>,<c|d>> . <v0,v1>
#ifdef DEBUG
    std::clog << matrix(a,b,c,d) << '.' << vector(v0,v1) << " - ";
#endif
        // < a v0 + b v1 , c v0 + d v1 >
    bn_t tmp; bn_null(tmp); bn_new(tmp);

    bn_mul(u0, a, v0);
    bn_mul(tmp, b, v1);
    bn_add(u0, u0, tmp);
    bn_mod(u0, u0, mod);

    bn_mul(u1, c, v0);
    bn_mul(tmp, d, v1);
    bn_add(u1, u1, tmp);
    bn_mod(u1, u1, mod);

    bn_free(tmp);
#ifdef DEBUG
    std::clog << vector(u0,u1) << " mod " << mod << ';' << std::endl;
#endif
}



inline void PMGS(vector& u, const matrix& M,
          const uint64_t k,
          const vector& v, const bn_t& mod) {
        // Projected matrix geometric sum
        // u <-- sum_{i=0}^k M^i . v
#ifdef DEBUG
    std::clog << "[PMGS] " << k << " BEG" << std::endl;
#endif
        // Characteristic polynomial
    bn_t p0, p1; bn_null(p0); bn_null(p1); bn_new(p0); bn_new(p1);
    modcharpoly(p0, p1, M.a, M.b, M.c, M.d, mod);

        // (A-I2)^{-1} --> 1/det < < d-1 | -b >, < -c | a-1 > >
    bn_t dm, am;
    bn_null(dm); bn_null(am);
    bn_new(dm); bn_new(am);
    bn_sub_dig(dm, M.d, 1); // d-1
    bn_sub_dig(am, M.a, 1); // a-1

        // Determinant( A-I2 ) --> am dm - b c
    bn_t idet; bn_null(idet); bn_new(idet);
    moddet(idet, am, M.b, M.c, dm, mod);
    bn_mod_inv(idet, idet, mod);


    bn_t t0, t1; bn_null(t0); bn_null(t1); bn_new(t0); bn_new(t1);
    bn_t w0, w1; bn_null(w0); bn_null(w1); bn_new(w0); bn_new(w1);
        // w <-- (idet)^{-1}(A-I2)^{-1} v
        // Direct computation: w = < (d-1)v0 - b v1 , (a-1)v1 -c v0 >
    bn_mul(t0, dm, v.first);	// (d-1)v0
    bn_mul(t1, M.b, v.second);	// b v1
    bn_sub(w0, t0, t1);			// (d-1)v0 - b v1
    bn_mod(w0, w0, mod);

    bn_mul(t0, am, v.second);	// (a-1)v1
    bn_mul(t1, M.c, v.first);	// c v0
    bn_sub(w1, t0, t1);			// (a-1)v1 - c v0
    bn_mod(w1, w1, mod);

        // Apply the inversed determinant
    bn_mul(w0, w0, idet);
    bn_mod(w0, w0, mod);
    bn_mul(w1, w1, idet);
    bn_mod(w1, w1, mod);

    bn_free(idet);

#ifdef DEBUG
    std::clog << "1/(" << M << "-IdentityMatrix(2))." << v << " - " << vector(w0,w1) << " mod " << mod << ';' << std::endl;
#endif


        // Polynomial Geometric sum
    bn_t f0, f1; bn_null(f0); bn_null(f1); bn_new(f0); bn_new(f1);
    TMMP(f0, f1, p0, p1, k+1, mod);

    bn_free(p0);bn_free(p1);

        // F(A) - I2 --> << a*f1+f0-1 | f1*b >, < f1*c | d*f1+f0-1 >>
    bn_sub_dig(dm, f0, 1);	// f0-1
    bn_mul(am, M.a, f1);	// a*f1
    bn_mod(am, am, mod);
    bn_add(am, am, dm);		// a*f1+f0-1

    bn_mul(t0, M.d, f1);	// d*f1
	bn_mod(t0, t0, mod);
    bn_add(dm, t0, dm);		// d*f1+f0-1

    bn_mul(t0, M.b, f1);	// f1*b
    bn_mod(t0, t0, mod);

    bn_mul(t1, M.c, f1);	// f1*c
    bn_mod(t1, t1, mod);

    bn_free(f0); bn_free(f1);

    matvectmod(u.first, u.second, am, t0, t1, dm, w0, w1, mod);

    bn_free(dm);bn_free(t0);bn_free(t1);bn_free(am);
    bn_free(t0);bn_free(t1);
    bn_free(w0);bn_free(w1);

#ifdef VESPO_CHECKERS
    {
        vector r, s, t;
        t.copy(v); s.copy(v);
        for(uint64_t i = 0; i<k; ++i) {
            matvectmod(r, M, t, mod);
            s.addin(r);
            t.copy(r);
        }
        s.modin(mod);
        if (u.areEqual(s)) {
            std::clog << "[ PMGS ] \033[1;32mOK\033[0m" << std::endl;
        } else {
            std::cerr << "[ PMGS ] \033[1;31m****** FAIL ******\033[0m" << std::endl;
        }
    }
#endif
#ifdef DEBUG
    std::clog << "[PMGS] " << k << " END" << std::endl;
#endif
}


inline void mat_geometric_sum(vector& c,
                       const bn_t& r, const matrix& A, const uint64_t k,
                       const vector& b,
                       const bn_t& mod) {

        // Projected matrix geometric sum
        // c <-- sum_{i=0}^k (rA)^i . b
    matrix M(A);
    M.mulin(r);
    M.modin(mod);

    PMGS(c, M, k, b, mod);
}

//=====================================================================
// VeSPo data structures
//       Client/Server
//=====================================================================

inline client_t::client_t(const int64_t degree) : d(degree) {
    g1_null(this->g1);g1_new(this->g1);
    g2_null(this->g2);g2_new(this->g2);
    gt_null(this->e_T);gt_new(this->e_T);
    paillier_newprvkey(this->prv);
    paillier_newpubkey(this->pub);
    bn_null(this->s);bn_new(this->s);
    gt_null(this->K1_bT);gt_new(this->K1_bT);
    gt_null(this->K2_bT);gt_new(this->K2_bT);
}

inline client_t::~client_t() {
    g1_free(this->g1);
    g2_free(this->g2);
    gt_free(this->e_T);
    paillier_freeprvkey(this->prv);
    paillier_freepubkey(this->pub);
    bn_free(this->s);
    gt_free(this->K1_bT);
    gt_free(this->K2_bT);
}

inline uint64_t client_t::private_size(uint64_t elements_size) const {
        // Groups, prime order and pubkey are supposed public data
        // The degree of the polynomial must be stored
    uint64_t size=0;
    size += sizeof(this->d)*8;	// degree
        // s,K_bT,valpha,vbeta,msigma + MerkleTreeRoot
    size += elements_size*(1+2+2+2+4 + 1);
    size += paillier_sizeprvkey( this->prv );
    return size;
}

inline uint64_t client_t::privpub_size(uint64_t elements_size) const {
        // adding keys
    return private_size(elements_size)
        + paillier_sizepubkey( this->pub );

}

    //=====================================================================
    // VeSPo: key generation
    //=====================================================================
inline void client_t::keygen(uint64_t pailliersize) {
        pc_get_ord( prv->a );

        cp_shpe_gen( pub.rlc, prv, bn_bits( prv->a ), pailliersize);
        bn_sqr(pub.nsq, pub.rlc->crt->n); // store n^2 too
}


inline int64_t setup_block_chunks(
    std::vector<Polynomial<paillier_ciphertext_t>>& W,
    const int64_t degree, const int64_t blocks, const bn_t& modulus) {
        // Cut W into 'blocks' chunks
        // last ones of size [ (degree+1)/Blocks ]
        // first ones of that size + 1
        // so that sum of sizes is 'degree+1'
    int64_t sW(degree+1);
    W.reserve(blocks);
    int64_t sdeg(sW/blocks);
    int64_t bigblocks = (sW-sdeg*blocks);
    const int64_t ldeg(sdeg); --sdeg;
    int64_t i=0; for(; i<bigblocks;++i) {
        W.emplace_back(ldeg, modulus);
    }
    for( ; i<blocks; ++i){
        W.emplace_back(sdeg, modulus);
    }
        // If no bigblocks, actually all are of same (big) size
    if (bigblocks == 0) bigblocks=blocks;
    return bigblocks;
}



inline server_t::server_t(const int64_t degree, const int64_t blocks,
                          const bn_t& modulus) :
        H1_b(degree-1,modulus), H2_b(degree-1,modulus),
        S(degree, modulus) {
        // For parallelism:
        //	Cut W into 'blocks' chunks
        //	last ones of size [ (degree+1)/Blocks ]
        //	first ones of that size + 1
        //	so that sum of sizes is 'degree+1'
    this->bigblocks=setup_block_chunks(this->W, degree, blocks, modulus);
}



//=====================================================================
// VeSPo: checkers implementations
//=====================================================================

inline void check_dp(g1_t& dp, const Polynomial<g1_t>& gS, const Polynomial<bn_t>& P) {
        // dp <-- \oplus [ p_i ] gS_i
    g1_t tmp; g1_new(tmp);
    g1_mul(dp, gS[0], P[0]);
    for(int64_t i=1; i <= P.degree(); ++i) {
        g1_mul(tmp, gS[i], P[i]);
        g1_add(dp, dp, tmp);
    }
    g1_free(tmp);
}


inline bool check_orders(const client_t& client, const bn_t& pairing_r) {
        // pairing sanity checks
#ifdef DEBUG
    printer(std::clog << "g1 : ", client.g1) << std::endl;
    printer(std::clog << "g2 : ", client.g2) << std::endl;
#endif

    g1_t o1; g1_new(o1);
    g1_mul_gen(o1, pairing_r);

    g2_t o2; g2_new(o2);
    g2_mul_gen(o2, pairing_r);

    gt_t ot; gt_new(ot);
    gt_exp_gen(ot, pairing_r);

    gt_t generatorT;
    gt_new(generatorT);
    pc_map(generatorT, client.g1, client.g2); // g_t

    const bool pass( (! g1_is_infty(client.g1)) &&
                     (! g2_is_infty(client.g2)) &&
                     (! gt_is_unity(client.e_T)) &&
                     g1_is_infty(o1) && g2_is_infty(o2) && gt_is_unity(ot) &&
                     (gt_cmp(client.e_T, generatorT) == RLC_EQ)
                     );

    if(pass) {
        std::clog << "[GOrders ] \033[1;32mOK\033[0m" << std::endl;

    } else {
        std::clog << "[GOrders ] \033[1;31m****** FAIL ******\033[0m" << std::endl;
#ifdef DEBUG
        printf("o1 (should be zero): "); g1_print(o1); printf("\n");
        printf("o2 (should be zero): "); g2_print(o2); printf("\n");
        printf("gt : "); gt_print(generatorT); printf("\n");
        printf("et : "); gt_print(client.e_T); printf("\n");
#endif
    }

    g1_free(o1); g2_free(o2);

    return pass;
}


inline bool check_pubkey(const client_t& client, const server_t& server,
                  const Polynomial<bn_t>& P_b,
                  const bn_t& eval, const g2_t& K_b) {
        // K_b sanity checks
	g1_t test,chet;

	g1_new(test);
	check_dp(test, server.S, P_b);

	g1_new(chet);
	g1_mul_gen(chet, eval);		// g1^{P_b(s)}

	if(g1_cmp(test, chet) == RLC_EQ) {
	    std::clog << "[ PowerP ] \033[1;32mOK\033[0m" << std::endl;
	} else {
#   ifdef DEBUG
        {
            printf("g1 test : "); g1_print(test); printf("\n");
            printf("g1 chet : "); g1_print(chet); printf("\n");
        }
#   endif
	    std::cerr << "[ PowerP ] \033[1;31m****** FAIL ******\033[0m" << std::endl;
	}

	gt_t left, right;
	gt_new(left);
	gt_new(right);
	pc_map(left, client.g1, K_b);	// e(g1;K_b)
	pc_map(right, test, client.g2);	// e(dp;g2)

#   ifdef DEBUG
    {
        printf("left  : "); gt_print(left); printf("\n");
        printf("right : "); gt_print(right);printf("\n");
    }

#   endif

    const bool pass( gt_cmp(left,right) == RLC_EQ );
    if(pass) {
        std::clog << "[ KeyPub  ] \033[1;32mOK\033[0m" << std::endl;
    } else {
        std::cerr << "[ KeyPub ] \033[1;31m****** FAIL ******\033[0m" << std::endl;
    }
    g1_free(test);
    gt_free(left); gt_free(right);
    g1_free(chet);

    return pass;
}

inline bool check_g2_omxp(const g2_t& T, const Polynomial<g2_t>& S,
                   const bn_t& r, const int64_t deg, const bn_t& mod) {
    bool pass(true);
    g2_t U;
    g2_copy(U, S[0]);
    for (int64_t i = 1; i <= deg; ++i) {
        g2_mul(U, U, r);		// t^r
        g2_add(U, U, S[i]);		// S_{i-1} t^r
    }
    if ( (g2_cmp(U,T) != RLC_EQ)) {
        printer(printer(std::cerr << "check_g2_omxp:", U) << " != ", T) << std::endl;
        pass = false;
    }
    return pass;
}


inline bool check_g1_hmxp(const Polynomial<g1_t>& T, const Polynomial<g1_t>& S,
                   const bn_t& r, const int64_t from, const int64_t length,
                   const bn_t& mod) {
    const int64_t deg(from+length-1);
    bool pass(true);
    g1_t U;
    g1_copy(U, S[0]);
    if (g1_cmp(U, T[0]) != RLC_EQ) {
        printer(printer(std::cerr << "check_g1_hmxp: U[0] != T[0]: ", U) << ' ', T[0]) << std::endl;
        pass = false;
    }
    for (int64_t i = 1; i <= deg; ++i) {
        g1_mul(U, U, r);		// t^r
        g1_add(U, U, S[i]);		// S_{i-1} t^r
        if ( (i>= from) && (g1_cmp(U, T[i]) != RLC_EQ)) {
            printer(printer(std::cerr << "check_g1_hmxp: U[" << i << "] != T[" << i << "]: ", U) << ' ', T[i]) << std::endl;
            pass = false;
        }
    }
    return pass;
}


inline void check_g1_msl(g1_t& r, const g1_t* P, const bn_t* K, int N) {
    g1_t q; g1_null(q); g1_new(q);
    g1_set_infty(r);
    for (int j = 0; j < N; ++j) {
        g1_mul(q, P[j], K[j]);
        g1_add(r, r, q);
        g1_mul_sim_lot(q,P,K,j+1);
        if ( ! (g1_cmp(q,r) == RLC_EQ) ) {
            std::cerr << "\033[1;31m****** FAIL: g1_mul_sim_lot " << (j+1) << " ******\033[0m" << std::endl;
        }
    }
    g1_free(q);
}

inline void check_g1_mul_sim_lot(g1_t& RES, const g1_t* P, const bn_t* K, int N) {
#ifdef VESPO_CHECKERS
    Chrono c_step; c_step.start();
#endif
    g1_mul_sim_lot(RES,P,K,N);

#ifdef VESPO_CHECKERS
    g1_t r; g1_null(r); g1_new(r);
    check_g1_msl(r, P, K, N);
    if ( ! (g1_cmp(RES,r) == RLC_EQ) ) {
        std::cerr << "\033[1;31m****** FAIL: g1_mul_sim_lot " << N << " last ******\033[0m" << std::endl;
    }
    g1_free(r);
#endif
#ifdef VESPO_CHECKERS
    std::clog << "    Group 1 MulSim : " << c_step.stop() << " (" << N << " operations)" << std::endl;
#endif
}

inline bool check_ciph_horner(gt_t sxi_b, const int64_t deg, const server_t& server,
                       const Polynomial<g2_t>& H_b, const bn_t& rmod) {
        // Server computation: check xi
    gt_t xi_b, tmpT;
    g1_t tmpG1;
    gt_new(tmpT);
    gt_new(xi_b);
    g1_new(tmpG1);
    gt_set_unity(xi_b);

    bn_t zero; vespo_init_set_ui(zero,0);
    g1_rand(tmpG1); g1_mul(tmpG1,tmpG1,zero);

    for (int64_t i = 0; i < deg; ++i) {
        g1_mul(tmpG1,tmpG1, rmod);                     // t^r
        g1_add(tmpG1,server.S[i], tmpG1);      // S_{i-1} t^r
#ifdef DEBUG
        printf("H[%lu]: ", i); g2_print(H_b[i]); printf("\n");
        printf("T[%lu]: ", i); g1_print(tmpG1); printf("\n");
#endif
        pc_map(tmpT, tmpG1, H_b[i]);   // e(H_i;t)
        gt_mul(xi_b, xi_b, tmpT);                      // xi e(H_i;t)
    }

    const bool pass(gt_cmp(xi_b, sxi_b) == RLC_EQ);

    if (pass) {
        std::clog << "[CiHorner] \033[1;32mOK\033[0m" << std::endl;
    } else {
        std::cerr << "[CiHorner] \033[1;31m****** FAIL ******\033[0m" << std::endl;
#ifdef DEBUG
        printf("xi_b : "); gt_print(xi_b); printf("\n");
        printf("sxi_b: "); gt_print(sxi_b); printf("\n");
#endif
    }

    gt_free(xi_b);
    g1_free(tmpG1);
    gt_free(tmpT);

    return pass;
}


//=====================================================================
// VeSPo: masking
//=====================================================================

inline void vhide_poly(Polynomial<bn_t> & P1, Polynomial<bn_t>& P2,
                       const vector & valpha, const vector& vbeta,
                       const matrix& msigma, const Polynomial<bn_t> & P) {
        // <P1,P2> <-- sum X^i ( p_i valpha + msigma^i vbeta)
    std::clog << "[VHide P] BEG: d°" << P.degree() << std::endl;
#ifdef VESPO_SUB_TIMINGS
    Chrono c_vhide; c_vhide.start();
#endif
#ifdef DEBUG
    std::clog << P << std::endl;
#endif
    vector sigmapow, tmp; sigmapow.copy(vbeta);
        // p_0 alpha + beta
    bn_mul(P1[0], valpha.first, P[0]);
    bn_mul(P2[0], valpha.second, P[0]);
    bn_add(P1[0], P1[0], vbeta.first);
    bn_add(P2[0], P2[0], vbeta.second);
    bn_mod(P1[0], P1[0], P.modulus());
    bn_mod(P2[0], P2[0], P.modulus());
    for(int64_t i = 1; i <= P.degree(); ++i) {
        bn_mul(P1[i], valpha.first, P[i]);	// alpha p_i
        bn_mul(P2[i], valpha.second, P[i]);	// alpha p_i
        matvectmod(tmp, msigma, sigmapow, P.modulus());	// Sigma^i beta
        sigmapow.copy(tmp);
        bn_add(P1[i], P1[i], sigmapow.first);
        bn_mod(P1[i], P1[i], P.modulus()); // alpha p_i + beta sigma^i
        bn_add(P2[i], P2[i], sigmapow.second);
        bn_mod(P2[i], P2[i], P.modulus()); // alpha p_i + beta sigma^i
    }

#ifdef VESPO_SUB_TIMINGS
    std::clog << "  SETUP VHide polynom. : " << c_vhide.stop() << " (" << (P.degree()+1) << " operations)" << std::endl;
#endif

#ifdef VESPO_CHECKERS
    bn_t ttmp, etmp; bn_new(ttmp); bn_new(etmp);
    vector check, htmp;
    bn_rand_mod(etmp, P.modulus()); bn_mod(etmp, etmp, P.modulus());
    P.eval(ttmp, etmp);
    mat_geometric_sum(check, etmp, msigma, P.degree(), vbeta, P.modulus());
    vespo_addinmul(check.first, valpha.first, ttmp);
    vespo_addinmul(check.second, valpha.second, ttmp);
    check.modin(P.modulus());

    P1.eval(htmp.first, etmp);
    P2.eval(htmp.second, etmp);

    if( htmp.areEqual(check) ) {
        std::clog << "[ VHidden ] \033[1;32mOK\033[0m" << std::endl;
    } else {
        std::cerr << "[ VHidden ] \033[1;31m****** FAIL ******\033[0m" << std::endl;
    }

    bn_free(ttmp); bn_free(etmp);
#endif

    std::clog << "[VHide P] END" << std::endl;
}

//=====================================================================
// VeSPo: Paillier
//=====================================================================

inline void geo_progression(Polynomial<bn_t>& pows_r,
                            const int64_t bigsdegree, const bn_t& r1,
                            const int64_t totaldegree, const bn_t& r2,
                            const bn_t& mod, double& time_r) {
        // Computes 1, r1, ..., r1^bigsdegree,
        // r1^bigsdegree r2, ..., r1^bigsdegree r2^(totaldegree-bigsdegree)
#ifdef VESPO_CHECKERS
    std::clog << "[PGeoProg] BEG: d°" << totaldegree << " (1 then " << bigsdegree << " then  " << (totaldegree-bigsdegree) << ')' << std::endl;
#endif
#ifdef VESPO_TIMINGS
    Chrono c_step; c_step.start();
#endif
    bn_set_dig(pows_r[0],1);
    int64_t i = 1; for (; i <= bigsdegree; ++i) {
        bn_mul(pows_r[i], pows_r[i-1], r1);
        bn_mod(pows_r[i], pows_r[i], mod);
    }
    for (; i <= totaldegree; ++i) {
        bn_mul(pows_r[i], pows_r[i-1], r2);
        bn_mod(pows_r[i], pows_r[i], mod);
    }

#ifdef VESPO_TIMINGS
    time_r = c_step.stop();
    std::clog << "  SERVER geom. prog: " << time_r << " (" << (totaldegree+1) << " Group mul operations)" << std::endl;
#endif
#ifdef VESPO_CHECKERS
    std::clog << "[PGeoProg] END" << std::endl;
#endif
}


inline void paillier_hom_dp(paillier_ciphertext_t& eval,
                     const paillier_pubkey_t& kpub,
                     const Polynomial<paillier_ciphertext_t>& Wi,
                     const Polynomial<bn_t>& pows_r,
                     const int64_t degree) {
        // eval <-- prod Wi[i]^pows_r[i]
#ifdef VESPO_SUB_TIMINGS
    Chrono c_step; c_step.start();
#endif
#ifdef VESPO_CHECKERS
    std::clog << "[PailCdot] BEG: d°" << degree << std::endl;
#endif

    if (degree>=0)
        bn_mxp_sim_lot(eval.c, reinterpret_cast<const bn_t*>( &(Wi[0])),
                       reinterpret_cast<const bn_t*>(&(pows_r[0])),
                       kpub.nsq, degree+1);
    else
        bn_set_dig(eval.c, 1);

#ifdef VESPO_TIMINGS
    double time_r = c_step.stop();
    std::clog << "    Paillier mxpsim: " << time_r << " (" << (degree+1) << " homomorphic ops)\n";
#endif
#ifdef VESPO_CHECKERS
    std::clog << "[PailCdot] END" << std::endl;
#endif
}

inline void encrypt_poly(std::vector<Polynomial<paillier_ciphertext_t>>& W,
                  const Polynomial<bn_t>& P,
                  const paillier_pubkey_t& pub, const paillier_prvkey_t& prv) {
        // W <-- E(P), and W is cut into chunks
    std::clog << "[PAILLIER ENC] BEG: d°" << P.degree() << std::endl;
#ifdef VESPO_SUB_TIMINGS
    Chrono c_paillier; c_paillier.start();
#endif
    PAILLIER_PRECOMPU(pub,prv,P.degree()+1);

    int64_t i = 0;

    for(size_t j=0; j<W.size(); ++j) {
        for(int64_t k=0; k<=W[j].degree(); ++k, ++i) {
            PAILLIER_ENCIPHER(W[j][k], P[i], pub, prv, i);
        }
    }

#ifdef VESPO_SUB_TIMINGS
    std::clog << "  SETUP Paillier Cipher: " << c_paillier.stop() << " (" << (P.degree()+1) << " operations)" << std::endl;
#endif
    std::clog << "[PAILLIER ENC] END" << std::endl;
}

//=====================================================================
// VeSPo: in the exponents
//=====================================================================

inline Polynomial<g2_t>& powered_poly_gen2(Polynomial<g2_t>& H,
                                    const Polynomial<bn_t>& P) {
        // hi[i] <-- g2^{p_(i+1)}
    std::clog << "[Powered] BEG: d°" << P.degree() << std::endl;
#ifdef VESPO_SUB_TIMINGS
    Chrono c_powpol; c_powpol.start();
#endif

#pragma omp parallel for shared(P,H)
    for (int64_t i = 0; i < P.degree(); ++i) {
        g2_mul_gen(H[i], P[i+1]);	// g2^{p_(i+1)}
    }

#ifdef VESPO_SUB_TIMINGS
    std::clog << "  SETUP Enpowered poly.: " << c_powpol.stop() << " (" << (P.degree()+1) << " operations)" << std::endl;
#endif
    std::clog << "[Powered] END" << std::endl;
    return H;
}

inline Polynomial<g1_t>& power_progression_gen(Polynomial<g1_t>& S,
                                               const bn_t& s,
                                               const int64_t d,
                                               const bn_t& mod) {
        // S[i] <-- g1^{s^i}
    std::clog << "[PowProg] BEG: d°" << d << std::endl;
#ifdef VESPO_SUB_TIMINGS
    Chrono c_powprog; c_powprog.start();
#endif
    double time_s(0.0);
    Polynomial<bn_t> pows_s(d,mod);
        // 1, s, s^2, ..., s^d
    geo_progression(pows_s, 0, s, d, s, mod, time_s);

    g1_get_gen(S[0]);					// Generator in G1
#pragma omp parallel for shared(pows_s,S)
    for (int64_t i = 1; i <= d; ++i) {
        g1_mul_gen(S[i], pows_s[i]);	// g1^{s^i}
    }

#ifdef VESPO_SUB_TIMINGS
    std::clog << "  G1 Power progress.: " << c_powprog.stop() << " (" << d << " operations)" << std::endl;
#endif
    std::clog << "[PowProg] END" << std::endl;
    return S;
}



//========================================================
// VeSPo: Init protocol
//========================================================
inline void setup(client_t& client, server_t& server,
           uint64_t pailliersize,
           const Polynomial<bn_t>& P, double& t_horner) {

    std::clog << "[ Setup  ] BEG: d°" << P.degree() << std::endl;

    vector veval;
    bn_t group_mod; bn_new(group_mod); pc_get_ord(group_mod);

        // Random valpha, vbeta
    client.valpha.randomize(group_mod);
    client.vbeta.randomize(group_mod);

        // Find s and msigma, s.t. (s Sigma - I_2) is invertible
    bn_t det; bn_new(det);
    do {
        bn_rand_mod(client.s, group_mod);
        client.msigma.randomize(group_mod);
        matrix T(client.msigma); T.mulin(client.s);
        bn_sub_dig(T.a,T.a,1);
        bn_sub_dig(T.d,T.d,1);
        T.modin(group_mod);      // T <-- s Sigma - I_2
        moddet(det, T, group_mod);
        std::clog << "[Setup det(sSigma-I2): " << det << ']' << std::endl;
    } while(bn_is_zero(det));

    std::clog << "[Setup group secret " << client.s << " (mod " << group_mod << ")]" << std::endl;

        // Generate Paillier's keys
    client.keygen(pailliersize);

#ifdef VESPO_CHECKERS
    bn_sqr(det,group_mod);
    bn_mul_dig(det,det,P.degree()+1);
    if ( bn_cmp(client.pub.rlc->crt->n, det) != RLC_GT) {
        std::cerr << "\033[1;31m****** FAIL: Paillier too small ******\033[0m" << std::endl;
        bn_print(det);
        std::cerr << std::endl;
        bn_print(client.pub.rlc->crt->n);
        exit(-1);
    }
#endif
    bn_free(det);

    std::clog << "[Setup Paillier mod "<< client.pub.rlc->crt->n << ']' << std::endl;

    server.pub = client.pub;
        // Pairing Generators
    g1_get_gen(client.g1);		//     g1_rand(client.g1);	// Generator in G1
    g2_get_gen(client.g2);		//     g2_rand(client.g2);	// Generator in G2
    gt_get_gen(client.e_T);

    Polynomial<bn_t>
        P1_b(P.degree(), P.modulus()),
        P2_b(P.degree(), P.modulus());

        // parallel store S on server
    power_progression_gen(server.S, client.s, P.degree(), P.modulus());
        // parallel store W on server
    encrypt_poly(server.W, P, client.pub, client.prv);
        // parallel hide the polynomial
    vhide_poly(P1_b, P2_b, client.valpha, client.vbeta, client.msigma, P);
        // parallel store H_b on server
    powered_poly_gen2(server.H1_b, P1_b);
    powered_poly_gen2(server.H2_b, P2_b);

#ifdef VESPO_TIMINGS
    Chrono c_horner; c_horner.start();
#endif

    P1_b.eval(veval.first, client.s);	// P1_b(s)
    P2_b.eval(veval.second, client.s);	// P2_b(s)

#ifdef VESPO_TIMINGS
    t_horner = c_horner.stop() / 2.;	// Two "pure" Horner evaluations
#endif


#ifdef VESPO_CHECKERS
    check_orders(client, group_mod);
#endif

        // Compute and store public key K on client
    g2_t K1_b, K2_b;
    g2_new(K1_b);g2_new(K2_b);

    g2_mul_gen(K1_b, veval.first);	// g2^{P1_b(s)}
    pc_map(client.K1_bT, client.g1, K1_b);
    g2_mul_gen(K2_b, veval.second);	// g2^{P2_b(s)}
    pc_map(client.K2_bT, client.g1, K2_b);


#ifdef VESPO_CHECKERS
    check_pubkey(client, server, P1_b, veval.first, K1_b);
    check_pubkey(client, server, P2_b, veval.second, K2_b);
#endif

    g2_free(K1_b); g2_free(K2_b);
    bn_free(group_mod);
    std::clog << "[ Setup  ] END" << std::endl;
}

//========================================================
// VeSPo: Update
//========================================================
inline bool update(client_t& client, server_t& server,
            Polynomial<bn_t>& P, const bn_t& delta, const size_t index,
            double& time_u) {
#ifdef VESPO_TIMINGS
    std::clog << "[Update " << index  << ':' << delta << "] BEG" << std::endl;
#endif
    Chrono t_upd; t_upd.start();

        // Polynomial update
    bn_add(P[index], P[index], delta);
    bn_mod(P[index], P[index], P.modulus());

        // Client updates
        // Client update (1): E(delta)
    paillier_ciphertext_t edelta;
    cp_shpe_enc_prv(edelta.c,delta,client.prv);

#ifdef VESPO_CHECKERS
    bn_t Pr; bn_new(Pr);
    cp_shpe_dec(Pr,edelta.c,client.prv);
    if(bn_cmp(Pr,delta) == RLC_EQ) {
        std::clog << "[UPD Poly] \033[1;32mOK\033[0m " << delta << std::endl;
    } else {
        std::cerr << "[UPD Poly] \033[1;31m****** FAIL ******\033[0m" << std::endl;
    }
#endif

        // Client update (2): W
    bn_t da1, da2; bn_null(da1); bn_new(da1); bn_null(da2); bn_new(da2);
    g2_t Hda1, Hda2; g2_new(Hda1); g2_new(Hda2);

    bn_mul(da1, delta, client.valpha.first);
    bn_mod(da1, da1, P.modulus());
    g2_mul_gen(Hda1, da1);	// g2^{da1}
    bn_mul(da2, delta, client.valpha.second);
    bn_mod(da2, da2, P.modulus());
    g2_mul_gen(Hda2, da2);	// g2^{da1}

        // Client update (3): K
    bn_t si; bn_new(si);
    bn_mxp_dig(si, client.s, index, P.modulus()); // s^i mod p

    g2_t Deltaj; g2_null(Deltaj); g2_new(Deltaj);
    gt_t DeltagT; gt_new(DeltagT);

    g2_mul(Deltaj, Hda1, si);
    pc_map(DeltagT, client.g1, Deltaj);
    gt_mul(client.K1_bT, client.K1_bT, DeltagT);

    g2_mul(Deltaj, Hda2, si);
    pc_map(DeltagT, client.g1, Deltaj);
    gt_mul(client.K2_bT, client.K2_bT, DeltagT);

        // Server updates
        // Server update (1): W[index] <- W[index] * delta
    size_t jloc(0);
    int64_t kloc(index);
    for(; jloc<server.W.size() && kloc>server.W[jloc].degree(); ++jloc) {
        kloc -= (server.W[jloc].degree()+1);
    }
    if (0<jloc && jloc >= server.W.size()) --jloc;

#ifdef VESPO_CHECKERS
    cp_shpe_dec(Pr, server.W[jloc][kloc].c, client.prv);
#endif

    paillier_mul(server.W[jloc][kloc], server.W[jloc][kloc], edelta, server.pub);
        // Server update (2): H1[index] <- H[index] * Delta
    if (index>0) {
        g2_add(server.H1_b[index-1], server.H1_b[index-1], Hda1);
        g2_add(server.H2_b[index-1], server.H2_b[index-1], Hda2);
    }

    bn_free(da1); bn_free(da2); bn_free(si);
    g2_free(Hda1); g2_free(Hda2); g2_free(Deltaj);
    gt_free(DeltagT)

        time_u = t_upd.stop();

#ifdef VESPO_CHECKERS
    bn_t Pn; bn_new(Pn);
    cp_shpe_dec(Pn,server.W[jloc][kloc].c,client.prv);
    bn_add(Pr, Pr, delta);
    if(bn_cmp(Pr,Pn) == RLC_EQ) {
        std::clog << "[UPD (" << jloc << ',' << kloc << ')' << " Hom.] \033[1;32mOK\033[0m " << Pn << std::endl;
    } else {
        std::cerr << "[UPD (" << jloc << ',' << kloc << ')' << " Hom.] \033[1;31m****** FAIL ******\033[0m" << std::endl;
    }
    bn_free(Pr); bn_free(Pn);
#endif
#ifdef VESPO_TIMINGS
    std::clog << "[Update " << index << "] END" << std::endl;
#endif
    return true;
}


//========================================================
// VeSPo: homomorphic operations
//========================================================

// Threshold for small block
#define VESPO_G1HM_SMALL 32

inline void g1_mul_sim_lot_par(g1_t& RES, const g1_t* P, const bn_t* K, int N,
                        const bn_t& mod, const int nbtasks) {
#ifdef VESPO_SUB_TIMINGS
    Chrono c_step; c_step.start();
#endif

    int64_t sizeloop(N);
    if (nbtasks <= 1 || N < VESPO_G1HM_SMALL || N < nbtasks) {
        g1_mul_sim_lot(RES,P,K,N);
    } else {
        int64_t taskfactor(nbtasks<<1);
        do {
            sizeloop = std::ceil((double)N/(double(taskfactor)));
            taskfactor <<= 1;
        } while( sizeloop >= VESPO_RELIC_LIMIT_MAX_ALLOC );

        const int64_t nbblocks( (int64_t)(std::ceil((double)N/(double(sizeloop)) )));
#ifdef VESPO_CHECKERS
        std::clog << "[G1MuSiLp] taskfactor : " << taskfactor << std::endl;
        std::clog << "[G1MuSiLp] sizeloop   : " << sizeloop << std::endl;
        std::clog << "[G1MuSiLp] nbblocks   : " << nbblocks << std::endl;
#endif
        Polynomial<g1_t> RR(nbblocks-1,mod);
#pragma omp parallel for shared(RR,P,K)
        for(int64_t i=0; i<nbblocks; ++i) {
            const int64_t ipp(i*sizeloop);
            g1_mul_sim_lot(RR[i],&(P[ipp]),&(K[ipp]),
                           std::min(static_cast<int64_t>(sizeloop), N-ipp));
        }

        g1_set_infty(RES);
        for(int64_t i=0; i<nbblocks; ++i)
            g1_add(RES, RES, RR[i]);

    }
#ifdef VESPO_CHECKERS
    g1_t r; g1_null(r); g1_new(r);
    check_g1_msl(r, P, K, N);
    if ( ! (g1_cmp(RES,r) == RLC_EQ) ) {
        std::cerr << "\033[1;31m****** FAIL: g1_mul_sim_lot " << N << " last ******\033[0m" << std::endl;
    }
    g1_free(r);
#endif
#ifdef VESPO_SUB_TIMINGS
    std::clog << "    Group 1 MulSimP: " << c_step.stop() << " (" << N << " operations, by " << sizeloop << ")" << std::endl;
#endif
}

inline void g1_horner_mxp_seq(Polynomial<g1_t>& U, const Polynomial<g1_t>& S,
                       const bn_t& r, const int64_t from, const int64_t length) {
    for (int64_t i = 1; i < length; ++i) {
        g1_mul(U[from+i], U[from+i-1], r);			// t^r
        g1_add(U[from+i], U[from+i], S[from+i]);	// S_{i-1} t^r
    }
}

inline void g1_horner_mxp_iter(Polynomial<g1_t>& U, const Polynomial<g1_t>& S,
                        const bn_t& r, const int64_t from, const int64_t length,
                        const Polynomial<bn_t>& revpow, const int64_t step,
                        const bn_t& mod, const int64_t nbtasks) {
        // Server xi computation
        //   first step : Horner-like prefix computations of S_j^{x_k}
#ifdef VESPO_CHECKERS
    std::clog << "[G1HorIte] BEG: len " << length << ", tasks=" << nbtasks
#	if defined(_OPENMP)
              << " on " << omp_get_thread_num()
#	endif
              << std::endl;
#endif
#ifdef VESPO_SUB_TIMINGS
    Chrono c_cho; c_cho.start();
#endif
    g1_copy(U[0], S[0]);

    if (nbtasks<=1 || step<=1 || length<=VESPO_G1HM_SMALL || length<=nbtasks) {
        g1_horner_mxp_seq(U,S,r,from,length);
    } else {
        const int64_t smallloops(step*nbtasks-length);
        const int64_t largeloops(nbtasks-smallloops);
        const int64_t sstep(step-1);
        const int64_t lstep(step);

        int64_t smallprecomp(smallloops);
        int64_t largeprecomp(largeloops);
		// one precomputation less than tasks
        if (smallloops) {
            --smallprecomp;
        } else {
            --largeprecomp;
        }

        g1_t tmp; g1_null(tmp); g1_new(tmp);
        for(int64_t i=0; i<largeprecomp; ++i) {
            const int64_t bipp(from+i*lstep);
            const int64_t eipp(bipp+lstep);
            g1_mul_sim_lot_par(tmp, &(S[bipp+1]), &(revpow[1]), lstep, mod, nbtasks);
            g1_mul(U[eipp], U[bipp], revpow[0]);
            g1_add(U[eipp], U[eipp], tmp);
        }

        for(int64_t i=0; i<smallprecomp; ++i) {
            const int64_t bipp(from+largeprecomp*lstep+i*sstep);
            const int64_t eipp(bipp+sstep);
            g1_mul_sim_lot_par(tmp, &(S[bipp+1]), &(revpow[2]), sstep, mod, nbtasks);
            g1_mul(U[eipp], U[bipp], revpow[1]);
            g1_add(U[eipp], U[eipp], tmp);

        }


#pragma omp parallel for shared(U,S,r)
        for(int64_t i=0; i<nbtasks; ++i) {
            if (i<largeloops) {
                const int64_t bipp(from+i*step);
                const int64_t mstep(lstep);
                g1_horner_mxp_seq(U, S, r, bipp, mstep);
            } else {
                const int64_t bipp(from+largeloops*lstep+(i-largeloops)*sstep);
                const int64_t mstep(sstep);
                g1_horner_mxp_seq(U, S, r, bipp, mstep);
            }
        }
    }
#ifdef VESPO_SUB_TIMINGS
    std::clog << "    Group 1 Horner : " << c_cho.stop() << " (" << length << " operations)" << std::endl;
#endif
#ifdef VESPO_CHECKERS
    if (check_g1_hmxp(U, S, r, from, length, mod))
        std::clog << "OK check_g1_hmxp" << std::endl;
    else
        std::cerr << "********* ERROR ********** check_g1_hmxp" << std::endl;
    std::clog << "[G1Horte] END" << std::endl;
#endif
}


#define VESPO_PC_MAP_SIM pc_map_sim

inline void pairing_sim(gt_t& xi_b, const g1_t* Ti, const g2_t* Hi, const int64_t deg) {
#ifdef DEBUG
    std::clog << "[PairiSim] BEG: d°" << deg << std::endl;
#endif
    gt_set_unity(xi_b);
    if (deg <= VESPO_RELIC_LIMIT_MAX_ALLOC) {
        VESPO_PC_MAP_SIM(xi_b, Ti, Hi, deg);
    } else {
        const size_t endlcm = ( (deg-1)/VESPO_RELIC_LIMIT_MAX_ALLOC)*VESPO_RELIC_LIMIT_MAX_ALLOC;
        gt_t txi_b; gt_new(txi_b);

        size_t i=0; for( ; i<endlcm; i+=VESPO_RELIC_LIMIT_MAX_ALLOC) {
            VESPO_PC_MAP_SIM(txi_b, &(Ti[i]), &(Hi[i]), VESPO_RELIC_LIMIT_MAX_ALLOC);
            gt_mul(xi_b, xi_b, txi_b);
        }
        VESPO_PC_MAP_SIM(txi_b, &(Ti[i]), &(Hi[i]), deg-endlcm);
        gt_mul(xi_b, xi_b, txi_b);

        gt_free(txi_b);
    }
#ifdef DEBUG
    std::clog << "[PairiSim] END" << std::endl;
#endif
}

inline void pairing_double(gt_t& xi_b, const int64_t length,
                    const int64_t nbblocks, const bn_t& mod,
                    const Polynomial<g1_t>& Ti, const Polynomial<g2_t>& Hi) {
        // Server xi computation
        //   second step: applying the pairing maps as a dotproduct
#ifdef VESPO_CHECKERS
    std::clog << "[PairDblH] BEG: len. " << length << std::endl;
#endif
#ifdef VESPO_SUB_TIMINGS
    Chrono c_cho; c_cho.start();
#endif
    if (length) {
        int64_t dnbblocks(std::min(nbblocks, length));
        const int64_t sizeloop( (int64_t)(std::ceil((double)length/(double(dnbblocks)))));
        dnbblocks = (int64_t)(std::ceil((double)length/(double(sizeloop))));

        Polynomial<gt_t> Txi(dnbblocks-1,mod);
#pragma omp parallel for shared(Txi,Ti,Hi)
        for(int64_t i=0; i<dnbblocks; ++i) {
            const int64_t ipp(i*sizeloop);
            pairing_sim(Txi[i], &(Ti[ipp]), &(Hi[ipp]),
                        std::min(static_cast<int64_t>(sizeloop), length-ipp));
        }

        gt_copy(xi_b, Txi[0]);
        for(int64_t i=1; i<dnbblocks; ++i)
            gt_mul(xi_b, xi_b, Txi[i]);
    }
//     pairing_sim(xi_b, &(Ti[0]), &(Hi[0]), length);

#ifdef VESPO_SUB_TIMINGS
    std::clog << "    Pairings prods.: " << c_cho.stop() << " (" << length << " operations)" << std::endl;
#endif
#ifdef VESPO_CHECKERS
    std::clog << "[PairDblH] END" << std::endl;
#endif
}


//========================================================
// VeSPo: Audit
//========================================================
inline bool eval(paillier_plaintext_t& z, const client_t& client,
          const server_t& server,
          const bn_t& r, const int64_t nb_tasks,
          double& time_c, double& time_s
#ifdef VESPO_TIMINGS
          , double& time_e, double& time_p,
          double& time_ce, double& time_cg, double& time_cp
#endif
          ) {

#ifdef VESPO_TIMINGS
    std::clog << "[Audit " << r << "] BEG" << std::endl;
#endif
    const int64_t nbblocks(server.W.size());
    const int64_t dblocks(nbblocks-1);
    double time_r(0.0);

    bn_t group_mod; bn_new(group_mod); pc_get_ord(group_mod);

    paillier_ciphertext_t zeta;
    Polynomial<paillier_ciphertext_t> bzeta(dblocks, group_mod);

    gt_t sxi1_b, sxi2_b; gt_new(sxi1_b); gt_new(sxi2_b);
    gt_set_unity(sxi1_b); gt_set_unity(sxi2_b);

    bn_t tmpz; bn_new(tmpz);
    bn_t rpowsm; bn_new(rpowsm);
    bn_t rpows; bn_new(rpows);

    const int64_t step( std::ceil( (double)client.d/(double)nb_tasks) );
    const int64_t powsr_deg(std::max(server.W[0].degree(),step));
    const int64_t apows_deg(std::max(powsr_deg, std::max(server.bigblocks,dblocks)));

#ifdef DEBUG
    std::clog << "pows_r: d°" << powsr_deg << '|' << apows_deg << std::endl;
#endif

    {	// server computation : compute xi & zeta
        Chrono c_server; c_server.start();

            // server computation : first powers of eval. point for zeta
#ifdef VESPO_TIMINGS
        Chrono c_step; c_step.start();
#endif

        Polynomial<bn_t> pows_r(apows_deg,group_mod);
        geo_progression(pows_r, 0, r, powsr_deg, r, group_mod, time_r);

        bn_copy(rpowsm, pows_r[server.W[0].degree()]);	// last power of r;
        bn_mul(rpows, rpowsm, r);						// next power of r;
        bn_mod(rpows, rpows, group_mod);

        Polynomial<bn_t> revpow(step,group_mod);
        for(int64_t i=0; i<(step+1); ++i) {
            bn_copy(revpow[i],pows_r[step-i]);
        }

// server computation : compute xi
        if (client.d) {
            Polynomial<g1_t> Ti(client.d-1, group_mod);

            g1_horner_mxp_iter(Ti, server.S, r, 0, client.d, revpow, step,
                               group_mod, nb_tasks );

            pairing_double(sxi1_b, client.d, nb_tasks, group_mod, Ti, server.H1_b);
            pairing_double(sxi2_b, client.d, nb_tasks, group_mod, Ti, server.H2_b);

#ifdef VESPO_TIMINGS
            time_p = c_step.stop();
            std::clog << "  SERVER dotprod xi: " << time_p << " (" << (3*client.d) << " Horner/pairings operations)" << std::endl;
            c_step.start();
#endif
        }

            // server computation : compute zeta by blocks
#pragma omp parallel for shared(bzeta,server,pows_r)
        for(int64_t i=0; i<nbblocks; ++i)
                // server computation : zeta baby-steps
            paillier_hom_dp(bzeta[i], server.pub, server.W[i], pows_r, server.W[i].degree());

            // server computation : zeta giant-step
        geo_progression(pows_r, server.bigblocks, rpows, dblocks, rpowsm,
                        group_mod, time_r);
        paillier_hom_dp(zeta, server.pub, bzeta, pows_r, dblocks);

#ifdef VESPO_TIMINGS
        time_e = c_step.stop();
        std::clog << "  SERVER homo. zeta: " << time_e << " (" << (client.d+1) << " Paillier operations on ciphers)" << std::endl;
#endif
#ifdef VESPO_CHECKERS
        bn_mxp_dig(tmpz, r, server.W[0].degree()+1, group_mod);
        if(bn_cmp(tmpz, rpows) == RLC_EQ) {
            std::clog << "[BASE PAI] \033[1;32mOK\033[0m" << std::endl;
        } else {
            std::clog << "[BASE PAI] \033[1;31m****** FAIL ******\033[0m" << std::endl;
        }
#endif

        time_s = c_server.stop();
    } // End of server computations

#ifdef VESPO_CHECKERS
    check_ciph_horner(sxi1_b, client.d, server, server.H1_b, r);
    check_ciph_horner(sxi2_b, client.d, server, server.H2_b, r);
#endif


    gt_t verif1; gt_new(verif1);
    gt_t verif2; gt_new(verif2);

    {	// Client verification
        vector vc;
        Chrono c_client; c_client.start();

#ifdef VESPO_TIMINGS
        Chrono c_step; c_step.start();
#endif

        mat_geometric_sum(vc, r, client.msigma, client.d, client.vbeta, group_mod);

#ifdef VESPO_TIMINGS
        time_cg = c_step.stop();
        std::clog << "  CLIENT mat geosum: " << time_cg << std::endl;
        c_step.start();
#endif

            // Client: Deciphering Paillier to get z
        PAILLIER_DECIPHER(z,client.pub, client.prv, zeta);
        bn_mod(z.m,z.m,group_mod);

#ifdef VESPO_TIMINGS
        time_cp = c_step.stop();
        std::clog << "  CLIENT paillier D: " << time_cp << std::endl;
        c_step.start();
#endif

            // Client: dot-product verification
            // s-r
        bn_sub(tmpz, client.s, r);
        bn_mod(tmpz, tmpz, group_mod);

            // c <-- alpha D(zeta) + c
        vespo_addinmul(vc.first, client.valpha.first, z.m);
        bn_mod(vc.first, vc.first, group_mod);

            // c <-- alpha D(zeta) + c
        vespo_addinmul(vc.second, client.valpha.second, z.m);
        bn_mod(vc.second, vc.second, group_mod);

            // xi^{s-r} * e^{alpha D(zeta) + c}
            //       First possiblity: via exp_sim
//     gt_exp_sim(verif1, sxi1_b, tmpz, client.e_T, vc.first);
//     gt_exp_sim(verif2, sxi2_b, tmpz, client.e_T, vc.second);

            //       Second possiblity: hand made Shamir trick
//     gt_exp_shamir(verif1, sxi1_b, tmpz, client.e_T, vc.first);
//     gt_exp_shamir(verif2, sxi2_b, tmpz, client.e_T, vc.second);

            //       Third possiblity: no Shamir trick, but with gen, seems faster with BN254 Gt group ...
        gt_exp(sxi1_b, sxi1_b, tmpz);	// xi^{s-r}
        gt_exp(sxi2_b, sxi2_b, tmpz);	// xi^{s-r}
        gt_exp_gen(verif1, vc.first);				// e^{alpha D(zeta) + c}
        gt_exp_gen(verif2, vc.second);				// e^{alpha D(zeta) + c}
        gt_mul(verif1, verif1, sxi1_b);
        gt_mul(verif2, verif2, sxi2_b);

        time_c = c_client.stop();

#ifdef VESPO_TIMINGS
        time_ce = c_step.stop();
        std::clog << "  CLIENT simpow gt : " << time_ce << std::endl;
        time_c = time_cg+time_cp+time_ce; // Not count std::clog
#endif


    } // End of client computations

    const bool pass(
        ( gt_cmp(verif1, client.K1_bT) == RLC_EQ ) &&
        ( gt_cmp(verif2, client.K2_bT) == RLC_EQ )
        );

    if (! pass) {
        printf("\033[1;31mK1_bT :\033[0m "); gt_print(client.K1_bT); printf("\n");
        printf("\033[1;31mverif1:\033[0m "); gt_print(verif1); printf("\n");
        printf("\033[1;31mK2_bT :\033[0m "); gt_print(client.K2_bT); printf("\n");
        printf("\033[1;31mverif2:\033[0m "); gt_print(verif2); printf("\n");
    }

    bn_free(tmpz); bn_free(rpowsm); bn_free(rpows);
    gt_free(sxi1_b); gt_free(sxi2_b);
    gt_free(verif1); gt_free(verif2);
    bn_free(group_mod);
#ifdef VESPO_TIMINGS
    std::clog << "[Audit (Eval)] END" << std::endl;
#endif
    return pass;
}
