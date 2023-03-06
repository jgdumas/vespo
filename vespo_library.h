// ==========================================================================
// VESPo: Verified Evaluation of Secret Polynomials
// Reference: [ https://arxiv.org/abs/2110.02022
//              J-G. Dumas, A. Maignan, C. Pernet, D. S. Roche ]
// Authors: J-G Dumas
// Time-stamp: <06 Mar 23 18:30:21 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/****************************************************************
 * VESPO Library definitions
 ****************************************************************/

#ifndef _VESPO_LIBRARY_H_
#define _VESPO_LIBRARY_H_

#include <time.h>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <numeric>
#include <omp.h>

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "relic/relic.h"
#include "relic/relic_conf.h"
#ifdef __cplusplus
}
#endif

#ifdef DEBUG
# ifndef VESPO_CHECKERS
    // Automatic verification of subroutines
#  define VESPO_CHECKERS
# endif
# ifndef VESPO_SUB_TIMINGS
    // Detailed timings
#  define VESPO_SUB_TIMINGS
# endif
#endif

#ifdef VESPO_SUB_TIMINGS
# ifndef VESPO_TIMINGS
    // Overall timings
#  define VESPO_TIMINGS
# endif
#endif

#ifndef VESPO_RELIC_LIMIT_MAX_ALLOC
    // With ALLOC=AUTO RELIC can't allocate too much
#define VESPO_RELIC_LIMIT_MAX_ALLOC 8192
#endif

/****************************************************************
 * Median deviation and chronograph
 ****************************************************************/
#include <cassert>
#include <vector>
#include <algorithm>

#ifndef CLOCKTYPE
#  ifdef CLOCK_PROCESS_CPUTIME_ID
    /*  cpu time in the current process */
#    define CLOCKTYPE  CLOCK_PROCESS_CPUTIME_ID
    /* real time in the current process */
// #    define CLOCKTYPE  CLOCK_REALTIME
#  else
    /* this one should be appropriate to
               avoid errors on multiprocessors systems */
#    define CLOCKTYPE  CLOCK_MONOTONIC
#  endif
#endif

#define VESPO_NANO_FACTOR 1.0e9

struct Chrono {
    struct timespec begin_time,end_time;
    void start() { clock_gettime(CLOCKTYPE, &begin_time); }
    double stop() {
        clock_gettime(CLOCKTYPE, &end_time);
        double ttime(difftime(end_time.tv_sec, begin_time.tv_sec));
        return ttime += ((double) (end_time.tv_nsec - begin_time.tv_nsec) )/ VESPO_NANO_FACTOR;
    }
};

// Returns maximal deviation in percent (v is supposed sorted)
template<typename Vect> double mediandeviation(const Vect& v);


/****************************************************************
 * relic bn_t utilities
 ****************************************************************/

#include <iostream>
std::ostream& operator<<(std::ostream& out, const bn_t& numb);


// Init and set RELIC bignum
void vespo_init_set(bn_t& a, const bn_t b);
void vespo_init_set_ui(bn_t& a, const dig_t b);


// Set a to a + b * c.
void vespo_addinmul(bn_t& a, const bn_t b, const bn_t c);



//=====================================================================
// VeSPo data structures
//       Polynomial types
//=====================================================================
template<typename Element> void freeer(Element& e);

template<typename Element> void newer(Element& e);

template<typename Element> std::ostream& printer(std::ostream& out, const Element& e);


template<typename Element> class Polynomial;
template<typename Element> std::ostream& operator<<(std::ostream& out, const Polynomial<Element>& P);

template<typename Element> class Polynomial {
    Polynomial() = delete;
    Polynomial(const Polynomial<Element>&) = delete;
protected:
    const bn_t & mod;
    int64_t d;
    Element * data;
public:
    Polynomial(Polynomial<Element>&&) = default;
    Polynomial(const int64_t degree, const bn_t& modulus);
    ~Polynomial();

    const bn_t& modulus() const { return this->mod; }
    const int64_t& degree() const { return this->d; }

    const Element& operator[](size_t i) const ;
    Element& operator[](size_t i);

    friend std::ostream& operator<< <Element>(std::ostream& out, const Polynomial<Element>& P);
    void random();

    void eval(Element&, const Element&) const;
};


//=====================================================================
// VeSPo data structures
//       2 x 2 linear algebra
//=====================================================================

struct vector {
    bn_t first, second;
    void init();
    vector();
    ~vector();
    vector(const bn_t& f, const bn_t& s);
    vector(const vector& v) ;

    void randomize(const bn_t mod);
    void modin(const bn_t& r) ;
    void copy(const vector& v) ;
    void addin(const vector& v);
    bool areEqual(const vector& v) const ;
};

struct matrix {
    bn_t a, b, c, d;
    void init() ;
    ~matrix();
    matrix();
    matrix(const bn_t& Aa, const bn_t& Ab, const bn_t& Ac, const bn_t& Ad);
    matrix(const matrix& A);
    void randomize(const bn_t& mod);
    void mulin(const bn_t& r);
    void modin(const bn_t& r);
    void copy(const matrix& A);
};

//=====================================================================
// VeSPo: linear algebra routines
//=====================================================================

    // (t0 + t1 X)^2 mod p0 + p1 X + X^2
    // --> t0^2 + 2t0t1 X + t1^2 ( -p0 -p1 X)
    // --> (t0^2 -p0t1^2) + (2t0t1 -p1t1^2) X
void sqr_modcharp(bn_t& s0, bn_t& s1,
                  const bn_t& p0, const bn_t& p1,
                  const bn_t& t0, const bn_t& t1,
                  const bn_t& mod);


    // (s0 + s1 X) X mod p0 + p1 X + X^2
    // --> s0 X + s1 (-p0 -p1 X)
    // --> (-p0 s1) + (s0 - s1p1) X
void mulinx_modcharp(bn_t& s0, bn_t& s1,
                     const bn_t& p0, const bn_t& p1,
                     const bn_t& mod);

    // Computes s0+s1Z=Z^d modulo P(X)=p0+p1Z+Z^2
void TMMP(bn_t& s0, bn_t& s1,
          const bn_t& p0, const bn_t& p1,
          const int64_t d, const bn_t& mod);

    // Determinant of <<a|b>,<c|d>>
void moddet(bn_t& det,
            const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
            const bn_t& mod);

void moddet(bn_t& det, const matrix& m, const bn_t& mod) {
    return moddet(det,m.a,m.b,m.c,m.d,mod);
}


    // Characteristic polynomial of <<a|b>,<c|d>>
    // (ad-bc) - (a+d) X + X^2
void modcharpoly(bn_t& p0, bn_t& p1,
                 const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
                 const bn_t& mod);

    // <u0,u1> <-- <<a|b>,<c|d>> . <v0,v1>
void matvectmod(bn_t& u0, bn_t& u1,
                const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
                const bn_t& v0, const bn_t& v1,
                const bn_t& mod);

void matvectmod(vector& u, const matrix& m, const vector& v, const bn_t& mod) {
    matvectmod(u.first, u.second, m.a, m.b, m.c, m.d, v.first, v.second, mod);
}


    // Projected matrix geometric sum
    // u <-- sum_{i=0}^k M^i . v
void PMGS(vector& u, const matrix& M,
          const uint64_t k,
          const vector& v, const bn_t& mod);

    // Projected matrix geometric sum
    // c <-- sum_{i=0}^k (rA)^i . b
void mat_geometric_sum(vector& c,
                       const bn_t& r, const matrix& A, const uint64_t k,
                       const vector& b,
                       const bn_t& mod);


/****************************************************************
 * Paillier utilities
 ****************************************************************/
typedef struct {
    shpe_t rlc;
    bn_t nsq;
} paillier_pubkey_t;

#define paillier_prvkey_t shpe_t

#define paillier_freeprvkey(a) { shpe_free(a); }
#define paillier_freepubkey(a) { shpe_free((a).rlc); bn_free((a).nsq); }
#define paillier_sizepubkey(a) ((bn_size_bin((a).rlc->crt->n)+bn_size_bin((a).nsq))<<3)
#define paillier_sizeprvkey(a) ((bn_size_bin((a)->crt->p)+bn_size_bin((a)->crt->q)+bn_size_bin((a)->crt->dp)+bn_size_bin((a)->crt->dq)+bn_size_bin((a)->crt->qi))<<3)

struct paillier_plaintext_t {
    bn_t m;
    paillier_plaintext_t() { bn_null(m); bn_new(m); }
    paillier_plaintext_t(bn_t r) : paillier_plaintext_t() { bn_copy(m,r); }
    ~paillier_plaintext_t() { bn_free(m); }
};

struct paillier_ciphertext_t {
    bn_t c;
    paillier_ciphertext_t() { bn_null(c); bn_new(c);}
    paillier_ciphertext_t(bn_t r) : paillier_ciphertext_t() { bn_copy(c,r); }
    ~paillier_ciphertext_t() { bn_free(c); }
};

	// a*b mod kpub.nsq
int paillier_mul(paillier_ciphertext_t& res,
                 const paillier_ciphertext_t& a,
                 const paillier_ciphertext_t& b,
                 const paillier_pubkey_t& kpub);

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

void random_precompute(Polynomial<bn_t>& reusedrands,
                       const paillier_pubkey_t& pub,
                       const paillier_prvkey_t& prv,
                       const int64_t psize);


#define PAILLIER_DECIPHER(res,pub,prv,cip) cp_shpe_dec((res).m,(cip).c,prv)

#ifndef VESPO_NOTSECURE
#  define PAILLIER_PRECOMPU(pub,prv,psize)
#  define PAILLIER_ENCIPHER(res,mes,pub,prv,ii) cp_shpe_enc_prv((res).c,mes,prv)
#else
#  define PAILLIER_PRECOMPU(pub,prv,psize)									\
	 const size_t paillierands(std::min((int64_t)VESPO_NOTSECURE,psize));	\
	 Polynomial<bn_t> reusedrands(paillierands-1,pub.nsq);					\
	 random_precompute(reusedrands,pub,prv,psize);

// Compute directly (1+ mes * N) mod N^2
// Instead of (1+N)^m mod N^2  (would be: bn_mxp(res.c, g, mes, pub->nsq);)
// Therefore almost no exponentiations at setup
#  define PAILLIER_ENCIPHER(res,mes,pub,prv,ii)				\
     bn_mul(res.c, prv->crt->n, mes);						\
     bn_mod(res.c, res.c, pub.nsq);							\
     bn_mul(res.c, res.c, prv->b);							\
     bn_add_dig(res.c, res.c, 1);							\
     bn_mod(res.c, res.c, pub.nsq);							\
	 bn_mul(res.c, res.c, reusedrands[ii % VESPO_NOTSECURE]);\
	 bn_mod(res.c, res.c, pub.nsq);

#endif
//=====================================================================

//=====================================================================
// VeSPo data structures
//       Client/Server
//=====================================================================

struct client_t {
    g1_t g1;
    g2_t g2;
    gt_t e_T;

    int64_t d;
    paillier_pubkey_t pub;
    paillier_prvkey_t prv;
    bn_t s;
    gt_t K1_bT, K2_bT;
    vector valpha, vbeta;
    matrix msigma;

    client_t(const int64_t degree);

    ~client_t();

    uint64_t private_size(uint64_t elements_size) const ;
    uint64_t privpub_size(uint64_t elements_size) const ;

    void keygen(uint64_t pailliersize);
};

struct server_t {
    paillier_pubkey_t pub;
    std::vector<Polynomial<paillier_ciphertext_t>> W;
    Polynomial<g2_t> H1_b, H2_b;
    Polynomial<g1_t> S;
    int64_t bigblocks;
    server_t(const int64_t degree, const int64_t blocks, const bn_t& modulus);
    ~server_t() {}
};


//=====================================================================
// VeSPo: checkers
//=====================================================================

#ifdef VESPO_CHECKERS
void check_dp(g1_t& dp, const Polynomial<g1_t>& gS, const Polynomial<bn_t>& P);

bool check_orders(const client_t& client, const bn_t& pairing_r);

bool check_pubkey(const client_t& client, const server_t& server,
                  const Polynomial<bn_t>& P_b,
                  const bn_t& eval, const g2_t& K_b);

bool check_g2_omxp(const g2_t& T, const Polynomial<g2_t>& S,
                   const bn_t& r, const int64_t deg, const bn_t& mod) ;

bool check_g1_hmxp(const Polynomial<g1_t>& T, const Polynomial<g1_t>& S,
                   const bn_t& r, const int64_t from, const int64_t length,
                   const bn_t& mod);

void check_g1_msl(g1_t& r, const g1_t* P, const bn_t* K, int N);

void check_g2_mul_sim_lot(g2_t& RES, const g2_t* P, const bn_t* K, int N);

bool check_ciph_horner(gt_t sxi_b, const int64_t deg, const server_t& server,
                       const Polynomial<g1_t>& H_b, const bn_t& rmod);
#endif

//=====================================================================
// VeSPo: masking
//=====================================================================

    // <P1,P2> <-- sum X^i ( p_i valpha + msigma^i vbeta)
void vhide_poly(Polynomial<bn_t> & P1, Polynomial<bn_t>& P2,
                const vector & valpha, const vector& vbeta,
                const matrix& msigma, const Polynomial<bn_t> & P);

//=====================================================================
// VeSPo: Paillier
//=====================================================================

    // Computes 1, r1, ..., r1^bigsdegree,
    // r1^bigsdegree r2, ..., r1^bigsdegree r2^(totaldegree-bigsdegree)
void geo_progression(Polynomial<bn_t>& pows_r,
                     const int64_t bigsdegree, const bn_t& r1,
                     const int64_t totaldegree, const bn_t& r2,
                     const bn_t& mod, double& time_r);

	// Paillier homomorphic dot-product
void paillier_hom_dp(paillier_ciphertext_t& eval,
                     const paillier_pubkey_t& kpub,
                     const Polynomial<paillier_ciphertext_t>& Wi,
                     const Polynomial<bn_t>& pows_r,
                     const int64_t degree);

    // W <-- E(P), and W is cut into chunks
void encrypt_poly(std::vector<Polynomial<paillier_ciphertext_t>>& W,
                  const Polynomial<bn_t>& P,
                  const paillier_pubkey_t& pub, const paillier_prvkey_t& prv);


//=====================================================================
// VeSPo: in the exponents
//=====================================================================

    // hi[i] <-- g2^{p_(i+1)}
Polynomial<g2_t>& powered_poly_gen2(Polynomial<g2_t>& H,
                                    const Polynomial<bn_t>& P);

    // S[i] <-- g1^{s^i}
Polynomial<g1_t>& power_progression_gen(Polynomial<g1_t>& S, const bn_t& s,
                                        const int64_t d, const bn_t& mod);


//========================================================
// VeSPo: homomorphic operations
//========================================================

	// Parallel simultaneous cipher/clear multiplication
void g1_mul_sim_lot_par(g1_t& RES, const g1_t* P, const bn_t* K, int N,
                        const bn_t& mod, const int nbtasks);

    // Sequential Server xi computation
    //   first step : Horner-like prefix computations of S_j^{x_k}
void g1_horner_mxp_seq(Polynomial<g1_t>& U, const Polynomial<g1_t>& S,
                       const bn_t& r, const int64_t from, const int64_t length);

    // Parallel Server xi computation
    //   first step : Horner-like prefix computations of S_j^{x_k}
void g1_horner_mxp_iter(Polynomial<g1_t>& U, const Polynomial<g1_t>& S,
                        const bn_t& r, const int64_t from, const int64_t length,
                        const Polynomial<bn_t>& revpow, const int64_t step,
                        const bn_t& mod, const int64_t nbtasks);

	// Simultaneous pairing computation (and final accumulation)
void pairing_sim(gt_t& xi_b, const g1_t* H, const g2_t* T, const int64_t deg);

    // Server xi computation
    //   second step: applying the pairing maps as a dotproduct
void pairing_double(gt_t& xi_b, const int64_t length,
                    const int64_t nbblocks, const bn_t& mod,
                    const Polynomial<g1_t>& H, const Polynomial<g2_t>& Ti);


//========================================================
// VeSPo: Init protocol
//========================================================
void setup(client_t& client, server_t& server,
           uint64_t pailliersize,
           const Polynomial<bn_t>& P, double& t_horner);

//========================================================
// VeSPo: Update protocol
//========================================================
bool update(client_t& client, server_t& server,
            Polynomial<bn_t>& P, const bn_t& delta, const size_t index,
            double& time_u);

//========================================================
// VeSPo: Audit protocol
//========================================================
bool eval(paillier_plaintext_t& z, const client_t& client,
          const server_t& server,
          const bn_t& r, const int64_t nb_tasks,
          double& time_c, double& time_s
#ifdef VESPO_TIMINGS
          , double& time_e, double& time_p,
          double& time_ce, double& time_cg, double& time_cp, double &time_cc
#endif
          );

#endif
