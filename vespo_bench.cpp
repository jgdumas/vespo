/****************************************************************
 * Paramétrés RELIC: ccmake relic-target
 * ARITH                            x64-asm-4l
 * BN_PRECI                         4096
 * FP_PRIME                         254
 * see relic-git-src/cmake/bn.cmake: set(BN_PRECI 4096)
 * see relic-git-src/cmake/fp.cmake: set(FP_PRIME 254)
 ****************************************************************/

/****************************************************************
 * Parameters RELIC: ccmake relic-target
 * ARITH                            gmp
 * BN_PRECI                         4096
 * FP_PRIME                         256
 * FP_QNRES                         OFF
 * see relic-git-src/cmake/bn.cmake: set(BN_PRECI 4096)
 * see relic-git-src/cmake/fp.cmake: set(FP_PRIME 256)
 ****************************************************************/

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
#ifdef __cplusplus
}
#endif

// If using a version of relic without simultaneous bn_mxp
#ifndef BN_XPWDT
			// Size of Generalized Shamir trick
# define BN_XPWDT 7
# include "relic_bn_mxp_sim.c"
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
			// With ALLOC=AUTO can't allocate too much
#define VESPO_RELIC_LIMIT_MAX_ALLOC 8192
#endif

/****************************************************************
 * Median deviation and chronograph
 ****************************************************************/
#include <cassert>
#include <vector>
#include <algorithm>
template<typename Vect> double mediandeviation(const Vect& v) {
    assert(v.size()>0);
    typename Vect::value_type median(v[v.size()/2]);
    double t1( median-v.front() );
    double t2( v.back()-median );
    return 100.*std::max(t1,t2)/median;
}

#ifndef CLOCKTYPE
#  ifdef CLOCK_PROCESS_CPUTIME_ID
			/* cpu time in the current process */
#    define CLOCKTYPE  CLOCK_PROCESS_CPUTIME_ID
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

/****************************************************************
 * relic utilities
 ****************************************************************/

#include <iostream>
std::ostream& operator<<(std::ostream& out, const bn_t& numb) {
    size_t bits = bn_size_str(numb, 10);
    char str[RLC_BN_BITS + 2];
    bn_write_str(str, bits, numb, 10);
    return out << str;
}

void vespo_init_set(bn_t& a, const bn_t b) {
    bn_null(a); bn_new(a); bn_copy(a,b);
}

void vespo_init_set_ui(bn_t& a, const dig_t b) {
    bn_null(a); bn_new(a); bn_set_dig(a,b);
}

// Set a to a + b * c.
void vespo_addinmul(bn_t& a, const bn_t b, const bn_t c) {
    bn_t tmp; bn_new(tmp);
    bn_mul(tmp,b,c);
    bn_add(a,a,tmp);
    bn_free(tmp);
}


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
#define paillier_sizeprvkey(a) ((bn_size_bin((a)->crt->n)+bn_size_bin((a)->crt->p)+bn_size_bin((a)->crt->q)+bn_size_bin((a)->crt->dp)+bn_size_bin((a)->crt->dq)+bn_size_bin((a)->crt->qi))<<3)


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

int paillier_mul(paillier_ciphertext_t& res,
                 const paillier_ciphertext_t& a,
                 const paillier_ciphertext_t& b,
                 const paillier_pubkey_t& kpub) {
	int result = RLC_OK;

    bn_mul(res.c, a.c, b.c);
    bn_mod(res.c, res.c, kpub.nsq);

	return result;
}

void paillier_exp(paillier_ciphertext_t& res,
                  const bn_t w, const bn_t r,
                  const paillier_pubkey_t& kpub) {
    bn_mxp(res.c, w, r, kpub.nsq);
}

void paillier_exp_sim(paillier_ciphertext_t& res,
                      const bn_t w[BN_XPWDT], const bn_t r[BN_XPWDT],
                      const paillier_pubkey_t& kpub) {

    bn_mxp_sim(res.c, w, r, kpub.nsq);
}

//=====================================================================
// VeSPo data structures
//       Polynomial types
//=====================================================================
template<typename Element> void freeer(Element& e);
template<> void freeer<bn_t>(bn_t& e) { bn_free(e); }
template<> void freeer<g1_t>(g1_t& e) { g1_free(e); }
template<> void freeer<g2_t>(g2_t& e) { g2_free(e); }
template<> void freeer<gt_t>(gt_t& e) { gt_free(e); }
template<> void freeer<paillier_ciphertext_t>(paillier_ciphertext_t& e) { }

template<typename Element> void newer(Element& e);
template<> void newer<bn_t>(bn_t& e) { bn_null(e); bn_new(e); }
template<> void newer<paillier_ciphertext_t>(paillier_ciphertext_t& e) { }
template<> void newer<g1_t>(g1_t& e) { g1_null(e); g1_new(e); }
template<> void newer<g2_t>(g2_t& e) { g2_null(e); g2_new(e); }
template<> void newer<gt_t>(gt_t& e) { gt_null(e); gt_new(e); }

template<typename Element> std::ostream& printer(std::ostream& out, const Element& e);
template<> std::ostream& printer<bn_t>(std::ostream& o, const bn_t& e) { return o<<e;}
template<> std::ostream& printer<g1_t>(std::ostream& o, const g1_t& e) { g1_print(e); return o;}
template<> std::ostream& printer<g2_t>(std::ostream& o, const g2_t& e) { g2_print(e); return o;}
template<> std::ostream& printer<gt_t>(std::ostream& o, const gt_t& e) { gt_print(e); return o;}

template<typename Element> class Polynomial {
    Polynomial() = delete;
    Polynomial(const Polynomial<Element>&) = delete;
protected:
    const bn_t & mod;
    int64_t d;
    Element * data;
public:
    Polynomial(Polynomial<Element>&&) = default;
    Polynomial(const int64_t degree, const bn_t& modulus) : mod(modulus), d(degree), data(new Element[degree+1]) {
        for (int64_t i = 0; i <= degree; ++i) {
            newer(this->data[i]);
        }
    }
    ~Polynomial() {
        for (int64_t i = 0; i <= this->d; ++i) {
            freeer(this->data[i]);
        }
        delete [] data;
    }

    const bn_t& modulus() const { return this->mod; }
    const int64_t& degree() const { return this->d; }

    virtual const Element& operator[](size_t i) const {
#ifdef DEBUG
        if ( (int64_t)i > this->degree()) {
            std::cerr << "****************** ERROR ******************* Polynomial access at " << i << " out of " << this->degree() << std::endl;
            throw(1);
        }
#endif
        return this->data[i];
    }
    virtual Element& operator[](size_t i) {
#ifdef DEBUG
        if ( (int64_t)i > this->degree()) {
            std::cerr << "****************** ERROR ******************* Polynomial access at " << i << " out of " << this->degree() << std::endl;
            throw(1);
        }
#endif
        return this->data[i];
    }

    friend std::ostream& operator<<(std::ostream& out, const Polynomial<Element>& P) {
        printer(out,P.data[0]) << " + ";
        for (int64_t i = 1; i < P.d; ++i) {
            printer(out,P.data[i]) << "X^" << i << " + ";
        }
        return printer(out,P.data[P.d]) << "X^" << P.d;
    }

    void random();

    void eval(Element&, const Element&) const;
};

template<typename Element> class Reverse : public Polynomial<Element> {
public:
    using  Polynomial<Element>::Polynomial;

    virtual const Element& operator[](size_t i) const { return this->data[this->d-i]; }
    virtual Element& operator[](size_t i) { return this->data[this->d-i]; }
};




template<> void Polynomial<bn_t>::random() {
    std::clog << "[Rand P] BEG: " << this->d << ", mod " << this->modulus() << std::endl;
    for(int64_t i = 0; i <= this->d; ++i) {
        bn_rand_mod(this->data[i], this->modulus());
    }
    std::clog << "[Rand P] END" << std::endl;
}

template<> void Polynomial<bn_t>::eval(bn_t& eval, const bn_t& s) const {
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
                       const int64_t psize) {
    if (reusedrands.degree()<(psize-1)) {
        std::clog << "*WARNING WARNING* Paillier randomness reused\n"
                  << "*WARNING WARNING* " << (reusedrands.degree()+1)
                  << " instead of " << psize
                  << "\n*WARNING WARNING* insecure, only benchmarking"
                  << std::endl;
    }
#pragma omp parallel for
    for (int64_t i = 0; i <= reusedrands.degree(); ++i) {
		bn_rand_mod(reusedrands[i], prv->a);
		bn_mxp(reusedrands[i], prv->gn, reusedrands[i], pub.nsq);
    }
}




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
//       2 x 2 linear algebra
//=====================================================================

struct vector {
    bn_t first, second;
    void init() {
        bn_null(this->first);
        bn_null(this->second);
        bn_new(this->first);
        bn_new(this->second);
    }
    vector() { init(); }
    ~vector() {
        bn_free(this->first);
        bn_free(this->second);
    }
    vector(const bn_t& f, const bn_t& s) : vector() {
        bn_copy(this->first, f);
        bn_copy(this->second,s);
    }
    vector(const vector& v) : vector() {
        this->copy(v);
    }
    void randomize(const bn_t mod) {
        bn_rand_mod(this->first, mod);
        bn_rand_mod(this->second, mod);
    }
    void modin(const bn_t& r) {
        bn_mod(this->first, this->first, r);
        bn_mod(this->second, this->second, r);
    }
    void copy(const vector& v) {
        bn_copy(this->first, v.first);
        bn_copy(this->second,v.second);
    }
    void addin(const vector& v) {
        bn_add(this->first, this->first, v.first);
        bn_add(this->second, this->second, v.second);
    }
    bool areEqual(const vector& v) const {
        return ( (bn_cmp(this->first,v.first)== RLC_EQ) &&
                 (bn_cmp(this->second,v.second)== RLC_EQ) );
    }
};

struct matrix {
    bn_t a, b, c, d;
    void init() {
        bn_null(a);bn_null(b);bn_null(c);bn_null(d);
        bn_new(a); bn_new(b); bn_new(c); bn_new(d);
    }
    ~matrix() { bn_free(a); bn_free(b); bn_free(c); bn_free(d); }

    matrix() { init(); }

    matrix(const bn_t& Aa, const bn_t& Ab, const bn_t& Ac, const bn_t& Ad) : matrix() {
        bn_copy(a, Aa);
        bn_copy(b, Ab);
        bn_copy(c, Ac);
        bn_copy(d, Ad);
    }
    matrix(const matrix& A) : matrix(A.a,A.b,A.c,A.d) {}

    void randomize(const bn_t& mod) {
        bn_rand_mod(a, mod);
        bn_rand_mod(b, mod);
        bn_rand_mod(c, mod);
        bn_rand_mod(d, mod);
    }
    void mulin(const bn_t& r) {
        bn_mul(a, a, r);
        bn_mul(b, b, r);
        bn_mul(c, c, r);
        bn_mul(d, d, r);
    }
    void modin(const bn_t& r) {
        bn_mod(a, a, r);
        bn_mod(b, b, r);
        bn_mod(c, c, r);
        bn_mod(d, d, r);
    }
    void copy(const matrix& A) {
        bn_copy(a, A.a);
        bn_copy(b, A.b);
        bn_copy(c, A.c);
        bn_copy(d, A.d);
    }
};

std::ostream& operator<<(std::ostream& out, const matrix& m) {
    return out << "<<" << m.a << '|' << m.b << ">,<" << m.c << '|' << m.d << ">>";
}

std::ostream& operator<<(std::ostream& out, const vector& v) {
    return out << '<' << v.first << ',' << v.second << '>';
}


//=====================================================================
// VeSPo: linear algebra
//=====================================================================

void sqr_modcharp(bn_t& s0, bn_t& s1,
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


void mulinx_modcharp(bn_t& s0, bn_t& s1,
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

void TMMP(bn_t& s0, bn_t& s1,
          const bn_t& p0, const bn_t& p1,
          const int64_t d, const bn_t& mod) {
        // Computes s0+s1Z=Z^d modulo P(X)=p0+p1Z+Z^2
#ifdef DEBUG
    std::clog << "[2-MMP] " << d << " BEG" << std::endl;
#endif
    if (d<=1) {
        if (d==1) {
            vespo_init_set_ui(s0,0);
            vespo_init_set_ui(s1,1);
            return;
        }
        if (d==0) {
            vespo_init_set_ui(s0,1);
            vespo_init_set_ui(s1,0);
            return;
        }
        if (d<0) {
            std::cerr << "**** ERROR **** negative exponent not implemented\n";
            return;
        }
    }

    bn_t t0, t1; bn_new(t0); bn_new(t1);
    TMMP(t0, t1, p0, p1, d>>1, mod);			// Recursive call

    sqr_modcharp(s0, s1, p0, p1, t0, t1, mod);	// T^2		mod P
    if (d & 0x1) {
        mulinx_modcharp(s0, s1, p0, p1, mod);	// X T^2	mod P
    }

    bn_free(t0);bn_free(t1);
#ifdef DEBUG
    std::clog << "[2MMP] " << d << " END" << std::endl;
#endif
}


void moddet(bn_t& det,
            const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
            const bn_t& mod) {
        // Determinant of <<a|b>,<c|d>>
    bn_t tmp; bn_new(tmp);

    bn_mul(det, a, d);		// ad
    bn_mul(tmp, b, c);		// bc
    bn_sub(det, det, tmp);	// ad-bc
    bn_mod(det, det, mod);

    bn_free(tmp);
}

void moddet(bn_t& det, const matrix& m, const bn_t& mod) {
    return moddet(det,m.a,m.b,m.c,m.d,mod);
}


void modcharpoly(bn_t& p0, bn_t& p1,
                 const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
                 const bn_t& mod) {
        // Characteristic polynomial of <<a|b>,<c|d>>
        // (ad-bc) - (a+d) X + X^2

    bn_t tmp; bn_new(tmp);

    bn_add(p1, a, d);	// a+d
    bn_neg(p1, p1);
    bn_mod(p1, p1, mod);

    bn_free(tmp);

    moddet(p0, a, b, c, d, mod);
}

void matvectmod(bn_t& u0, bn_t& u1,
                const bn_t& a, const bn_t& b, const bn_t& c, const bn_t& d,
                const bn_t& v0, const bn_t& v1,
                const bn_t& mod) {
        // <u0,u1> <-- <<a|b>,<c|d>> . <v0,v1>
#ifdef DEBUG
    std::clog << matrix(a,b,c,d) << '.' << vector(v0,v1) << " - ";
#endif
	// < a v0 + b v1 , c v0 + d v1 >
    bn_t tmp; bn_new(tmp);

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

void matvectmod(vector& u, const matrix& m, const vector& v, const bn_t& mod) {
    matvectmod(u.first, u.second, m.a, m.b, m.c, m.d, v.first, v.second, mod);
}


void PMGS(vector& u, const matrix& M,
          const uint64_t k,
          const vector& v, const bn_t& mod) {
        // Projected matrix geometric sum
        // u <-- sum_{i=0}^k M^i . v
#ifdef DEBUG
    std::clog << "[PMGS] " << k << " BEG" << std::endl;
#endif
	// Characteristic polynomial
    bn_t p0, p1; bn_new(p0); bn_new(p1);
    modcharpoly(p0, p1, M.a, M.b, M.c, M.d, mod);

	// (A-I2)^{-1} --> 1/det < < d-1 | -b >, < -c | a-1 > >
    bn_t dm, nb, nc, am; bn_new(dm); bn_new(nb); bn_new(nc); bn_new(am);
    bn_sub_dig(dm, M.d, 1);
    bn_neg(nb, M.b);
    bn_neg(nc, M.c);
    bn_sub_dig(am, M.a, 1);

	// Determinant( A-I2 ) --> am dm - b c
    bn_t idet; bn_new(idet);
    moddet(idet, am, M.b, M.c, dm, mod);
    bn_mod_inv(idet, idet, mod);

#ifdef DEBUG
    std::clog << "(" << M << "-IdentityMatrix(2))." << matrix(dm, nb, nc, am) << " mod " << mod << ';' << std::endl;
#endif

	// (A-I2)^{-1} v
    bn_t w0, w1; bn_new(w0); bn_new(w1);
    matvectmod(w0, w1, dm, nb, nc, am, v.first, v.second, mod);
    bn_mul(w0, w0, idet);
    bn_mod(w0, w0, mod);
    bn_mul(w1, w1, idet);
    bn_mod(w1, w1, mod);

#ifdef DEBUG
    std::clog << "1/(" << M << "-IdentityMatrix(2))." << v << " - " << vector(w0,w1) << " mod " << mod << ';' << std::endl;
#endif


	// Polynomial Geometric sum
    bn_t f0, f1; bn_new(f0); bn_new(f1);
    TMMP(f0, f1, p0, p1, k+1, mod);

	// F(A) - I2 --> << a*f1+f0-1 | f1*b >, < f1*c | d*f1+f0-1 >>
    bn_sub_dig(dm, f0, 1);	// f0-1
    bn_mul(am, M.a, f1);	// a*f1
    bn_mod(am, am, mod);
    bn_add(am, am, dm);		// a*f1+f0-1

    bn_mul(nb, M.d, f1);	// d*f1
	bn_mod(nb, nb, mod);
    bn_add(dm, nb, dm);		// d*f1+f0-1

    bn_mul(nb, M.b, f1);	// f1*b
    bn_mod(nb, nb, mod);

    bn_mul(nc, M.c, f1);	// f1*c
    bn_mod(nc, nc, mod);

    matvectmod(u.first, u.second, am, nb, nc, dm, w0, w1, mod);


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
            std::clog << "[ PGMS ] \033[1;32mOK\033[0m" << std::endl;
        } else {
            std::cerr << "[ PGMS ] \033[1;31m****** FAIL ******\033[0m" << std::endl;
        }
    }
#endif

    bn_free(p0);bn_free(p1);
    bn_free(dm);bn_free(nb);bn_free(nc);bn_free(am);
    bn_free(idet);
    bn_free(w0);bn_free(w1);
    bn_free(f0);bn_free(f1);
#ifdef DEBUG
    std::clog << "[PMGS] " << k << " END" << std::endl;
#endif
}


void mat_geometric_sum(vector& c,
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

    client_t(const int64_t degree) : d(degree) {
        g1_null(this->g1);g1_new(this->g1);
        g2_null(this->g2);g2_new(this->g2);
        gt_null(this->e_T);gt_new(this->e_T);
        bn_null(this->pub.nsq); bn_new(this->pub.nsq);
        shpe_null( this->pub.rlc ); shpe_new( this->pub.rlc );
        shpe_null( this->prv ); shpe_new( this->prv );
        bn_null(this->s);bn_new(this->s);
        gt_null(this->K1_bT);gt_new(this->K1_bT);
        gt_null(this->K2_bT);gt_new(this->K2_bT);
    }
    ~client_t() {
        g1_free(this->g1);
        g2_free(this->g2);
        gt_free(this->e_T);
        paillier_freeprvkey(this->prv);
        paillier_freepubkey(this->pub);
        bn_free(this->s);
        gt_free(this->K1_bT);
        gt_free(this->K2_bT);
    }

    uint64_t private_size(uint64_t elements_size) const {
            // Groups, prime order and pubkey are supposed public data
            // The degree of the polynomial must be stored
        uint64_t size=0;
        size += sizeof(this->d)*8;	// degree
        size += paillier_sizepubkey( this->pub ); // pubkey
        size += paillier_sizeprvkey( this->prv ); // prvkey
            // s,K_bT,valpha,vbeta,msigma + MerkleTreeRoot
        size += elements_size*(1+2+2+2+4 + 1);
        return size;
    }
	//=====================================================================
	// VeSPo: key generation
	//=====================================================================
    void keygen(uint64_t pailliersize) {
        pc_get_ord( prv->a );

        cp_shpe_gen( pub.rlc, prv, bn_bits( prv->a ), pailliersize);
        bn_sqr(pub.nsq, pub.rlc->crt->n); // store n^2 too
    }
};

struct server_t {
    paillier_pubkey_t pub;
    std::vector<Polynomial<paillier_ciphertext_t>> W;
    Polynomial<g2_t> H1_b, H2_b;
    Polynomial<g1_t> S;
    int64_t bigblocks;
    server_t(const int64_t degree, const int64_t blocks, const bn_t& modulus) :
            H1_b(degree-1,modulus), H2_b(degree-1,modulus),
            S(degree, modulus) {
            // Cut W into 'blocks' chunks
            // last ones of size [ (degree+1)/Blocks ]
            // first ones of that size + 1
            // so that sum of sizes is 'degree+1'
        int64_t sW(degree+1);
        W.reserve(blocks);
        int64_t sdeg(sW/blocks);
        bigblocks = (sW-sdeg*blocks);
        const int64_t ldeg(sdeg); --sdeg;
        int64_t i=0; for(; i<bigblocks;++i) {
            W.emplace_back(ldeg, modulus);
        }
        for( ; i<blocks; ++i){
            W.emplace_back(sdeg, modulus);
        }
            // If no bigblocks, actually all are of same (big) size
        if (bigblocks == 0) bigblocks=blocks;
    }
    ~server_t() {}
};


//=====================================================================
// VeSPo: checkers
//=====================================================================

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



//=====================================================================
// VeSPo: checkers implementations
//=====================================================================

void check_dp(g1_t& dp, const Polynomial<g1_t>& gS, const Polynomial<bn_t>& P) {
        // dp <-- \oplus [ p_i ] gS_i
    g1_t tmp; g1_new(tmp);
    g1_mul(dp, gS[0], P[0]);
    for(int64_t i=1; i <= P.degree(); ++i) {
        g1_mul(tmp, gS[i], P[i]);
        g1_add(dp, dp, tmp);
    }
    g1_free(tmp);
}


bool check_orders(const client_t& client, const bn_t& pairing_r) {
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


bool check_pubkey(const client_t& client, const server_t& server,
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

bool check_g2_omxp(const g2_t& T, const Polynomial<g2_t>& S,
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


bool check_g1_hmxp(const Polynomial<g1_t>& T, const Polynomial<g1_t>& S,
                   const bn_t& r, const int64_t from, const int64_t length,
                   const bn_t& mod) {

    const int64_t deg(from+length-1);
    bool pass(true);
    Polynomial<g1_t> U(deg,mod);
    g1_copy(U[0], S[0]);
    if (g1_cmp(U[0], T[0]) != RLC_EQ) {
        printer(printer(std::cerr << "check_g1_hmxp: U[0] != T[0]: ", U[0]) << ' ', T[0]) << std::endl;
        pass = false;
    }
    for (int64_t i = 1; i <= deg; ++i) {
        g1_mul(U[i], U[i-1], r);		// t^r
        g1_add(U[i], U[i], S[i]);		// S_{i-1} t^r
        if ( (i>= from) && (g1_cmp(U[i], T[i]) != RLC_EQ)) {
            printer(printer(std::cerr << "check_g1_hmxp: U[" << i << "] != T[" << i << "]: ", U[i]) << ' ', T[i]) << std::endl;
            pass = false;
        }
    }
    return pass;
}


void check_g1_msl(g1_t& r, const g1_t* P, const bn_t* K, int N) {
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

void check_g1_mul_sim_lot(g1_t& RES, const g1_t* P, const bn_t* K, int N) {
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

bool check_ciph_horner(gt_t sxi_b, const int64_t deg, const server_t& server,
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

void vhide_poly(Polynomial<bn_t> & P1, Polynomial<bn_t>& P2,
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

void geo_progression(Polynomial<bn_t>& pows_r,
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


void paillier_hom_dp(paillier_ciphertext_t& eval,
                     const paillier_pubkey_t& kpub,
                     const Polynomial<paillier_ciphertext_t>& Wi,
                     const Polynomial<bn_t>& pows_r,
                     const int64_t degree) {
        // eval <-- prod Wi[i]^pows_r[i]
#ifdef VESPO_CHECKERS
    std::clog << "[PailCdot] BEG: d°" << degree << std::endl;
#endif

    bn_mxp_sim_lot(eval.c, reinterpret_cast<const bn_t*>( &(Wi[0])),
                   reinterpret_cast<const bn_t*>(&(pows_r[0])),
                   kpub.nsq, degree+1);

#ifdef VESPO_CHECKERS
    std::clog << "[PailCdot] END" << std::endl;
#endif
}

void encrypt_poly(std::vector<Polynomial<paillier_ciphertext_t>>& W,
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

Polynomial<g2_t>& powered_poly_gen2(Polynomial<g2_t>& H,
                                    const Polynomial<bn_t>& P) {
        // hi[i] <-- g2^{p_(i+1)}
    std::clog << "[Powered] BEG: d°" << P.degree() << std::endl;
#ifdef VESPO_SUB_TIMINGS
    Chrono c_powpol; c_powpol.start();
#endif

#pragma omp parallel for
    for (int64_t i = 0; i < P.degree(); ++i) {
        g2_mul_gen(H[i], P[i+1]);	// g2^{p_(i+1)}
    }

#ifdef VESPO_SUB_TIMINGS
    std::clog << "  SETUP Enpowered poly.: " << c_powpol.stop() << " (" << (P.degree()+1) << " operations)" << std::endl;
#endif
    std::clog << "[Powered] END" << std::endl;
    return H;
}

Polynomial<g1_t>& power_progression_gen(Polynomial<g1_t>& S, const bn_t& s,
                                        const int64_t d, const bn_t& mod) {
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
#pragma omp parallel for
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
void setup(client_t& client, server_t& server,
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
bool update(client_t& client, server_t& server,
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
    bn_t da1, da2; bn_new(da1); bn_null(da1); bn_new(da2); bn_null(da2);
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

    g2_t Deltaj; g2_new(Deltaj); g2_new(Deltaj);
    gt_t DeltagT; gt_new(DeltagT);

    g2_mul(Deltaj, Hda1, si);
    pc_map(DeltagT, client.g1, Deltaj);
    gt_mul(client.K1_bT, client.K1_bT, DeltagT);

    g2_mul(Deltaj, Hda2, si);
    pc_map(DeltagT, client.g1, Deltaj);
    gt_mul(client.K2_bT, client.K2_bT, DeltagT);

        // Server updates
        // Server update (1): W[index] <- W[index] * delta
    size_t jloc(0), kloc(index);
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

#define VESPO_G1HM_SMALL 32

void g1_mul_sim_lot_par(g1_t& RES, const g1_t* P, const bn_t* K, int N,
                        const bn_t& mod, const int nbtasks) {
#ifdef VESPO_SUB_TIMINGS
    Chrono c_step; c_step.start();
#endif

    if (nbtasks <= 1 || N < VESPO_G1HM_SMALL || N < nbtasks) {
        g1_mul_sim_lot(RES,P,K,N);
    } else {
        int64_t taskfactor(nbtasks<<1), sizeloop;
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
#pragma omp parallel for
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
    std::clog << "    Group 1 MulSimP: " << c_step.stop() << " (" << N << " operations)" << std::endl;
#endif
}

void g1_horner_mxp_seq(Polynomial<g1_t>& U, const Polynomial<g1_t>& S,
                       const bn_t& r, const int64_t from, const int64_t length) {
        for (int64_t i = 1; i < length; ++i) {
            g1_mul(U[from+i], U[from+i-1], r);			// t^r
            g1_add(U[from+i], U[from+i], S[from+i]);	// S_{i-1} t^r
        }
}

void g1_horner_mxp_iter(Polynomial<g1_t>& U, const Polynomial<g1_t>& S,
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

    if (nbtasks<=1 || length <= VESPO_G1HM_SMALL || length <= nbtasks) {
        g1_horner_mxp_seq(U,S,r,from,length);
    } else {

        for(int64_t i=0; i<(nbtasks-1); ++i) {
            g1_t tmp; g1_null(tmp); g1_new(tmp);

            g1_mul_sim_lot_par(tmp, &(S[from+i*step+1]), &(revpow[1]), step, mod, nbtasks);

            g1_mul(U[from+(i+1)*step], U[from+i*step], revpow[0]);
            g1_add(U[from+(i+1)*step], U[from+(i+1)*step], tmp);

        }


#pragma omp parallel for
        for(int64_t i=0; i<nbtasks; ++i) {
            const int64_t ipp(i*step);
            g1_horner_mxp_seq(U, S, r, from+ipp, std::min(step,length-ipp));
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

void pairing_sim(gt_t& xi_b, const g1_t* H, const g2_t* T, const int64_t deg) {
#ifdef DEBUG
    std::clog << "[PairiSim] BEG: d°" << deg << std::endl;
#endif
    gt_set_unity(xi_b);
    if (deg <= VESPO_RELIC_LIMIT_MAX_ALLOC) {
        VESPO_PC_MAP_SIM(xi_b, H, T, deg);
    } else {
        const size_t endlcm = ( (deg-1)/VESPO_RELIC_LIMIT_MAX_ALLOC)*VESPO_RELIC_LIMIT_MAX_ALLOC;
        gt_t txi_b; gt_new(txi_b);

        size_t i=0; for( ; i<endlcm; i+=VESPO_RELIC_LIMIT_MAX_ALLOC) {
            VESPO_PC_MAP_SIM(txi_b, &(H[i]), &(T[i]), VESPO_RELIC_LIMIT_MAX_ALLOC);
            gt_mul(xi_b, xi_b, txi_b);
        }
        VESPO_PC_MAP_SIM(txi_b, &(H[i]), &(T[i]), deg-endlcm);
        gt_mul(xi_b, xi_b, txi_b);

        gt_free(txi_b);
    }
#ifdef DEBUG
    std::clog << "[PairiSim] END" << std::endl;
#endif
}

void pairing_double(gt_t& xi_b, const int64_t deg,
                    const int64_t nbblocks, const bn_t& mod,
                    const Polynomial<g1_t>& H, const Polynomial<g2_t>& Ti) {
        // Server xi computation
        //   second step: applying the pairing maps as a dotproduct
#ifdef VESPO_CHECKERS
    std::clog << "[PairDblH] BEG: d°" << deg << std::endl;
#endif
#ifdef VESPO_SUB_TIMINGS
    Chrono c_cho; c_cho.start();
#endif

    int64_t dnbblocks(std::min(nbblocks, deg));
    const int64_t sizeloop( (int64_t)(std::ceil((double)deg/(double(dnbblocks)))));
    dnbblocks = (int64_t)(std::ceil((double)deg/(double(sizeloop))));
    Polynomial<gt_t> Txi(dnbblocks-1,mod);
#pragma omp parallel for
    for(int64_t i=0; i<dnbblocks; ++i) {
        const int64_t ipp(i*sizeloop);
        pairing_sim(Txi[i], &(H[ipp]), &(Ti[ipp]),
                    std::min(static_cast<int64_t>(sizeloop), deg-ipp));
    }

    gt_set_unity(xi_b);
    for(int64_t i=0; i<dnbblocks; ++i)
        gt_mul(xi_b, xi_b, Txi[i]);

//     pairing_sim(xi_b, &(H[0]), &(Ti[0]), deg);

#ifdef VESPO_SUB_TIMINGS
    std::clog << "    Pairings prods.: " << c_cho.stop() << " (" << deg << " operations)" << std::endl;
#endif
#ifdef VESPO_CHECKERS
    std::clog << "[PairDblH] END" << std::endl;
#endif
}


//========================================================
// VeSPo: Audit
//========================================================
bool eval(paillier_plaintext_t& z, const client_t& client,
          const server_t& server,
          const bn_t& r, const int64_t nb_tasks,
          double& time_c, double& time_s
#ifdef VESPO_TIMINGS
          , double& time_e, double& time_p,
          double& time_ce, double& time_cg, double& time_cp, double &time_cc
#endif
          ) {

#ifdef VESPO_TIMINGS
    std::clog << "[Audit " << r << "] BEG" << std::endl;
#else
    double time_e(0.0), time_p(0.0);
#endif
    const int64_t nbblocks(server.W.size());
    const int64_t dblocks(nbblocks-1);
    double time_r(0.0);

    bn_t group_mod; bn_new(group_mod); pc_get_ord(group_mod);

    paillier_ciphertext_t zeta;
    Polynomial<paillier_ciphertext_t> bzeta(dblocks, group_mod);

    gt_t sxi1_b, sxi2_b; gt_new(sxi1_b); gt_new(sxi2_b);
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
        bn_null(revpow[i]); bn_new(revpow[i]);
        bn_copy(revpow[i],pows_r[step-i]);
    }

        // server computation : compute xi
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

        // server computation : compute zeta by blocks
#pragma omp parallel for
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


#ifdef VESPO_TIMINGS
    time_cc = c_step.stop();
    c_step.start();
#endif

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


#ifdef VESPO_TIMINGS
    time_ce = c_step.stop();
    std::clog << "  CLIENT simpow gt : " << time_ce << std::endl;
    std::clog << "  CLIENT mult      : " << time_cc << std::endl;
#endif

    time_c = c_client.stop();

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

    bn_free(tmpz);
    gt_free(sxi1_b); gt_free(sxi2_b);
    gt_free(verif1); gt_free(verif2);
    bn_free(group_mod);
#ifdef VESPO_TIMINGS
    std::clog << "[Audit (Eval)] END" << std::endl;
#endif
    return pass;
}

//========================================================
// Main
//========================================================
int main(int argc, char * argv[]) {
    //argv[1]: size of the vector (polynomial degree + 1)
	//argv[2]: modulus size in bits for Paillier
    //argv[3]: number of experiments
    //argv[4]: number of parallel tasks (default is number of threads)

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <Vector dimension> <Paillier bitsize> <Verification iterations> [tasks (default=threads)]" << std::endl;
        return 1;
    }

    std::srand(std::time(nullptr));
    int numthreads(1);

#if defined(_OPENMP)
#pragma omp parallel
#pragma omp single
{
    numthreads = omp_get_num_threads();
}
#endif

    const int64_t degree = atoi(argv[1]) - 1;	// Polynomial degree
    const uint64_t pailliersize = atoi(argv[2]);// Paillier bitsize
    const uint64_t nb_iter = atoi(argv[3]);		// Verification iterations
    const int nb_tasks = (argc>4?atoi(argv[4]):std::min(numthreads,(int)degree));	// tasks

        /********************************************************************
         * RELIC SETUP
         *******************************************************************/
    if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}
	conf_print();

    pc_core_init();

	if (pc_param_set_any() != RLC_OK) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}
	pc_param_print();

    const uint64_t group_bits = pc_param_level()<<1;
	std::clog << "  security level: " << pc_param_level() << std::endl;

        /********************************************************************
         * VESPo: parameters
         *******************************************************************/
	util_banner("VESPO:", 0);

    std::clog << "  MAX_ALLOC: " << VESPO_RELIC_LIMIT_MAX_ALLOC << std::endl;
    std::clog << "  SHAMIR_WD: " << BN_XPWDT << std::endl;
    std::clog << "  CLOCK TYP: " << CLOCKTYPE << std::endl;
	std::clog << "  OMP CORES: " << numthreads << std::endl;
	std::clog << "  NUM ITERS: " << nb_iter << std::endl;
	std::clog << "  NUM TASKS: " << nb_tasks << std::endl;

#ifdef VESPO_NOTSECURE
    std::clog << "  REUSEDPRD: " << VESPO_NOTSECURE << std::endl;
#endif

    std::vector<double> time_c(nb_iter), time_s(nb_iter), time_u((nb_iter*(nb_iter+1))/2);
    double time_i=0., time_eval=0.;

        /********************************************************************
         * VESPo: group order
         *******************************************************************/
    bn_t group_mod; bn_null(group_mod); bn_new(group_mod); pc_get_ord(group_mod);

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
        // For now 2 threads for xi1 and x2 the rest (at least 2) for Paillier
    server_t server(degree, std::max(2,nb_tasks), group_mod);

    setup(client, server, pailliersize, P, time_eval);

    time_i = c_setup.stop();

#ifdef VESPO_TIMINGS
    std::vector<double> time_e(nb_iter), time_p(nb_iter), time_ce(nb_iter), time_cg(nb_iter), time_cp(nb_iter), time_cc(nb_iter);
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
                         , time_ce[i], time_cg[i], time_cp[i], time_cc[i]
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
            std::clog << "[PAILLIER] \033[1;32mOK\033[0m" << std::endl;
        } else {
            std::cerr << "[PAILLIER] \033[1;31m****** FAIL ******\033[0m" << std::endl;
        }
        bn_free(Pr);
#endif
    }


         /********************************************************************
         * VESPo: clean-up
         *******************************************************************/
    bn_free(r);
    bn_free(group_mod);
    pc_core_clean();
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
    std::sort(time_cc.begin(),time_cc.end());
    printf("[Eval DEVI] %lu | zeta : %f%% | xi : %f%% | C-powm : %f%% | C-gsum : %f%% | C-H_dec: %f%% | C-pairings : %f%%\n", degree+1, mediandeviation(time_e), mediandeviation(time_p), mediandeviation(time_ce), mediandeviation(time_cg), mediandeviation(time_cp), mediandeviation(time_cc));
    printf("[Eval TIME] %lu | zeta : %f | xi : %f | C-powm : %f | C-gsum : %f | C-H_dec: %f | C-pairings : %f | pure-horner : %f\n", degree+1, time_e[nb_iter/2], time_p[nb_iter/2], time_ce[nb_iter/2], time_cg[nb_iter/2], time_cp[nb_iter/2], time_cc[nb_iter/2], time_eval);
#endif

    printf("[DEVIATIO] | %lu %lu %lu | audit-client : %f%% | audit-server : %f%%\n", degree+1, group_bits, pailliersize, mediandeviation(time_c), mediandeviation(time_s));

    printf("[TIMINGS ] | %lu %lu %lu | setup : %f | audit-client : %f | audit-server : %f | pure-horner : %f | update : %f  \n=== end ===\n\n", degree+1, group_bits, pailliersize, time_i, time_c[time_c.size()/2], time_s[time_s.size()/2], time_eval, time_u[time_u.size()/2]);

	return 0;
}
