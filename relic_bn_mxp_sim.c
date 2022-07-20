/*
 * Copyright (c) 2022 Jean-Guillaume Dumas
 *
 */

/**
 * @file
 *
 * Implementation of the simultaneous multiple precision exponentiation.
 *
 * Available also in, e.g., https://github.com/jgdumas/relic/tree/bnMXPsim
 *
 * @ingroup bn
 */


/****************************************************************
 * bn_t simultaneous exponentiation generalized Shamir trick
 * utility subroutine overwritting u, with precomputed table t
 ****************************************************************/
void _bn_mxp_sim(bn_t S, const bn_t P[BN_XPWDT], bn_t u[BN_XPWDT],
                 const bn_t T[BN_XPWDT], const bn_t mod) {
        // WARNING: overwrites u
    int iszeroexp = 0x1;
    for(unsigned int j=0; j<BN_XPWDT; ++j) {
        iszeroexp &= bn_is_zero(u[j]);
    }

    if (iszeroexp) {
        bn_set_dig(S,1); // all exponents zero, just return 1
        return;
    }

        // Select odd exponents
    unsigned int parities = !bn_is_even(u[0]);
    for(unsigned int j=1; j<BN_XPWDT; ++j) parities |= (unsigned int)(!bn_is_even(u[j]))<<j;

        // Halving exponents
    for(unsigned int j=0; j<BN_XPWDT; ++j) {
        bn_hlv(u[j],u[j]);			// WARNING: u overwriten
    }

        // Recursive Power up to the halves
    _bn_mxp_sim(S, P, u, T, mod);

        // One Squaring
    bn_sqr(S,S); bn_mod(S,S,mod);

        // One multiplication by the odd exponents
    if (parities) {
        bn_mul(S, S, T[parities]);
        bn_mod(S,S,mod);
    }
}

#define BN_XPWDT_TABLE_SIZE (1u<<BN_XPWDT)

/****************************************************************
 * bn_t simultaneous exponentiation generalized Shamir trick and fixed width
 ****************************************************************/
void bn_mxp_sim(bn_t S, const bn_t P[BN_XPWDT], const bn_t u[BN_XPWDT], const bn_t mod) {
    bn_t hu[BN_XPWDT];
    bn_t T[BN_XPWDT_TABLE_SIZE];

    RLC_TRY {

        bn_null(T[0]); bn_new(T[0]); bn_set_dig(T[0],1);

            // Precompute all 2^{BN_XPWDT} combinations of points P
        for(unsigned int i=0; i<BN_XPWDT; ++i) {
            const unsigned int star = 1<<i; const unsigned int stars = star<<1;
            if (! bn_is_zero(u[i])) { // Otherwise will never need P[i]
                bn_null(T[star]); bn_new(T[star]); bn_copy(T[star], P[i]);
                for(unsigned int j=star+1; j<stars; ++j) {
                    bn_null(T[j]); bn_new(T[j]);
                    bn_mul(T[j], T[star], T[j-star]); bn_mod(T[j],T[j],mod);
                }
            }
        }

            // copy u, as hu will be overwritten by the subroutine
        for(unsigned int j=0; j<BN_XPWDT; ++j) {
            bn_null(hu[j]); bn_new(hu[j]); bn_copy(hu[j],u[j]);
        }

            // Call the exponentiation subroutine
        _bn_mxp_sim(S,P,hu,T,mod);

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
        for(unsigned int i=0; i<BN_XPWDT_TABLE_SIZE; ++i) { bn_free(T[i]); }
        for(unsigned int j=0; j<BN_XPWDT; ++j) { bn_free(hu[j]); }
    }
}

/****************************************************************
 * bn_t simultaneous exponentiation generalized Shamir trick any width
 ****************************************************************/
void bn_mxp_sim_lot(bn_t S, const bn_t P[], const bn_t u[], const bn_t mod, int n) {
    bn_t wP[BN_XPWDT], wu[BN_XPWDT], tmp;
	RLC_TRY {
            // Will use blocks of size BN_XPWDT
        bn_null(tmp); bn_new(tmp);
        for(unsigned int j=0; j<BN_XPWDT; ++j) {
            bn_null(wP[j]); bn_new(wP[j]);
            bn_null(wu[j]); bn_new(wu[j]);
        }

            // Largest multiple of BN_XPWDT lower than n
        const int endblockingloop = ( (n/BN_XPWDT)*BN_XPWDT );
        bn_set_dig(S, 1);
            // Exponentiate by blocks of size BN_XPWDT
        int i = 0; for(; i<endblockingloop; ) {
            for(unsigned int j=0; j<BN_XPWDT; ++j, ++i) {
                bn_copy(wP[j], P[i]);
                bn_copy(wu[j], u[i]);
            }
            bn_mxp_sim(tmp, wP, wu, mod);
            bn_mul(S, S, tmp);
            bn_mod(S, S, mod);
        }

            // Remaining (n-endblockingloop) exponentiations
        const int r=n-i;
        if (r) {
            if (r>1) {
                unsigned int j=0; for(; i<n; ++j, ++i) {
                    bn_copy(wP[j], P[i]);
                    bn_copy(wu[j], u[i]);
                }
                for(; j<BN_XPWDT; ++j) {	// Set remaining to exponent zero
                    bn_set_dig(wu[j], 0);
                }
                bn_mxp_sim(tmp, wP, wu, mod);
            } else { // A single exponent
                bn_mxp(tmp, P[i], u[i], mod);
            }
            bn_mul(S, S, tmp);
            bn_mod(S, S, mod);
        }

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
        for(unsigned int j=0; j<BN_XPWDT; ++j) {
            bn_free(wP[j]); bn_free(wu[j]);
        }
        bn_free(tmp);
    }
}
