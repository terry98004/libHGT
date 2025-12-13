// -------------------------------------------------------------------
// Program last modified October 7, 2025. 
// Copyright (c) 2025 Terrence P. Murphy
// MIT License -- see hgt.h for details.
// -------------------------------------------------------------------

#include <quadmath.h>
#define MPFR_WANT_FLOAT128 1
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <mpfr.h>

#include "hgt.h"

extern struct	HGT_INIT	hgt_init;

// -------------------------------------------------------------------
// We compute \theta(t) as used in the Riemann-Siegel formula.  
// The formula from our book is:
//	
// Theta = ((t/2) * (log(t / (2 * PI)))) - PI/8 - t/2
//		+ 1/(48 * t) + 7/(5760 * powl(t, 3));
//
// Using the variable names below, the formula is:
//    Theta = tOver2 * [LogOftOver2Pi] - PiOver8 - tOver2
//         + Recip48t + Power3Term
// 
// Or, equally (as implemented below):
//    Theta = tOver2 * [LogOftOver2Pi - 1] 
//         + Recip48t - PiOver8 + Power3Term
// -------------------------------------------------------------------
int ThetaOfT(mpfr_t *Theta, mpfr_t t)
{
mpfr_t		tOver2, PiOver8, LogOftOver2Pi;
mpfr_t		Recip48t, Power3Term, Temp1, MinorTerms;

// -------------------------------------------------------------------
// initialize all mpfr_t variables used in computing Theta
// -------------------------------------------------------------------
mpfr_inits2 (hgt_init.DefaultBits, tOver2, PiOver8, LogOftOver2Pi, 
	Recip48t, Power3Term, Temp1, MinorTerms, (mpfr_ptr) 0);

// set tOver2
mpfr_div_ui (tOver2, t, 2, MPFR_RNDN);

// set LogOftOver2Pi
mpfr_div (Temp1, tOver2, hgt_init.myPi, MPFR_RNDN);
mpfr_log (LogOftOver2Pi, Temp1, MPFR_RNDN);

// set PiOver8
mpfr_div_ui (PiOver8, hgt_init.myPi, 8,  MPFR_RNDN);

// set Recip48t
mpfr_mul_ui (Temp1, t, 48, MPFR_RNDN);
mpfr_ui_div (Recip48t, 1, Temp1, MPFR_RNDN);

// -------------------------------------------------------------------
// Compute the minor terms in the \theta(t) formula:
//    MinorTerms = Recip48t - PiOver8 + Power3Term
// -------------------------------------------------------------------
mpfr_set (MinorTerms, Recip48t, MPFR_RNDN);
mpfr_sub (MinorTerms, MinorTerms, PiOver8, MPFR_RNDN);

// -------------------------------------------------------------------
// Calculate and add the Powers3Term UNLESS t is so large that the 
// computed value of this term will be too small to matter.
// -------------------------------------------------------------------
if(mpfr_cmp_d (t, THETA_MAX_T_POWER3) < 0)
	{
	mpfr_pow_si (Temp1, t, -3, MPFR_RNDN);
	mpfr_mul_ui (Temp1, Temp1, 7, MPFR_RNDN);
	mpfr_div_ui (Power3Term, Temp1, 5760, MPFR_RNDN);
	mpfr_add (MinorTerms, MinorTerms, Power3Term, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// Now calculate the major term = tOver2 * [LogOftOver2Pi - 1]
// -------------------------------------------------------------------
mpfr_sub_ui (Temp1, LogOftOver2Pi, 1, MPFR_RNDN);
mpfr_mul (Temp1, tOver2, Temp1, MPFR_RNDN);
mpfr_add (*Theta, Temp1, MinorTerms, MPFR_RNDN);	

// -------------------------------------------------------------------
// Free the space used by the local mpfr (constant) variables
// -------------------------------------------------------------------	
mpfr_clears ( tOver2, PiOver8, LogOftOver2Pi, 
	Recip48t, Power3Term, Temp1, MinorTerms, (mpfr_ptr) 0);
return(1);
}