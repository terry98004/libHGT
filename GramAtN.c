// -------------------------------------------------------------------
// Program last modified January 9, 2026. 
// Copyright (c) 2025-2026 Terrence P. Murphy
// MIT License -- see hgt.h for details.
// -------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <mpfr.h>

#include "hgt.h"

extern struct	HGT_INIT	hgt_init;


const char 	chGramSmall[4][50] = {
	{"3.4360829361"},
	{"17.8455995404108608168263384125190970356932874"},
	{"23.1702827012463092789966435383015320517470983"}, 
	{"27.6701822178483449617890"} };
	
// #####################################################################
// For our given N, compute the Gram number.
// #####################################################################
int GramAtN(mpfr_t *Result, mpfr_t N, mpfr_t Accuracy)
{
mpfr_t		tAbove, tBelow;
mpfr_t		LogN, NOverLogN, k, x, xN, Temp1, Temp2;
mpfr_t		tNewBelow, tNewAbove, tMid, tGap, tHalfGap, thetaMid, thetaDelta, nPi;
long int	i;
bool		bFinished = false;

// -------------------------------------------------------------------
// initialize all mpfr_t variables used in computing the Gram point
// -------------------------------------------------------------------
mpfr_inits2 (hgt_init.DefaultBits, tAbove, tBelow, (mpfr_ptr) 0);
mpfr_inits2 (hgt_init.DefaultBits, LogN, NOverLogN, k, x, xN, Temp1, Temp2, (mpfr_ptr) 0);
mpfr_inits2 (hgt_init.DefaultBits, tNewBelow, tNewAbove, tMid, tGap, tHalfGap, 
	thetaMid, thetaDelta, nPi, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Step 1. For n = -1, 0, 1, 2 we use a lookup table (to avoid "edge" 
// conditions with those small numbers) to set their final Gram number 
// value.
// Otherwise, we initialize tBelow and tAbove with values known to
// ensure tBelow < g_n < tAbove. 
// -------------------------------------------------------------------
if( mpfr_cmp_si (N, 3.0) < 0 ) {
	i = mpfr_get_si (N, MPFR_RNDN);
	mpfr_set_str (tMid, chGramSmall[i+1], 10, MPFR_RNDN);
	bFinished = true;
	}
else {
	// ---------------------------------------------------------------
	// Set tBelow to max (30, n / log n).
	// ---------------------------------------------------------------
	mpfr_log (LogN, N, MPFR_RNDN);
	mpfr_div (NOverLogN, N, LogN, MPFR_RNDN);
	 mpfr_set_ui (Temp1, 30, MPFR_RNDN);
	mpfr_max (tBelow, Temp1, NOverLogN, MPFR_RNDN);

	// ---------------------------------------------------------------
	// Set tAbove (a bit more complicated formula).  We
	//
	// Set Temp1 = 6.64 / (LogN - 2.85)
	// Set k = max (0.5, Temp1)
	// Set x = min (3.0, k)
	// Set xN = x * N
	// Set tAbove = max (xN, 250)
	// ---------------------------------------------------------------	
	mpfr_sub_d (Temp1, LogN, 2.85, MPFR_RNDN);
	mpfr_d_div (Temp1, 6.64, Temp1, MPFR_RNDN); 
	
	mpfr_set_d (Temp2, 0.5, MPFR_RNDN);
	mpfr_max (k, Temp1, Temp2, MPFR_RNDN);
	
	mpfr_set_d (Temp2, 3.0, MPFR_RNDN);
	mpfr_min (x, Temp2, k, MPFR_RNDN);
	
	mpfr_mul (xN, x, N, MPFR_RNDN);
	mpfr_set_d (Temp2, 250.0, MPFR_RNDN);
	mpfr_max (tAbove, xN, Temp2, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// Step 2. Final calculation of the Gram point.
// -------------------------------------------------------------------
	// ---------------------------------------------------------------
	// Set the starting values for tNewBelow and tNewAbove, plus the 
	// value of nPi.
	// ---------------------------------------------------------------
if(bFinished == false) {
	mpfr_set (tNewBelow, tBelow, MPFR_RNDN);
	mpfr_set (tNewAbove, tAbove, MPFR_RNDN);
	mpfr_mul (nPi, hgt_init.myPi, N, MPFR_RNDN);
	}

	// ---------------------------------------------------------------
	// Loop up to GRAM_LOOP_MAX number of times.  In each case, we find 
	// the midpoint tMid between tNewBelow and tNewAbove, then compute 
	// \theta(tMid) and compare that result to nPi.  Depending on the 
	// result, we remove the top half or bottom half of the interval 
	// [setting either tNewBelow or tNewAbove to tMid] and loop again.  
	// We stop the loop if \theta(tMid) is sufficiently close to nPi.
	// ---------------------------------------------------------------
for(i = 0; i < HGT_GRAM_LOOP_MAX && bFinished == false; i++) {	
	// ---------------------------------------------------------------
	// Set tMid to halfway between tAbove and tBelow.
	// ---------------------------------------------------------------
	mpfr_sub (tGap, tNewAbove, tNewBelow, MPFR_RNDN);
	mpfr_div_ui (tHalfGap, tGap, 2, MPFR_RNDN);
	mpfr_add (tMid, tHalfGap, tNewBelow, MPFR_RNDN);
	
	// ---------------------------------------------------------------
	// Next compute \theta(tMid).
	// Then compute thetaDelta = [\theta(tMid) - nPi].  
	// We want to know two things: (1) is thetaDelta small enough so 
	// that we are done? and (2) if not, is thetaDelta positive (so we 
	// can set tAbove to tMid) or negative (tBelow = tMid).
	// ---------------------------------------------------------------	
	ThetaOfT(&thetaMid, tMid);
	mpfr_sub (thetaDelta, thetaMid, nPi, MPFR_RNDN);
	mpfr_abs (Temp1, thetaDelta, MPFR_RNDN);

	if( mpfr_cmp (Temp1, Accuracy) < 0){
		bFinished = true;
		}
	else if (mpfr_sgn (thetaDelta) > 0) {
		mpfr_set (tNewAbove, tMid, MPFR_RNDN);
		}
	else {
		mpfr_set (tNewBelow, tMid, MPFR_RNDN);
		}
	} // end of for loop 

	// ---------------------------------------------------------------
	// Did we fail to find the required tMid? If so, display an error
	// message.  In all cases, save the last tMid value in Result.
	// ---------------------------------------------------------------
if(bFinished == false) {
	mpfr_printf("FAILURE!!!: Gram = %.20Rf, (Gram - n * pi) = %.20Rf \n", tMid, thetaDelta);
	}
mpfr_swap (*Result, tMid);

// -------------------------------------------------------------------
// Clear our local MPFR variables.
// -------------------------------------------------------------------
mpfr_clears (tAbove, tBelow, (mpfr_ptr) 0);
mpfr_clears (LogN, NOverLogN, k, x, xN, Temp1, Temp2, (mpfr_ptr) 0);
mpfr_clears (tNewBelow, tNewAbove, tMid, tGap, tHalfGap, 
	thetaMid, thetaDelta, nPi, (mpfr_ptr) 0);
return(1);
}
