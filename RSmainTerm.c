// -------------------------------------------------------------------
// Program last modified December 4, 2025. 
// Copyright (c) 2024-2025 Terrence P. Murphy
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


// *******************************************************************
// We compute the main term of the Riemann-Siegel formula.  
// *******************************************************************
int RS_MainTerm(mpfr_t *Result, mpfr_t t, uint64_t N, int iFloatBits)
{
mpfr_t		Theta, Temp1, Temp2;
mpfr_t		Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne;
uint64_t	n;

// -------------------------------------------------------------------
// Check the case N < 1 (nothing to do so return 0 in Result).
// -------------------------------------------------------------------
if(N < 1)
	{
	mpfr_set_ui (*Result, 0, MPFR_RNDN);
	return(1);
	}

// -------------------------------------------------------------------
// Initialize the MPFR variables.
// -------------------------------------------------------------------
mpfr_inits2 (iFloatBits, Theta, Temp1, Temp2, 
	Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Compute Theta.
// -------------------------------------------------------------------
ThetaOfT(&Theta, t);

// -------------------------------------------------------------------
// Loop n = 1 to N.  
// For each n, we are calculating: sqrt(1/n) * cos[theta(t) - t log n]
//
// We will handle the n = 1 case separately before entering the loop.
// ------------------------------------------------------------------
// -------------------------------------------------------------------
// For the n = 1 term, we set the initial value of Main to cos(theta).
// -------------------------------------------------------------------
mpfr_cos (Main, Theta, MPFR_RNDN); 

// -------------------------------------------------------------------
// Now process the n = 2 through n = N terms
// -------------------------------------------------------------------
for (n = 2; n <= N; ++n) { 
	// ---------------------------------------------------------------
	// We need an mpfr_t version of n, so save in Temp1
	// ---------------------------------------------------------------	
	mpfr_set_uj (Temp1, n, MPFR_RNDN);
	
	// ---------------------------------------------------------------
	// First, compute the square root of 1/n
	// ---------------------------------------------------------------	
	mpfr_rec_sqrt (RecipSqrtn, Temp1, MPFR_RNDN);	

	// ---------------------------------------------------------------
	// Second, compute the argument to the cosine term.
	// That is, CosArg = [theta(t) - t log n].  
	// Then (further below) compute cos(CosArg).
	// ---------------------------------------------------------------	
	mpfr_log (Temp2, Temp1, MPFR_RNDN);			// log n
	mpfr_mul (Temp2, t, Temp2, MPFR_RNDN); 		// t * log n
	mpfr_sub (CosArg, Theta, Temp2, MPFR_RNDN); // theta(t) - [t * log n]	

	//----------------------------------------------------------------
	// We are now ready to compute the cosine value = CosCalc.
	//----------------------------------------------------------------
	mpfr_cos (CosCalc, CosArg, MPFR_RNDN);
	//----------------------------------------------------------------
	// For the full term, multiply CosCalc by RecipSqrtn, then
	// add to Main.
	//----------------------------------------------------------------
	mpfr_mul (FullTerm, RecipSqrtn, CosCalc, MPFR_RNDN);
	mpfr_add (Main, Main, FullTerm, MPFR_RNDN);
	} // end of for loop

// -------------------------------------------------------------------
// We have calculated Main.  Now, multiply by 2 and return result.
// -------------------------------------------------------------------
mpfr_mul_2ui (*Result, Main, 1, MPFR_RNDN);

// -------------------------------------------------------------------
// Free the space used by the local mpfr (constant) variables
// -------------------------------------------------------------------	
mpfr_clears (Theta, Temp1, Temp2,  
  Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne, (mpfr_ptr) 0);

return(1);
}
