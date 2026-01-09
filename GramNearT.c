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


// #####################################################################
// For the given 't', we returns the largest 'n' value where g_n <= T.
// #####################################################################
int GramNearT(mpfr_t *Result, mpfr_t T)
{
mpfr_t			Theta, Temp1;

// -------------------------------------------------------------------
// Initiate our local MPFR variables.
// -------------------------------------------------------------------
mpfr_inits2 (hgt_init.DefaultBits, Theta, Temp1, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// With \theta(g_n) = \pi * n, we are looking for the largest Gram
// Point g_n <= 'T'. Just compute \theta(T) / \pi and round down.
// -------------------------------------------------------------------
ThetaOfT(&Theta, T);
mpfr_div (Temp1, Theta, hgt_init.myPi, MPFR_RNDN);
mpfr_floor (Temp1, Temp1);		// round down to integer
mpfr_swap (*Result, Temp1);

// -------------------------------------------------------------------
// Clear our local MPFR variables.
// -------------------------------------------------------------------
mpfr_clears (Theta, Temp1, (mpfr_ptr) 0);
return(1);
}
