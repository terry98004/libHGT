// -------------------------------------------------------------------
// Program last modified January 14, 2026. 
// Copyright (c) 2024-2026 Terrence P. Murphy
// MIT License -- see hgt.h for details.
// -------------------------------------------------------------------

#include <stdio.h>
#include <math.h>			
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <pthread.h>
#include <mpfr.h>

#include "hgt.h"

extern struct	HGT_INIT	hgt_init;

// *******************************************************************
// We compute the Hardy Z values here.  The loop is over the count of
// different 't' values to compute (based on iCount). Inside the 
// loop, we call ComputeSingleHardyZ. We then printf the result,
// and then increase 't' by Incr and repeat iCount times.
// *******************************************************************
int HardyZWithCount(mpfr_t t, mpfr_t Incr, int Count, int CallerID, pHardyZCallback pCallbackHZ)
{
struct computeHZ 	comphz[HGT_THREADS_MAX];
mpfr_t				localT;	// to avoid overwriting the passed 't'
int					i;

// -------------------------------------------------------------------
// Initialize several MPFR variables.
// -------------------------------------------------------------------
for(i=0; i< hgt_init.MaxThreads; i++) { 
	mpfr_inits2 (hgt_init.DefaultBits, comphz[i].t, comphz[i].Result, (mpfr_ptr) 0);
	}
mpfr_init2 (localT, hgt_init.DefaultBits);
mpfr_set (localT, t, MPFR_RNDN);

// -------------------------------------------------------------------
// We loop Count times, processing 't' in the first (i = 1) loop.
// If Count is greater than one, we increment 't' by the Incr
// amount and loop again.  
//
// In the loop, we manage threads (if requested) and call HardyZSingle 
// (directly or indirectly) Count times.  We print the results using 
// pCallbackHZ.
// -------------------------------------------------------------------
int			j, m, iRemaining;
pthread_t	thread_id[HGT_THREADS_MAX];

i = 0;
while(i < Count)
	{
	iRemaining = Count - i;
	m = hgt_init.MaxThreads > iRemaining ? iRemaining : hgt_init.MaxThreads;
	
	for (j = 0; j < m; j++)
		{
		// ---------------------------------------------------------------
		// Update comphz values for the current 't'
		// ---------------------------------------------------------------	
		mpfr_set (comphz[j].t, localT, MPFR_RNDN);
		mpfr_set_ui (comphz[j].Result, 0, MPFR_RNDN);
		
		if(hgt_init.MaxThreads > 1) {
			pthread_create( &thread_id[j], NULL, HardyZSingleThreaded, &comphz[j] );
			}
		else {
			HardyZSingle(&comphz[0]);
			}
		mpfr_add (localT, localT, Incr, MPFR_RNDN);
		}
	for (j = 0; j < m; j++)
		{
		if(hgt_init.MaxThreads > 1) {
			pthread_join( thread_id[j], NULL); 
			}
		}		
	for (j = 0; j < m; j++)
		{	
		pCallbackHZ(comphz[j].t, comphz[j].Result, i, CallerID);
		i++;	// we need 'i' here so that it is incremented 'm' times 
		}			
	}	// end of outer for loop

// -------------------------------------------------------------------
// We are done.  Clear our local MPFR variables.
// -------------------------------------------------------------------
for(i=0; i< hgt_init.MaxThreads; i++) {
	mpfr_clears (comphz[i].t, comphz[i].Result, (mpfr_ptr) 0);
	}
mpfr_clear(localT);
return(1);
}

// *******************************************************************
// For correct "C type" on the pthread_create call, we need this
// pass-through function.
// *******************************************************************
void * HardyZSingleThreaded(void * comphz)
{
HardyZSingle((struct computeHZ *) comphz);
return(NULL);
}

// *******************************************************************
// We compute a single Hardy Z values here.  We use the following 
// passed variable:
// 		comphz->t and comphz->ptrResult
// and the following global variables:
// 		hgt_init.my2Pi and hgt_init.DebugFlags
// *******************************************************************
int HardyZSingle(struct computeHZ * comphz)
{
mpfr_t			tOver2Pi, T, N, P, Main, Remainder;
uint64_t		ui64N;
bool			nEven;

mpfr_inits2 (hgt_init.DefaultBits,  
				tOver2Pi, T, N, P, Main, Remainder, (mpfr_ptr) 0);

// ---------------------------------------------------------------
// Compute N and P for the given 't'. 
// NOTE: Because N is (currently) a uint64_t (we assume 64-bit),
// then 0 <= N <= 18,446,744,073,709,551,615.  
// Although this allows for a much larger 't', we will currently 
// assume that 't' does not exceed about 1.15 * 10^{20} in the 
// calcs below. (A value of 't' that was previously checked).
// ---------------------------------------------------------------	
mpfr_div (tOver2Pi, comphz->t, hgt_init.my2Pi, MPFR_RNDN);
mpfr_sqrt (T, tOver2Pi, MPFR_RNDN);
mpfr_modf (N, P, T, MPFR_RNDN);
ui64N = mpfr_get_uj (N, MPFR_RNDN);
nEven = (ui64N % 2 == 0) ? true : false;

// ---------------------------------------------------------------
// Compute the remainder term.
// ---------------------------------------------------------------		
RS_Remainder(&Remainder, tOver2Pi, nEven, P, hgt_init.DefaultBits);
	
// ---------------------------------------------------------------
// Now compute the Main term and add to Remainder to get HardyZ.
// ---------------------------------------------------------------	
RS_MainTerm(&Main, comphz->t, ui64N, hgt_init.DefaultBits);
mpfr_add (comphz->Result, Main, Remainder, MPFR_RNDN);

// -------------------------------------------------------------------
// Clear our local MPFR variables.
// -------------------------------------------------------------------
mpfr_clears (tOver2Pi, T, N, P, Main, Remainder, (mpfr_ptr) 0);
return(1);
}

