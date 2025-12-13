// -------------------------------------------------------------------
// Program last modified December 10, 2025. 
// Copyright (c) 2025 Terrence P. Murphy
// MIT License -- see hgt.h for details.
// -------------------------------------------------------------------

#include <errno.h>
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

struct	HGT_INIT	hgt_init;

// -------------------------------------------------------------------
// We call this function before using any MPFR functions.  We set the
// default MPFR precision and create global variables holding the
// values of Pi 2Pi and Log(2), 
// -------------------------------------------------------------------
int InitMPFR(int DefaultBits, int MaxThreads, int DebugFlags, bool CalcHardy)
{
hgt_init.DefaultBits 	= DefaultBits;
hgt_init.MaxThreads		= MaxThreads;
hgt_init.DebugFlags		= DebugFlags;

// -------------------------------------------------------------------
// Set default precision for MPFR
// -------------------------------------------------------------------
mpfr_set_default_prec (hgt_init.DefaultBits);

// -------------------------------------------------------------------
// Initialize the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_inits2 (hgt_init.DefaultBits, hgt_init.myPi, hgt_init.my2Pi, 
	hgt_init.myLog2, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Set the value of the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_const_pi (hgt_init.myPi, MPFR_RNDN); 
mpfr_mul_2ui (hgt_init.my2Pi, hgt_init.myPi, 1, MPFR_RNDN); /* 2Pi */
mpfr_const_log2 (hgt_init.myLog2, MPFR_RNDN);
if(CalcHardy == true){
	InitCoeffMPFR(DefaultBits);
	}
return(1);
}

// -------------------------------------------------------------------
// We call this function after we are done using all MPFR functions.  
// We clear the remaining memory and any cache used by MPFR. 
// -------------------------------------------------------------------
int CloseMPFR(void)
{
// -------------------------------------------------------------------
// Free the space used by the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_clears (hgt_init.myPi, hgt_init.my2Pi, 
			hgt_init.myLog2, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Clear the cache used by MPFR.
// -------------------------------------------------------------------
mpfr_free_cache ();
return(1);
}


// -------------------------------------------------------------------
// Validate the text string representing the ordinate 't' on the
// ctitical line.
// -------------------------------------------------------------------
int ValidateHardyT(const char *str)
{
return(GetLargeNumber(str, HGT_HARDY_T_MIN, HGT_HARDY_T_MAX, false));
}


// -------------------------------------------------------------------
// Validate the text string with the amount to increment 't'.
// -------------------------------------------------------------------
int ValidateIncr(const char *str)
{
return(GetLargeNumber(str, HGT_T_INCR_MIN, HGT_T_INCR_MAX, false));
}


// -------------------------------------------------------------------
// Validate the text string with a positive integer count value.
// -------------------------------------------------------------------
int ValidateCount(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_COUNT_MIN, HGT_COUNT_MAX));
}


// -------------------------------------------------------------------
// Validate the text string with a positive integer number of threads.
// -------------------------------------------------------------------
int ValidateThreads(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_THREADS_MIN, HGT_THREADS_MAX));
}


// -------------------------------------------------------------------
// Validate the text string with the requested debug flags.
// -------------------------------------------------------------------
int ValidateDebugFlags(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_DEBUG_MIN, HGT_DEBUG_MAX));
}


// -------------------------------------------------------------------
// Validate the text string with a positive integer number of MPFR bits.
// -------------------------------------------------------------------
int ValidatePrecisionMPFR(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_PRECISION_MIN, HGT_PRECISION_MAX));
}


// -------------------------------------------------------------------
// Validate the text string with a positive integer number of decimal
// places in a report.
// -------------------------------------------------------------------
int ValidateReportDecimalPlaces(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_RPT_DEC_PLACES_MIN, HGT_RPT_DEC_PLACES_MAX));
}


// -------------------------------------------------------------------
// We are given a validated string (numbers plus zero or one decimal
// point). We return the number of decimal digits (i.e., digits to
// the right of the decimal point).
// -------------------------------------------------------------------
int GetDecimalDigits(const char *str)
{
const char 	* ptr = strchr(str, '.'); // point to decimmal point, if any

if(!ptr) return(0); 	// no decimal point so no decimal digits
return(strlen(++ptr));	// found decimal point, count decimal digits
}


// -------------------------------------------------------------------
// Validate the text string representing the ordinate 't' on the
// ctitical line. (Larger value allowed in Gram than in Hardy calcs).
// -------------------------------------------------------------------
int ValidateGramT(const char *str)
{
return(GetLargeNumber(str, HGT_GRAM_T_MIN, HGT_GRAM_T_MAX, false));
}


// -------------------------------------------------------------------
// Validate the text string representing the Nth Gram Point.
// -------------------------------------------------------------------
int ValidateGramN(const char *str)
{
return(GetLargeNumber (str, HGT_GRAM_N_MIN, HGT_GRAM_N_MAX, true));
}


// -------------------------------------------------------------------
// Validate the text string representing the required accuracy of the
// Gram Point calculation (in decimal places).
// -------------------------------------------------------------------
int ValidateGramAccuracy(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_GRAM_ACCURACY_MIN, HGT_GRAM_ACCURACY_MAX));
}


// -------------------------------------------------------------------
// Validate the text string with a positive integer count value.
// -------------------------------------------------------------------
int ValidateTuringGramPoints(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_TUR_GRAM_PTS_MIN, HGT_TUR_GRAM_PTS_MAX));
}


// -------------------------------------------------------------------
// Validate the text string with a positive integer count value.
// -------------------------------------------------------------------
int ValidateTuringSubIntervals(const char *str)
{
return(GetSmallPositiveInteger(str, HGT_TUR_SUBINTVL_MIN, HGT_TUR_SUBINTVL_MAX));
}

// -------------------------------------------------------------------
// We check whether the value passed on the command line (using the
// -d debug parameter, and saved in hgt_init.DebugFlags "matches" the 
// DebugNum parameter passed in here.  There is a "match" if there is 
// no remainder when you divide DebugFlags by DebugNum.
// -------------------------------------------------------------------
bool DebugMode(int DebugFlags, int DebugNum)
{
return(DebugFlags % DebugNum == 0 ? true : false);
}

// -------------------------------------------------------------------
// We are given: (1) a string that must be validated as a large number,
// and (2) the Min or Max allowed values for the number.
// If validated, we return1, otherwise we return a negative value.
// As a side note, __float128 are accurate to about 34 digits.
// -------------------------------------------------------------------
int GetLargeNumber(const char *str, __float128 Min, __float128 Max, bool IntegerOnly)
{
char * 		endptr;
size_t		Len;
__float128 	Value;
int			Result = 1;
const char  sAllowed[] = "-.0123456789"; 

Len = strlen(str);
Value = strtoflt128 (str, &endptr);

if(Len < 1 || Len > HGT_MAX_CMDLINE_STRLEN || strspn(str, sAllowed) != Len){
	Result = -1; 
	}
else if (endptr == str){
	Result = -2;			// No digits were found
	} 
else if (*endptr != '\0'){
	Result = -3;			// Invalid character
	}
else if (Value < Min || Value > Max){
	Result = -4; 
	}
else if (IntegerOnly == true && strchr(str, '.') != NULL){
	Result = -5; 
	}	
else if(errno != 0) {
	Result = -6;
	}
return(Result);
}

// -------------------------------------------------------------------
// We are given: (1) a string that must be validated as a small positive
// integer, and (2) the Min or Max allowed values for the integer.
// If validated, we return the integer value; otherwise we return Min - 1.
// -------------------------------------------------------------------
int GetSmallPositiveInteger(const char *str, int Min, int Max)
{
char * 		ptr;
size_t		Len;
int			Value;
const char  sNum[] = "0123456789"; 

Len = strlen(str);

if(Len < 1 || Len > 6 || strspn(str, sNum) != Len){
	return(-1); 
	}

Value = (int) strtol(str, &ptr, 10);
if(Value < Min || Value > Max){
	return(-1); 
	}
return(Value);
}

// -------------------------------------------------------------------
// For testing, this does a printf of a __float128 value.
// -------------------------------------------------------------------
void Show128(__float128 x)
{
char buf[128];

quadmath_snprintf (buf, sizeof buf, "%.20Qf", x);
printf("%s \n", buf);	
}
