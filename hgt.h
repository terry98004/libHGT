
// -------------------------------------------------------------------
// File last modified December 8, 2025. 
// Copyright (c) 2025 Terrence P. Murphy
// MIT License -- see below for details.
// -------------------------------------------------------------------

/*
MIT License

Copyright (c) 2025 Terrence P. Murphy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

struct HGT_INIT {
	mpfr_t		myPi;
	mpfr_t		my2Pi;
	mpfr_t		myLog2;
	int			DefaultBits;
	int			MaxThreads;
	int			DebugFlags;
}; 

typedef int	(*pHardyZCallback)(mpfr_t, mpfr_t, int, int);

struct computeHZ {
	mpfr_t		t; 					// 't' value to compute
	mpfr_t		Result; 			// To hold mpfr computed value
}; 


#define		HGT_PRECISION_DEFAULT	256
#define		HGT_PRECISION_MIN		64
#define		HGT_PRECISION_MAX		1024

#define		HGT_THREADS_MIN			1
#define		HGT_THREADS_MAX			8
#define		HGT_MAX_CMDLINE_STRLEN	98

#define		HGT_DEBUG_MIN			2
#define		HGT_DEBUG_MAX			30030		// for up to 23 = 223,092,870
												// for up to 19 = 9,699,690
#define		HGT_HARDY_T_MIN			1	
#define		HGT_HARDY_T_MAX			1.15e20		// (otherwise N > 8-byte uint64_t)

#define		HGT_T_INCR_MIN			1e-32	
#define		HGT_T_INCR_MAX			1.15e10

#define		HGT_COUNT_MIN			1	
#define		HGT_COUNT_MAX			9999

#define		HGT_GRAM_T_MIN			100	
#define		HGT_GRAM_T_MAX			2e30	

#define		HGT_GRAM_N_MIN			-1											
#define		HGT_GRAM_N_MAX			2e32

#define		HGT_GRAM_LOOP_MAX		1000									

#define		HGT_GRAM_ACCURACY_MIN	1
#define		HGT_GRAM_ACCURACY_MAX	32								

#define		HGT_TUR_D_POINT_MIN		2
#define		HGT_TUR_D_POINT_MAX		16

#define		HGT_TUR_GRAM_PTS_MIN	1
#define		HGT_TUR_GRAM_PTS_MAX	16

#define		HGT_TUR_SUBINTVL_MIN	8
#define		HGT_TUR_SUBINTVL_MAX	32

#define		HGT_RPT_DEC_PLACES_MIN	2	
#define		HGT_RPT_DEC_PLACES_MAX	60

#define		GABCKE_COEFF_PER_Cj		44
#define		GABCKE_NUM_Cj_TERMS		5
#define		GABCKE_DECIMAL_PLACES	50
#define		GABCKE_NUM_POWERS_P		88

#define		THETA_MAX_T_POWER3		1.1e12

// -------------------------------------------------------------------
// The last 4 debug flaga are reserved for the code that uses the
// libhgt.a library.  The first two debug flags are reserved for the 
// files making up the libhgt.a library.
// -------------------------------------------------------------------
#define		PRINT_REMAINDER			2	// used in HardyZcalc.c
#define		PRINT_COEFF				3	// used in RSbuildcoeff.c
#define 	HGT_DEBUG_RESERVED1		5
#define		HGT_DEBUG_RESERVED2		7
#define		HGT_DEBUG_RESERVED3		11
#define		HGT_DEBUG_RESERVED4		13

bool DebugMode(int iDebug, int DebugNum);
int GetSmallPositiveInteger(const char *str, int Min, int Max);
int GetLargeNumber(const char *str, __float128 Min, __float128 Max, bool IntegerOnly);
int GetDecimalDigits(const char *str);

int InitMPFR(int DefaultBits, int MaxThreads, int DebugFlags, bool CalcHardy);
int CloseMPFR(void);

int	InitCoeffMPFR(int iFloatBits);
int	CloseCoeffMPFR(void);
int	BuildCoefficientsMPFR(void);
int	CoeffStrToMPFR(mpfr_t *Result, const char *strCoeff);

int ThetaOfT(mpfr_t *Theta, mpfr_t t);
int GramAtN(mpfr_t *Result, mpfr_t N, mpfr_t Accuracy);
int GramNearT(mpfr_t *Result, mpfr_t T);
int RS_MainTerm(mpfr_t *Result, mpfr_t t, uint64_t N, int iFloatBits);
int RS_Remainder(mpfr_t *Result, mpfr_t tOver2Pi, bool nEven, mpfr_t P, int iFloatBits);

int HardyZWithCount(mpfr_t t, mpfr_t Incr, int Count, int CallerID, pHardyZCallback pCallbackHZ);
void * HardyZSingleThreaded(void * comphz);
int HardyZSingle(struct computeHZ * comphz);

void Show128(__float128 x);

int ValidateHardyT(const char *str);
int ValidateIncr(const char *str);
int ValidateCount(const char *str);
int ValidateThreads(const char *str);
int ValidateDebugFlags(const char *str);
int ValidatePrecisionMPFR(const char *str);
int ValidateReportDecimalPlaces(const char *str);

int ValidateGramN(const char *str);
int ValidateGramT(const char *str);
int ValidateGramAccuracy(const char *str);

int ValidateTuringGramPoints(const char *str);
int ValidateTuringSubIntervals(const char *str);


