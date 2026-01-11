# libHGT
A library of functions that Calculate Values of the Hardy Z Function and locates Gram Points.

## Overview

Let 't' be a positive real number and the ordinate of a point (1/2 + it) along the critical line. The **libhgt.a** library is a static library, written in C, providing functions to: (1) calculate *Z(t)*, the Hardy Z function at 't' (using the Riemann-Siegel Formula), (2) locate the Gram Point for a given interger, (3) find the largest Gram Point less than or equal to a given 't', and (4) validate user input for parameters that will be passed to the calculating functions.  The library functions also support the use of Turing's Method for verifying the Riemann Hypothesis up to a given point 't' on the critical line.  For the above reasons, the HGT in the library name stands for Hardy-Gram-Turing.

The included PDF provides important additional information, and briefly describes the mathematics behind the calculations.  Much greater detail on the mathematics is given in my book: *A Study of Riemann's Zeta Function* by Terrence P. Murphy.  The book is available on Amazon.  Go [here][website-link] to view the Table of Contents and Preface for all of my books.

## Building the Static Library

For Windows 11 users, the static library **libhgt.a** is included with any release posted on GitHub.

For other operating systems, you will need to build the library, as follows.

*  You need the [**gcc**][gcc-gnu-link] C compiler installed on your system. That installation must include the **GNU GMP** and **GNU MPFR** (floating point) libraries -- you need **mpfr.h** here and you will need the **libmpfr.a** and **libgmp.a** static libraries when you link to **libhgt.a** to create an execuable. 

*  Following the build logic in the **makehgt.bat** file, we provide a **makefile**, in the form that should work with your operating system and the **gcc** compiler.

You can then build the static library **libhgt.a** from the provided source files.

## Files

This distribution consists of the following files:

  * [README.md][readme-link]. The file you are currently reading.
    
  * [libHGT.pdf][libHGT-pdf-link]. A PDF file with further discussion of our software, the building
  of our software program, and the mathematics behind our software. The PDF is created using LaTex, 
  which was needed to allow proper display / layout of the mathematics.
  
  * [hgtInit.c][hgtInit-c-link]. This source code file contains code to initialize the MPFR floating point system, and code to validate user input for parameters that will be used by the library functions.

  * [RSbuildcoeff.c][RSbuildcoeff-c-link]. This source code file builds an **MPFR** version of the Gabcke power series coefficients as part of the overall task of initializing the **MPFR** floating point system.

  * [RSmainTerm.c][RSmainTerm-c-link]. This source code file computes the main term of the Riemann-Siegel formula.

  * [RSremainder.c][RSremainder-c-link]. This source code file computes the remainder term of the Riemann-Siegel formula.

  * [ThetaOfT.c][ThetaOfT-c-link]. This source code file computes the theta value of the passed positive ordinate T.  That computed value is a factor in the main term of the Riemann-Siegel formula.
 
  * [GramAtN.c][GramAtN-c-link]. This source code file computes the Gram Point associated with the positive integer N. 
  
  * [GramNearT.c][GramNearT-c-link]. This source code file computes the positive integer N associated with the largest Gram Point less than or equal to the positive ordinate T.

  * [HardyZcalc.c][HardyZcalc-c-link]. This source code file contains the public facing library function used to compute one or more Hardy Z values.

  * [hgt.h][hgt-h-link]. The is the only (local) include file for the library.
  
  * [makefile][makefile-link]. This makefile is for use with the make program that is available with most development environments.
  
  * [makehgt.bat][makehgt-bat-link]. The is the "makefile" for the program.  Currently,
  this file is a Windows batch file (**not** an actual makefile), but can be easily converted to 
  a standard makefile.

## Examples Using the libHGT Library

See [HardyZ][HardyZ-link], [Gram][Gram-link] and [Turing][Turing-link] for sample code using this **libHGT** library.

## Terms of use

This **Hardy Z Function Calculator** is free and distributed under the
**MIT License** (MIT). 

We used the [**gcc**][gcc-gnu-link] compiler and the [**MPFR**][mpfr-link] floating point library.
We also used the [**msys2**][msys2-link] software distribution and building platform for windows.
See their respective links for theirs terms of license.  

[website-link]:			https://riemann1859.com
[license-link]:			https://github.com/terry98004/libHGT/blob/master/license.txt
[readme-link]:			https://github.com/terry98004/libHGT/blob/master/README.md
[libHGT-pdf-link]:		https://github.com/terry98004/libHGT/blob/master/libHGT.pdf
[hgtInit-c-link]:		https://github.com/terry98004/libHGT/blob/master/hgtInit.c
[RSbuildcoeff-c-link]:	https://github.com/terry98004/libHGT/blob/master/RSbuildcoeff.c
[RSmainTerm-c-link]:	https://github.com/terry98004/libHGT/blob/master/RSmainTerm.c
[RSremainder-c-link]:	https://github.com/terry98004/libHGT/blob/master/RSremainder.c
[ThetaOfT-c-link]:		https://github.com/terry98004/libHGT/blob/master/ThetaOfT.c
[GramAtN-c-link]:		https://github.com/terry98004/libHGT/blob/master/GramAtN.c
[GramNearT-c-link]:		https://github.com/terry98004/libHGT/blob/master/GramNearT.c
[HardyZcalc-c-link]:	https://github.com/terry98004/libHGT/blob/master/HardyZcalc.c
[hgt-h-link]:			https://github.com/terry98004/libHGT/blob/master/hgt.h
[makefile-link]:	https://github.com/terry98004/libHGT/blob/master/makefile
[makehgt-bat-link]:		https://github.com/terry98004/libHGT/blob/master/makehgt.bat
[HardyZ-link]:		https://github.com/terry98004/HardyZ/
[Gram-link]:		https://github.com/terry98004/Gram/
[Turing-link]:		https://github.com/terry98004/Turing/
[mpfr-link]:			https://www.mpfr.org/
[gcc-gnu-link]:			https://gcc.gnu.org/
[msys2-link]:			https://www.msys2.org/
