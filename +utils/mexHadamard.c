/* Hadamard Transform
   mex function to take hadamard transform

   Usage: w = hadamard(x)
   x must be a REAL VALUED COLUMN VECTOR or MATRIX
   m = size(x,1) must be a POWER OF TWO

   Notes:
   1) This implementation uses exactly m*log2(m) additions/subtractions.
   2) This is symmetric and orthogonal. To invert, apply again and
      divide by vector length.

   Written by: Peter Stobbe, Caltech
   Email: stobbe@acm.caltech.edu
   Created: August 2008
   Edits by Stephen Becker, 2009--2014

      Note: in R2008b, Matlab added "fwht" and "ifwht" (the Fast Walsh-
          Hadamart Transform and the inverse) to its Signal Processing
          Toolbox.  With the default ordering and scaling, it's not
          equivalent to this, but you can change this with the following:
          y = length(x) * fwht( x, [], 'hadamard' );
          Then y should be the same as hadamard(x) up to roundoff.
          However, it appears that this code is faster than fwht.

   Update Stephen Becker, Feb 27 2014, fix compiling issue for Mac OS X
   Update Stephen Becker, Mar  3 2014, issue error if input data is sparse
   https://github.com/stephenbeckr/SparsifiedKMeans
  
   Copyright (c) 2011, Peter Stobbe
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

       * Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.
       * Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in
         the documentation and/or other materials provided with the distribution

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.

*/

#include <stdlib.h>


/* SRB: Feb 27 2014, gcc-4.8 has problems with char16_t not being defined. 
 * This  seems to fix it
 * (and do this BEFORE including mex.h) */
/* See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56086#c4 
 (but for, e.g., Mac w/ Xcode and Clang, this fails, so test
 for gcc. more possibilities here:
  https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
 but clang defines GNUC too!
 http://nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros
 */
#ifndef NO_UCHAR
#define UCHAR_OK
#endif
#if defined(__GNUC__) && !(defined(__clang__)) && defined(UCHAR_OK)
#include <uchar.h>
#endif

#include "mex.h"

/* 
 y - output
 x - input
 m - length of vector
 */
void hadamard_apply_vector(double *y, double *x, unsigned m)
{
  unsigned bit, j, k;
  double temp;

  for (j = 0; j < m; j+=2) {
      k = j+1;
      y[j] = x[j] + x[k];
      y[k] = x[j] - x[k];
  }

  for (bit = 2; bit < m; bit <<= 1) {   
    for (j = 0; j < m; j++) {
        if( (bit & j) == 0 ) {
              k = j | bit;
              temp = y[j];
              y[j] = y[j] + y[k];
              y[k] = temp - y[k];
        }
    }
  }
}

/* 
 y - output
 x - input
 m - length of vectors (number of rows)
 n - number of vectors (number of columns)
 */
void hadamard_apply_matrix(double *y, double *x, unsigned m, unsigned n)
{
    unsigned j;
    for(j = 0; j < n; j++) {
        hadamard_apply_vector(y + j*m, x + j*m, m);
    }
}


/* check that the vector length is a power of 2,
   just using bitshifting instead of log */
void checkPowerTwo(unsigned m)
{
    /* check that it's not a degenerate 0 by 1 vector or singleton */
    if (m <= 1) {
        mexErrMsgTxt("Vector length must be greater than 1.");
    }
    /* keep dividing by two until result is odd */
    while( (m & 1) == 0 ){
        m >>= 1;
    }
    /* check that m is not a multiple of an odd number greater than 1 */
    if (m > 1) {
        mexErrMsgTxt("Vector length must be power of 2.");
    }
}


/* The gateway routine. */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double *x, *y;
  unsigned m, n;
    
  /* Check for the proper number of arguments. */
  if (nrhs != 1) {
    mexErrMsgTxt("One and only one input required; must be a column vector or matrix, with # rows a power of 2.");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  /* input size */
  m = mxGetM(prhs[0]);
  checkPowerTwo(m);
  n = mxGetN(prhs[0]);
  
  if (mxIsComplex(prhs[0])) {
    mexErrMsgTxt("Input must be real.");   
  } else if (mxIsSparse(prhs[0])) {
    mexErrMsgTxt("Input must be a full matrix, not sparse.");   
  } else if (!mxIsDouble(prhs[0])) {
    mexErrMsgTxt("Input must be of type double.");      
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  
  /* Assign pointers to each input and output. */
  x = mxGetPr(prhs[0]);
  y = mxGetPr(plhs[0]);
  
  /* Call the C subroutine. */
  hadamard_apply_matrix(y, x, m, n);
  return;
}
