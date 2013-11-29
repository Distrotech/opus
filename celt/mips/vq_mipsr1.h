/* Copyright (c) 2007-2008 CSIRO
   Copyright (c) 2007-2009 Xiph.Org Foundation
   Written by Jean-Marc Valin */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef __VQ_MIPSR1_H__
#define __VQ_MIPSR1_H__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mathops.h"
#include "cwrs.h"
#include "vq.h"
#include "arch.h"
#include "os_support.h"
#include "bands.h"
#include "rate.h"

static unsigned extract_collapse_mask(int *iy, int N, int B);
static void normalise_residual(int * OPUS_RESTRICT iy, celt_norm * OPUS_RESTRICT X, int N, opus_val32 Ryy, opus_val16 gain);
static void exp_rotation(celt_norm *X, int len, int dir, int stride, int K, int spread);

#define OVERRIDE_vq_exp_rotation1
static void exp_rotation1(celt_norm *X, int len, int stride, opus_val16 c, opus_val16 s)
{
   int i;
   celt_norm *Xptr;
   Xptr = X;
   for (i=0;i<len-stride;i++)
   {
      celt_norm x1, x2;

      x1 = Xptr[0];
      x2 = Xptr[stride];
      Xptr[stride] = EXTRACT16(MULT16_16_Q15_ADD(c,x2,s,x1));
      *Xptr++      = EXTRACT16(MULT16_16_Q15_SUB(c,x1,s,x2));
   }
   Xptr = &X[len-2*stride-1];
   for (i=len-2*stride-1;i>=0;i--)
   {
      celt_norm x1, x2;

      x1 = Xptr[0];
      x2 = Xptr[stride];
      Xptr[stride] = EXTRACT16(MULT16_16_Q15_ADD(c,x2,s,x1));
      *Xptr--      = EXTRACT16(MULT16_16_Q15_SUB(c,x1,s,x2));
   }
}


/** Decode pulse vector and combine the result with the pitch vector to produce
    the final normalised signal in the current band. */
#define OVERRIDE_alg_unquant
unsigned alg_unquant(celt_norm *X, int N, int K, int spread, int B,
      ec_dec *dec, opus_val16 gain)
{
   int i;
   opus_val32 Ryy;
   int X0;
   unsigned collapse_mask;
   VARDECL(int, iy);
   SAVE_STACK;

   celt_assert2(K>0, "alg_unquant() needs at least one pulse");
   celt_assert2(N>1, "alg_unquant() needs at least two dimensions");
   ALLOC(iy, N, int);
   decode_pulses(iy, N, K, dec);
   Ryy = 0;
   i=0;
   asm volatile("mult $ac1, $0, $0");
   do {
        X0 = (int)iy[i];
        asm volatile("MADD $ac1, %0, %1" : : "r" (X0), "r" (X0));
   } while (++i < N);
   asm volatile("MFLO %0, $ac1" : "=r" (Ryy));
   normalise_residual(iy, X, N, Ryy, gain);
   exp_rotation(X, N, -1, B, K, spread);
   collapse_mask = extract_collapse_mask(iy, N, B);
   RESTORE_STACK;
   return collapse_mask;
}

#define OVERRIDE_renormalise_vector
void renormalise_vector(celt_norm *X, int N, opus_val16 gain)
{
   int i;
#ifdef FIXED_POINT
   int k;
#endif
   opus_val32 E = EPSILON;
   opus_val16 g;
   opus_val32 t;
   celt_norm *xptr = X;

   int X0, X2, X3, X1;
   asm volatile("mult $ac1, $0, $0");
   asm volatile("MTLO %0, $ac1" : :"r" (E));
   /*if(N %4)
       printf("error");*/
   for (i=0;i<N;i+=2)
   {
      X0 = (int)*xptr++;
      asm volatile("MADD $ac1, %0, %1" : : "r" (X0), "r" (X0));

        X1 = (int)*xptr++;
      asm volatile("MADD $ac1, %0, %1" : : "r" (X1), "r" (X1));
      
   }
     asm volatile("MFLO %0, $ac1" : "=r" (E));
#ifdef FIXED_POINT
   k = celt_ilog2(E)>>1;
#endif
   t = VSHR32(E, 2*(k-7));
   g = MULT16_16_P15(celt_rsqrt_norm(t),gain);

   xptr = X;
   for (i=0;i<N;i++)
   {
      *xptr = EXTRACT16(PSHR32(MULT16_16(g, *xptr), k+1));
      xptr++;
   }
   /*return celt_sqrt(E);*/
}

#endif /* __VQ_MIPSR1_H__ */
