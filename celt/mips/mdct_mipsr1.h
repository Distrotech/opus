/* Copyright (c) 2007-2008 CSIRO
   Copyright (c) 2007-2008 Xiph.Org Foundation
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

/* This is a simple MDCT implementation that uses a N/4 complex FFT
   to do most of the work. It should be relatively straightforward to
   plug in pretty much and FFT here.

   This replaces the Vorbis FFT (and uses the exact same API), which
   was a bit too messy and that was ending up duplicating code
   (might as well use the same FFT everywhere).

   The algorithm is similar to (and inspired from) Fabrice Bellard's
   MDCT implementation in FFMPEG, but has differences in signs, ordering
   and scaling in many places.
*/
#ifndef __MDCT_MIPSR1_H__
#define __MDCT_MIPSR1_H__

#ifndef SKIP_CONFIG_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#endif

#include "mdct.h"
#include "kiss_fft.h"
#include "_kiss_fft_guts.h"
#include <math.h>
#include "os_support.h"
#include "mathops.h"
#include "stack_alloc.h"

/* Forward MDCT trashes the input array */
#define OVERRIDE_clt_mdct_forward
void clt_mdct_forward(const mdct_lookup *l, kiss_fft_scalar *in, kiss_fft_scalar * OPUS_RESTRICT out,
      const opus_val16 *window, int overlap, int shift, int stride)
{
   int i;
   int N, N2, N4;
   kiss_twiddle_scalar sine;
   VARDECL(kiss_fft_scalar, f);
   VARDECL(kiss_fft_scalar, f2);
   SAVE_STACK;
   N = l->n;
   N >>= shift;
   N2 = N>>1;
   N4 = N>>2;
   ALLOC(f, N2, kiss_fft_scalar);
   ALLOC(f2, N2, kiss_fft_scalar);
   /* sin(x) ~= x here */
#ifdef FIXED_POINT
   sine = TRIG_UPSCALE*(QCONST16(0.7853981f, 15)+N2)/N;
#else
   sine = (kiss_twiddle_scalar)2*PI*(.125f)/N;
#endif

   /* Consider the input to be composed of four blocks: [a, b, c, d] */
   /* Window, shuffle, fold */
   {
      /* Temp pointers to make it really clear to the compiler what we're doing */
      const kiss_fft_scalar * OPUS_RESTRICT xp1 = in+(overlap>>1);
      const kiss_fft_scalar * OPUS_RESTRICT xp2 = in+N2-1+(overlap>>1);
      kiss_fft_scalar * OPUS_RESTRICT yp = f;
      const opus_val16 * OPUS_RESTRICT wp1 = window+(overlap>>1);
      const opus_val16 * OPUS_RESTRICT wp2 = window+(overlap>>1)-1;
      for(i=0;i<((overlap+3)>>2);i++)
      {
          *yp++ = S_MUL_ADD(*wp2, xp1[N2],*wp1,*xp2);
          *yp++ = S_MUL_SUB(*wp1, *xp1,*wp2, xp2[-N2]);

         xp1+=2;
         xp2-=2;
         wp1+=2;
         wp2-=2;
      }
      wp1 = window;
      wp2 = window+overlap-1;
      for(;i<N4-((overlap+3)>>2);i++)
      {
         /* Real part arranged as a-bR, Imag part arranged as -c-dR */
         *yp++ = *xp2;
         *yp++ = *xp1;
         xp1+=2;
         xp2-=2;
      }
      for(;i<N4;i++)
      {
          *yp++ =  S_MUL_SUB(*wp2, *xp2, *wp1, xp1[-N2]);
          *yp++ = S_MUL_ADD(*wp2, *xp1, *wp1, xp2[N2]);

         xp1+=2;
         xp2-=2;
         wp1+=2;
         wp2-=2;
      }
   }
   /* Pre-rotation */
   {
      kiss_fft_scalar * OPUS_RESTRICT yp = f;
      const kiss_twiddle_scalar *t = &l->trig[0];
      for(i=0;i<N4;i++)
      {
         kiss_fft_scalar re, im, yr, yi;
         re = yp[0];
         im = yp[1];
         yr = -S_MUL_ADD(re,t[i<<shift],im,t[(N4-i)<<shift]);
         yi = S_MUL_SUB(re,t[(N4-i)<<shift],im,t[i<<shift]);
         /* works because the cos is nearly one */
         *yp++ = yr + S_MUL(yi,sine);
         *yp++ = yi - S_MUL(yr,sine);
      }
   }

   /* N/4 complex FFT, down-scales by 4/N */
   opus_fft(l->kfft[shift], (kiss_fft_cpx *)f, (kiss_fft_cpx *)f2);

   /* Post-rotate */
   {
      /* Temp pointers to make it really clear to the compiler what we're doing */
      const kiss_fft_scalar * OPUS_RESTRICT fp = f2;
      kiss_fft_scalar * OPUS_RESTRICT yp1 = out;
      kiss_fft_scalar * OPUS_RESTRICT yp2 = out+stride*(N2-1);
      const kiss_twiddle_scalar *t = &l->trig[0];
      /* Temp pointers to make it really clear to the compiler what we're doing */
      for(i=0;i<N4;i++)
      {
         kiss_fft_scalar yr, yi;
         yr = S_MUL_ADD(fp[1],t[(N4-i)<<shift],fp[0],t[i<<shift]);
         yi = S_MUL_SUB(fp[0],t[(N4-i)<<shift],fp[1],t[i<<shift]);
         /* works because the cos is nearly one */
         *yp1 = yr - S_MUL(yi,sine);
         *yp2 = yi + S_MUL(yr,sine);;
         fp += 2;
         yp1 += 2*stride;
         yp2 -= 2*stride;
      }
   }
   RESTORE_STACK;
}

#define OVERRIDE_clt_mdct_backward
void clt_mdct_backward(const mdct_lookup *l, kiss_fft_scalar *in, kiss_fft_scalar * OPUS_RESTRICT out,
      const opus_val16 * OPUS_RESTRICT window, int overlap, int shift, int stride)
{
   int i;
   int N, N2, N4;
   kiss_twiddle_scalar sine;
   VARDECL(kiss_fft_scalar, f2);
   SAVE_STACK;
   N = l->n;
   N >>= shift;
   N2 = N>>1;
   N4 = N>>2;
   ALLOC(f2, N2, kiss_fft_scalar);
   /* sin(x) ~= x here */
#ifdef FIXED_POINT
   sine = TRIG_UPSCALE*(QCONST16(0.7853981f, 15)+N2)/N;
#else
   sine = (kiss_twiddle_scalar)2*PI*(.125f)/N;
#endif

   /* Pre-rotate */
   {
      /* Temp pointers to make it really clear to the compiler what we're doing */
      const kiss_fft_scalar * OPUS_RESTRICT xp1 = in;
      const kiss_fft_scalar * OPUS_RESTRICT xp2 = in+stride*(N2-1);
      kiss_fft_scalar * OPUS_RESTRICT yp = f2;
      const kiss_twiddle_scalar *t = &l->trig[0];
      for(i=0;i<N4;i++)
      {
         kiss_fft_scalar yr, yi;
         yr = S_MUL_SUB(*xp1,t[(N4-i)<<shift],*xp2, t[i<<shift]);
         yi =  -S_MUL_ADD(*xp2, t[(N4-i)<<shift],*xp1,t[i<<shift]);
         /* works because the cos is nearly one */
         *yp++ = yr - S_MUL(yi,sine);
         *yp++ = yi + S_MUL(yr,sine);
         xp1+=2*stride;
         xp2-=2*stride;
      }
   }

   /* Inverse N/4 complex FFT. This one should *not* downscale even in fixed-point */
   opus_ifft(l->kfft[shift], (kiss_fft_cpx *)f2, (kiss_fft_cpx *)(out+(overlap>>1)));

   /* Post-rotate and de-shuffle from both ends of the buffer at once to make
      it in-place. */
   {
      kiss_fft_scalar * OPUS_RESTRICT yp0 = out+(overlap>>1);
      kiss_fft_scalar * OPUS_RESTRICT yp1 = out+(overlap>>1)+N2-2;
      const kiss_twiddle_scalar *t = &l->trig[0];
      /* Loop to (N4+1)>>1 to handle odd N4. When N4 is odd, the
         middle pair will be computed twice. */
      for(i=0;i<(N4+1)>>1;i++)
      {
         kiss_fft_scalar re, im, yr, yi;
         kiss_twiddle_scalar t0, t1;
         re = yp0[0];
         im = yp0[1];
         t0 = t[i<<shift];
         t1 = t[(N4-i)<<shift];
         yr = S_MUL_SUB(re,t0,im,t1);
         yi = S_MUL_ADD(im,t0,re,t1);
         re = yp1[0];
         im = yp1[1];
         /* works because the cos is nearly one */
         yp0[0] = -(yr - S_MUL(yi,sine));
         yp1[1] = yi + S_MUL(yr,sine);

         t0 = t[(N4-i-1)<<shift];
         t1 = t[(i+1)<<shift];
         yr = S_MUL_SUB(re,t0,im,t1);
         yi = S_MUL_ADD(im,t0,re,t1);
         /* works because the cos is nearly one */
         yp1[0] = -(yr - S_MUL(yi,sine));
         yp0[1] = yi + S_MUL(yr,sine);
         yp0 += 2;
         yp1 -= 2;
      }
   }

   /* Mirror on both sides for TDAC */
   {
      kiss_fft_scalar * OPUS_RESTRICT xp1 = out+overlap-1;
      kiss_fft_scalar * OPUS_RESTRICT yp1 = out;
      const opus_val16 * OPUS_RESTRICT wp1 = window;
      const opus_val16 * OPUS_RESTRICT wp2 = window+overlap-1;

      for(i = 0; i < overlap/2; i++)
      {
         kiss_fft_scalar x1, x2;
         x1 = *xp1;
         x2 = *yp1;
         *yp1++ = MULT16_32_Q15(*wp2, x2) - MULT16_32_Q15(*wp1, x1);
         *xp1-- = MULT16_32_Q15(*wp1, x2) + MULT16_32_Q15(*wp2, x1);
         wp1++;
         wp2--;
      }
   }
   RESTORE_STACK;
}


#endif /* __MDCT_MIPSR1_H__ */
