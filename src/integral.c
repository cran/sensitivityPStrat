/* This is C-code for calculating gamma2 from Stute, 1995.  */

#include <R.h>
#include <Rmath.h>

void integral(int *lengthH, double *junk, double *omega, double *dH0, 
              double *dH1, double *Hadj, double *gamma0, double *phi, 
              double *gamma2)
{

  int i, j, k;
  double *nu = junk;
  double stuff;

  /* pre run redundant math ops  and Initialize gamma2 *
     phi[i] * gamma0[i] * dH1[i] -> phi[i]
     dH0[i] / ((1-Hadj[i]) * (1-Hadj[i])) -> dH0[i]
   */
  for (i=0; i<*lengthH; i++){
    gamma2[i] = 0;
    /* if phi[i] is zero then the product of gamma0[i] and dH1[i]
       will always be zero thus leaving phi[i] unchanged */
    if (phi[i] != 0)
    /* equvalent to 
       phi[i] = phi[i] * gamma0[i] * dH1[i] */
       phi[i] *= gamma0[i] * dH1[i];
    
    /* leave dH0[i] as zero if zero because 0/anything = 0 */
    if (dH0[i] != 0)
    /* equivalent to 
       dH0[i] = dH0[i]/((1 - Hadj[i]) * (1 - Hadj[i])) */
       dH0[i] /= (1 - Hadj[i]) * (1 - Hadj[i]);
  }

  for (i=0; i<*lengthH; i++) {
    for (j=0; j<*lengthH; j++) {
      /* if phi[j] is zero then stuff * phi[j] will be zero
         so this can be skipped because adding zero is the same
         as doing nothing */
      if (phi[j] == 0) {
        continue;
      }

      stuff = 0;
      for (k=0; k<*lengthH; k++) {
         /* if nu[k] < junk[i] and nu[k] < omega[j] then this 
            can be skipped because indicator would be zero and
            would just be adding zero to stuff */
        if (nu[k]<junk[i] && nu[k]<omega[j]) {
          /* This is the same as
             stuff += dH0[k]/((1-Hadj[k])*(1-Hadj[k])) */
          stuff += dH0[k];
        }
      }
      /* skip to next j if stuff = 0 because adding 0 is the 
         same as doing nothing */
      if(stuff == 0) {
        continue;
      }

      /* equivalent to
         gamma2[i] = stuff * phi[j] * gamma0[j] * dH1[j] */
      gamma2[i] += stuff*phi[j];
    }
  }
}
