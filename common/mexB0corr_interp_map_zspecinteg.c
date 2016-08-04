/*You can include any C libraries that you normally use*/
#include "mex.h"   /*This C library is required*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * This code implements the B0 corrected CEST map calculation through interpolation of full z spectrum.
 * function [pimg,nimg] = mexB0correctedCESTimginterp(ppm,ppmlist,posimages,negimages,B0map,mask,b0ppmstep)
 * ppmlist should contain 0 ppm and be sorted.
*/

#define PI 3.141592653589793

void spline (double	x[], double	y[], int	n, double	yp1, double	ypn, double	y2[])
{
    
    int	i,k;
    double	p,qn,sig,un,*u;
    
    u   = (double *) mxCalloc( n, sizeof(double) );
    if(yp1 > 0.99e30)
        y2[0] = u[0] = 0.0;
    else
    {
        y2[0] = -0.5;
        u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for(i = 1; i < n-1; i++)
    {
        sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
        p = sig*y2[i-1] + 2.0;
        y2[i] = (sig - 1.0)/p;
        u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
        u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
    }
    if(ypn > 0.99e30)
        qn = un = 0.0;
    else
    {
        qn = 0.5;
        un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
    }
    y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);
    for(k = n-2; k >= 0; k--)
    {
        y2[k] = y2[k]*y2[k+1] + u[k];
    }
    
    mxFree(u);
}

void splint( double	xa[], double	ya[], double	y2a[],int	n, double	x, double	*y)
{
    
    int		klo,khi,k;
    double		h,b,a;
    static int	pklo=0,pkhi=1;
    
    /*
  Based on the assumption that sequential calls to this function are made
  with closely-spaced, steadily-increasing values of x, I first try using
  the same values of klo and khi as were used in the previous invocation.
  If that interval is no longer correct, I do a binary search for the
  correct interval.
     */
    if(xa[pklo] <= x && xa[pkhi] > x)
    {
        klo = pklo;
        khi = pkhi;
    }
    else
    {
        klo = 0;
        khi = n - 1;
        while(khi - klo > 1)
        {
            k = (khi + klo) >> 1;
            if(xa[k] > x) khi = k;
            else          klo = k;
        }
    }
    
    h = xa[khi] - xa[klo];
    if(h == 0)
    {
        mexPrintf("Bad xa input to function splint() setting y value to zero\n");
        *y = 0.0;
    }
    else
    {
        a = (xa[khi] - x)/h;
        b = (x - xa[klo])/h;
        *y = a*ya[klo] + b*ya[khi] +
                ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Inside mexFunction---*/
 /* function fmap = getB0correctedCESTmap(ppm,ppmlist,posimages,negimages,B0map,mask,b0ppmstep) */
    /*Declarations*/
    
    double  ppmval, *ppmlist, *pimg, *nimg, *B0map, *mask, *pimgout, *nimgout, *zimgout,*zppmout; /* Pointers to ppmlist,W,H,pos and negative images, and mask arrays */
    double maxppm, minppm, stepppm, B0val, maxB0, minB0,interp_ppmstep;
    double ppmpos, ppmneg,zsp2p,zsp2n;
    mwSize nppm, ndim, nelem, ielem, ippm1, ippm2, ndimimg, nelemimg,startindex,endindex,nfit;
    const mwSize *dim_array_mask, *dim_array_image;
    mwSize H,W,dim_array_zimg[3];
    double f,x,y,yp1,ypn,*xa,*y2,temp1,temp2;
    double *rawz, *interpz, *interpzs, *xinterp;

    mwSize i1, i2, i3, i4, nzppm, ninterp, nppmout,ninterp2;
    
    ppmval = mxGetScalar(prhs[0]);
    
    ppmlist = mxGetPr(prhs[1]);
    nppm=mxGetNumberOfElements(prhs[1]);
    nzppm = 2*nppm-1;
    
    pimg = mxGetPr(prhs[2]);
    nimg = mxGetPr(prhs[3]);
    nelemimg = mxGetNumberOfElements(prhs[2]);
    ndimimg=mxGetNumberOfDimensions(prhs[2]);
    dim_array_image=mxGetDimensions(prhs[2]);
    
    B0map = mxGetPr(prhs[4]);
    nelem = mxGetNumberOfElements(prhs[4]);
    maxB0 = -1000.0;
    minB0 = 1000.0;
    for (ielem = 0; ielem<nelem; ielem++)
    {
        if (B0map[ielem] > maxB0)
            maxB0 = B0map[ielem];
        if (B0map[ielem] < minB0)
            minB0 = B0map[ielem];
    }
    mexPrintf(" minB0 = %f maxB0 = %f  \n", minB0, maxB0);
    if (maxB0 < -minB0)
        maxB0 = -minB0;
    else
        minB0 = -maxB0;
    
    mask = mxGetPr(prhs[5]);
    /* Get the dimensions in the mask array */
    nelem = mxGetNumberOfElements(prhs[5]);
    ndim=mxGetNumberOfDimensions(prhs[5]);
    dim_array_mask=mxGetDimensions(prhs[5]);
    H = dim_array_mask[0];
    W = dim_array_mask[1];

    plhs[0]=mxCreateDoubleMatrix(H,W,mxREAL);
    pimgout=(double *)mxGetData(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(H,W,mxREAL);
    nimgout=(double *)mxGetData(plhs[1]);
        
    
    maxppm = -10000.0;
    minppm = 10000.0;
    for (ippm1 =0; ippm1 < nppm; ippm1++)
    {
        if (ppmlist[ippm1] > maxppm)
            maxppm = ppmlist[ippm1];
        if (ppmlist[ippm1] < minppm)
            minppm = ppmlist[ippm1];       
    }
    stepppm = (maxppm - minppm)/(nppm-1);
    
    if (maxB0 > (maxppm - 3*stepppm) )
    {
        maxB0 = (maxppm - 3*stepppm);
        minB0 = -maxB0;
    }
    mexPrintf(" minB0 = %f maxB0 = %f  \n", minB0, maxB0);
                
    xa   = (double *) mxCalloc( nzppm, sizeof(double) );
    for (ippm1 = 0; ippm1 < nzppm; ippm1++)
        xa[ippm1] = -maxppm  + ippm1*stepppm;
    
    rawz = (double *) mxCalloc( nzppm, sizeof(double) );
    y2   = (double *) mxCalloc( nzppm, sizeof(double) );
    
    interp_ppmstep = stepppm/100;
    if (minB0 < 0.0)
        temp1 = xa[0] + stepppm * (int) ( 1.0 - (minB0/stepppm) );
    else
        temp1 = xa[0] + stepppm * (int) ( 1.0 + (minB0/stepppm) );
    
    if (maxB0 > 0.0)
        temp2 = xa[nzppm-1] - stepppm * (int)( 1.0 + (maxB0/stepppm) );
    else
        temp2 = xa[nzppm-1] - stepppm * (int)( 1.0 - (maxB0/stepppm) );
    mexPrintf(" minB0 = %f maxB0 = %f minppmin = %f maxppmin %f stepppm = %f interpstep = %f \n", minB0, maxB0, minppm, maxppm, stepppm, interp_ppmstep);
    
    nppmout = 1 + (temp2-temp1+0.0001)/stepppm;
    dim_array_zimg[0] = H;
    dim_array_zimg[1] = W;
    dim_array_zimg[2] = nppmout;
    plhs[2] = mxCreateNumericArray(3, dim_array_zimg, mxDOUBLE_CLASS, mxREAL);
    zimgout = (double *)mxGetData(plhs[2]);
    
    plhs[3]=mxCreateDoubleMatrix(nppmout,1,mxREAL);
    zppmout=(double *)mxGetData(plhs[3]);
    for (ippm1 = 0; ippm1<nppmout;ippm1++)
        zppmout[ippm1] = temp1 + ippm1*stepppm;
    mexPrintf(" minppmout = %f maxppmout = %f inppmout - %i\n", temp1, temp2,nppmout);
    
          
    ninterp = 1 + (nzppm-1)*100;
    ninterp2 = ninterp/2;
    xinterp = (double *)mxCalloc( ninterp, sizeof(double) );
    interpz = (double *)mxCalloc( ninterp, sizeof(double) );
    interpzs = (double *)mxCalloc( ninterp, sizeof(double) );
    
    for ( ippm1 = 0; ippm1 < ninterp; ippm1++)
        xinterp[ippm1] = -maxppm + ippm1*interp_ppmstep;
    
    for (ielem = 0; ielem < nelem; ielem++) 
    {
        if (mask[ielem] > 0) 
        {            
            for (i2 = 1; i2 < nppm; i2++)
            {
                i3 = (nppm-i2)*nelem + ielem;
                rawz[i2-1] = nimg[i3];
                i1 = i2 * nelem + ielem;
                rawz[i2+nppm-1] = pimg[i1];
            }
            rawz[nppm-1] = 0.5 * ( pimg[ielem] + nimg[ielem] );
            yp1 = rawz[0];
            ypn = rawz[nzppm-1];
            
            
            /* Call spline to get second derivatives */
            spline(xa,rawz,nzppm,yp1,ypn,y2);
            
            
            /* Call splint for interpolations */
            for (i2 = 0; i2 < ninterp; i2++)
            {
                splint(xa,rawz,y2,nzppm,xinterp[i2],&interpz[i2]);
            }            
            B0val = B0map[ielem];
            i1 = (B0val/interp_ppmstep);     /* Index to shift z spectrum based on B0val */
            
            if ( i1 >= 0)
            {
                for (i2 = 0; i2 < ninterp; i2++)
                {
                    if ((i2+i1) < ninterp)
                        interpzs[i2] = interpz[i2+i1];
                    else
                        interpzs[i2] = 0.0;
                }
            }
            else
            {
                for (i2 = ninterp-1; i2 >= 0; i2--)
                {
                    if ((i2+i1) >= 0)
                        interpzs[i2] = interpz[i2+i1];
                    else
                        interpzs[i2] = 0.0;
                }
            }
            
            ppmpos = ppmval;
            ppmneg = -ppmval;
            nfit = (0.2/interp_ppmstep);  /* find mean around ppmval +/-0.1 ppm for final B0 corrected spectral values */
            zsp2p = 0.0;
            zsp2n = 0.0;
            for (i4 = 0; i4 <= nfit; i4++)
            {
                i2 = (ppmpos/interp_ppmstep) + ninterp2 + (i4-(nfit/2));
                zsp2p += interpzs[i2];
                i2 = (ppmneg/interp_ppmstep) + ninterp2 + (i4-(nfit/2));
                zsp2n += interpzs[i2];
            }
            
            /* Calculate B0 corrected  CEST map */
            pimgout[ielem] = zsp2p/nfit;
            nimgout[ielem] = zsp2n/nfit;
            
            for (i3 = 0; i3 < nppmout; i3++)
            {
                i1 = (i3*nelem)+ielem;
                i2 = ninterp2 + (zppmout[i3]/interp_ppmstep);
                zimgout[i1] = interpzs[i2];
            }
                
        }
        else
        {
            pimgout[ielem] = 0.0;
            nimgout[ielem] = 0.0;
            for (i3 = 0; i3 < nppmout; i3++)
            {
                i1 = (i3*nelem)+ielem;
                zimgout[i1] = 0.0;
            }
        }
    }
    
    mxFree(interpzs);
    mxFree(interpz);
    mxFree(xinterp);
    mxFree(rawz);
    mxFree(y2);
    mxFree(xa);
    
    return;
}
