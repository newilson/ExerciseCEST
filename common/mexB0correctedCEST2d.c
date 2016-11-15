/*
This code implements the 3D B0 corrected CEST map, positive and negative images at ppm of interest.
function [mapout,pimgout,nimgout] = getB0correctedCEST3d(ppm,ppmlist,posimages,negimages,B0map,mask,b0ppmstep)
*/

/* You can include any C libraries that you normally use */

#include "mex.h"   /*This C library is required*/
#include <math.h>
#include "functions/nrutil.h"

#define PI 3.141592653589793

/* Underlying Math routines */

extern void lfit(float *x,float *y, float *sig, int ndata,float *a, int ma, int *listan,int mfit,
                 float **covar,float *chisq, void (*)(float,float *,int)); 

void xpoly(float x,float *afunc,int mma)
{
	int i;

	afunc[1]=1.0;
	for (i = 2; i <= mma;i++) 
        afunc[i]=x*afunc[i-1];
}


void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;

    u = mxCalloc( n, sizeof(double) );
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else 
    {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) 
    {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
        qn=un=0.0;
	else
    {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	mxFree(u);
}
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
	int klo,khi,k;
	float h,b,a;

	klo=0;
	khi=n-1;
	while (khi-klo > 1) 
    {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Inside mexFunction---*/
    /*Declarations*/
    
    double  ppmval, *ppmlist, *ppmlistasc, *pimg, *nimg, *B0map, *mask, *pimgout, *nimgout, *mapout; /* Pointers to ppmlist,W,H,pos and negative images, and mask arrays */
    double maxppm, minppm, stepppm, B0val,temp;
    double zsp2p, zsp2n, *yp, *yn, ppmpos, ppmneg;
    mwSize nppm, nintppm, ndim, nelem, ielem, ippm1, ippm2, ndimimg, nelemimg;
    const mwSize *dim_array_mask, *dim_array_image;
    mwSize *psort;
    
    mwSize NDEG, i,j,k, intfactor=10;
    double *xint,*ypint,*ynint, *y2;
    
    double *xsc, mean, std;
    
    int i1, i2;
    int *listap, *listan, NTERM, NPTS;
    float *xzsp, *afunc;
	float *yp1, chisqp,*ap,*sigmap,**covarp;
	float *yn1, chisqn,*an,*sigman,**covarn;
    
    void (*fp) (float, float *, int);
    fp = xpoly;

    NDEG = 3;
    NTERM = NDEG+1;
    
    ppmval = mxGetScalar(prhs[0]);
    mexPrintf("ppmval = %f\n", ppmval);
    
    ppmlist = mxGetPr(prhs[1]);
    nppm=mxGetNumberOfElements(prhs[1]);
    nintppm = 1 + (nppm-1)*intfactor;
    NPTS = nintppm;
    
    for (ippm1 = 0; ippm1 < nppm; ippm1++)
        mexPrintf(" ppmlist[%i] = %f \n", ippm1, ppmlist[ippm1]);
    
    pimg = mxGetPr(prhs[2]);
    nimg = mxGetPr(prhs[3]);
    nelemimg = mxGetNumberOfElements(prhs[2]);
    ndimimg=mxGetNumberOfDimensions(prhs[2]);
    
    B0map = mxGetPr(prhs[4]);
    mask = mxGetPr(prhs[5]);
    
    /* Get the dimensions in the mask array */
    nelem = mxGetNumberOfElements(prhs[5]);
    ndim=mxGetNumberOfDimensions(prhs[5]);
    dim_array_mask=mxGetDimensions(prhs[5]);
            
    plhs[0] = mxCreateNumericArray(ndim,dim_array_mask,mxDOUBLE_CLASS,mxREAL);
    mapout  = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(ndim,dim_array_mask,mxDOUBLE_CLASS,mxREAL);
    pimgout = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(ndim,dim_array_mask,mxDOUBLE_CLASS,mxREAL);
    nimgout = mxGetPr(plhs[2]);
    
   
    yp = (double *)mxCalloc(nppm,sizeof(double)); /* positive yvalues for fit */
    yn = (double *)mxCalloc(nppm,sizeof(double)); /* negative yvalues for fit */
    y2 = (double *)mxCalloc( nppm, sizeof(double) );  /* storage for second derivative */

    xsc   = (double *)mxCalloc(nintppm,sizeof(double)); /* Scaled to mean and std of xin */
    xint  = (double *)mxCalloc(nintppm,sizeof(double)); /* Scaled to mean and std of xin */
    ypint = (double *)mxCalloc(nintppm,sizeof(double)); /* Scaled to mean and std of xin */
    ynint = (double *)mxCalloc(nintppm,sizeof(double)); /* Scaled to mean and std of xin */
                
    ppmlistasc = mxCalloc(nppm, sizeof(double));
    psort      = mxCalloc(nppm, sizeof(int));
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
    for (ippm1 =0; ippm1 < nppm; ippm1++)
    {
        ppmlistasc[ippm1] = minppm + ippm1*stepppm;
    }
    for (ippm1 = 0; ippm1 < nppm; ippm1++)
    {
        psort[ippm1] = ippm1;
        for (ippm2 = 0; ippm2 < nppm; ippm2++)
        {
            if ( fabs(ppmlistasc[ippm2] - ppmlist[ippm1]) < 0.001 )
                psort[ippm1] = ippm2;
        }
    }
        
    mean = 0.0;
    std = 0.0;    
    for (ippm1 = 0; ippm1 < nintppm; ippm1++)
    {
        xint[ippm1] = minppm + ippm1 * (stepppm/intfactor);
        mean += xint[ippm1];
    }
    mean /= nintppm;
    for (ippm1 = 0; ippm1 < nintppm; ippm1++)
    {
        std += pow((xint[ippm1]-mean),2);
    }
    std = sqrt( std/(nintppm-1) );
    mexPrintf("mean = %f, std = %f\n\n", mean,std);
    
    for (ippm1 = 0; ippm1 < nintppm; ippm1++)
    {
        xsc[ippm1] = (xint[ippm1]-mean)/std;
    }
    
    afunc  = vector(1,NTERM);	
	xzsp   = vector( 1, NPTS );
    for (i1 =0; i1 < NPTS; i1++)
    {
        xzsp[i1+1] = xint[i1];       
    }
    
	listap = ivector( 1,NTERM);
	ap     = vector ( 1, NTERM );
	covarp = matrix( 1, NTERM, 1, NTERM );
    
	listan = ivector( 1,NTERM);
    for (i2 = 1; i2 <= NTERM; i2++)
    {
        listap[i2] = i2;
        listan[i2] = i2;
    }
    
	an     = vector ( 1, NTERM );
	covarn = matrix( 1, NTERM, 1, NTERM );
    
	yp1     = vector( 1 , NPTS );
	sigmap = vector( 1, NPTS );
	yn1     = vector( 1 , NPTS );
	sigman = vector( 1, NPTS );
      
    for (ielem = 0; ielem < nelem; ielem++) 
    {
        pimgout[ielem] = 0.0;
        nimgout[ielem] = 0.0;
        mapout[ielem]  = 0.0;
        if (mask[ielem] > 0) 
        {
            B0val = B0map[ielem];
            ppmpos = ppmval + B0val;
            if (ppmpos > maxppm)
                ppmpos = maxppm;
            if (ppmpos < minppm)
                ppmpos = minppm;
            ppmneg = ppmval - B0val;
            if (ppmneg > maxppm)
                ppmneg = maxppm;
            if (ppmneg < minppm)
                ppmneg = minppm;

            for  (ippm1 = 0; ippm1 < nppm; ippm1++) 
            {
                i1 = psort[ippm1]*nelem + ielem;
                yp[ippm1] = pimg[i1];
                yn[ippm1] = nimg[i1];
            } 

            /* Call spline to get second derivatives and splint to interpolate Positive side */
            spline(ppmlistasc, yp, nppm, yp[0], yp[nppm-1], y2);
            for (i2 =0; i2< nintppm; i2++) 
            {
                splint(ppmlistasc, yp, y2, nppm, xint[i2], &ypint[i2]);
            }
            
            /* Call spline to get second derivatives and splint to interpolate Negative side */
            spline(ppmlistasc, yn, nppm, yn[0], yn[nppm-1], y2);
            for (i2 =0; i2< nintppm; i2++) 
            {
                splint(ppmlistasc, yn, y2, nppm, xint[i2], &ynint[i2]);
            }
            
            /* Now call polyfit routine for positive side to get polynomial coeffs */
            for (i2 = 0; i2 < NPTS; i2++)
            {
                yp1[i2+1] = ypint[i2];
                yn1[i2+1] = ynint[i2];
                sigmap[i2+1] = 1.0;
                sigman[i2+1] = 1.0;
            }
            lfit(xzsp,yp1,sigmap,NPTS,ap,NTERM,listap,NTERM,covarp,&chisqp,fp);
            xpoly(ppmpos,afunc,NTERM);
            zsp2p = 0.0;
            for (i2=1;i2<=NTERM;i2++)
                zsp2p += (afunc[i2]*ap[i2]);
            
            lfit(xzsp,yn1,sigman,NPTS,an,NTERM,listap,NTERM,covarp,&chisqp,fp);
            xpoly(ppmneg,afunc,NTERM);
            zsp2n = 0.0;
            for (i2=1;i2<=NTERM;i2++)
                zsp2n += (afunc[i2]*an[i2]);
            
            /* Load B0 corrected  CEST images and map */
            pimgout[ielem] = zsp2p;
            nimgout[ielem] = zsp2n;
            mapout[ielem] = 100.0 * (zsp2n -zsp2p)/(0.001+zsp2n);
        }
    }
    
    mxFree(psort);
    mxFree(xsc);
    mxFree(xint);
    mxFree(ypint);
    mxFree(ynint);
    mxFree(yp);
    mxFree(yn);
    
	free_vector(xzsp,1,NPTS);
    
    free_matrix(covarp,1,NTERM,1,NTERM);
	free_vector(sigmap,1,NPTS);
	free_vector(yp1,1,NPTS);
	free_vector(ap,1,NTERM);
	free_ivector(listap,1,NTERM);
    
    free_matrix(covarn,1,NTERM,1,NTERM);
	free_vector(sigman,1,NPTS);
	free_vector(yn1,1,NPTS);
	free_vector(an,1,NTERM);
	free_ivector(listan,1,NTERM);
    
    
    return;
}
