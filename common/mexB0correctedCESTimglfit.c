/*You can include any C libraries that you normally use*/
#include "mex.h"   /*This C library is required*/
#include <math.h>
#include "functions/nrutil.h"

#define PI 3.141592653589793


extern void lfit(float *x,float *y, float *sig, int ndata,float *a, int ma, int *listan,int mfit,
                 float **covar,float *chisq, void (*)(float,float *,int)); 

void xpoly(float x,float *afunc,int mma)
{
	int i;

	afunc[1]=1.0;
	for (i = 2; i <= mma;i++) 
        afunc[i]=x*afunc[i-1];
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Inside mexFunction---*/
 /* function [pimg,nimg] = getB0correctedCESTmap(ppm,ppmlist,posimages,negimages,B0map,mask,b0ppmstep) */
    /*Declarations*/
    
    double  ppmval, *ppmlist, *pimg, *nimg, *B0map, *mask, *b0ppmstepptr, b0ppmstep, *pimgout,*nimgout; /* Pointers to ppmlist,W,H,pos and negative images, and mask arrays */
    double maxppm, minppm, stepppm, B0val;
    double zsp2p, zsp2n;
    mwSize nppm, ndim, nelem, ielem, ippm1, ippm2, ndimimg, nelemimg;
    const mwSize *dim_array_mask, *dim_array_image;
    mwSize *psort,H,W;

    void (*fp) (float, float *, int);
    
    int i1, i2;
    int *listap, *listan, NTERM, NPTS;
    float *xzsp, *afunc, ppmpos, ppmneg;
	float *yp, chisqp,*ap,*sigmap,**covarp;
	float *yn, chisqn,*an,*sigman,**covarn;
    NTERM = 3;
    fp = xpoly;

    ppmval = mxGetScalar(prhs[0]);
    
    ppmlist = mxGetPr(prhs[1]);
    nppm=mxGetNumberOfElements(prhs[1]);
    NPTS = nppm;
    
    pimg = mxGetPr(prhs[2]);
    nimg = mxGetPr(prhs[3]);
    nelemimg = mxGetNumberOfElements(prhs[2]);
    ndimimg=mxGetNumberOfDimensions(prhs[2]);
    dim_array_image=mxGetDimensions(prhs[2]);
    
    B0map = mxGetPr(prhs[4]);
    mask = mxGetPr(prhs[5]);
    /* Get the dimensions in the mask array */
    nelem = mxGetNumberOfElements(prhs[5]);
    ndim=mxGetNumberOfDimensions(prhs[5]);
    dim_array_mask=mxGetDimensions(prhs[5]);
    H = dim_array_mask[0];
    W = dim_array_mask[1];

    if (nrhs > 6)
        b0ppmstep = mxGetScalar(prhs[6]);
    else
        b0ppmstep = 0.005;
        
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
    
    
    afunc  = vector(1,NTERM);	
	xzsp   = vector( 1, NPTS );
    for (i1 =0; i1 < NPTS; i1++)
    {
        xzsp[i1+1] = minppm + i1*stepppm;       
    }
    
	listap = ivector( 1,NTERM);
	ap     = vector ( 1, NTERM );
	covarp = matrix( 1, NTERM, 1, NTERM );
    
	listan = ivector( 1,NTERM);
	an     = vector ( 1, NTERM );
	covarn = matrix( 1, NTERM, 1, NTERM );
    
	yp     = vector( 1 , NPTS );
	sigmap = vector( 1, NPTS );
	yn     = vector( 1 , NPTS );
	sigman = vector( 1, NPTS );
      
    psort      = mxCalloc(nppm, sizeof(int));
    for (ippm1 = 0; ippm1 < nppm; ippm1++)
    {
        psort[ippm1] = floor( 0.001 + (ppmlist[ippm1] - minppm)/stepppm );
    }
        
    for (ielem = 0; ielem < nelem; ielem++) 
    {
        pimgout[ielem] = 0.0;
        nimgout[ielem] = 0.0;
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
                yp[ippm1+1] = pimg[i1];
                sigmap[ippm1+1] = 1.0;
                yn[ippm1+1] = nimg[i1];
                sigman[ippm1+1] = 1.0;
            } 
            for (i2 = 1; i2 <= NTERM; i2++)
            {
                listap[i2] = i2;
                listan[i2] = i2;
            }
/* Use polyfit + polyeval equivalent on yp and yn to get zsp2p and zsp2n  
 	lfit(x,y,sig,NPT,a,NTERM,lista,mfit,covar,&chisq,xpoly);
*/
            lfit(xzsp,yp,sigmap,NPTS,ap,NTERM,listap,NTERM,covarp,&chisqp,fp);
            xpoly(ppmpos,afunc,NTERM);
            zsp2p = 0.0;
            for (i2=1;i2<=NTERM;i2++)
                zsp2p += (afunc[i2]*ap[i2]);
            pimgout[ielem] = zsp2p;

            lfit(xzsp,yn,sigman,NPTS,an,NTERM,listan,NTERM,covarn,&chisqn,fp);
            xpoly(ppmneg,afunc,NTERM);
            zsp2n = 0.0;
            for (i2=1;i2<=NTERM;i2++)
                zsp2n += (afunc[i2]*an[i2]);
            nimgout[ielem] = zsp2n;
        }
    }
    
    mxFree(psort);
    
	free_vector(xzsp,1,nppm);
    
    free_matrix(covarp,1,NTERM,1,NTERM);
	free_vector(sigmap,1,nppm);
	free_vector(yp,1,nppm);
	free_vector(ap,1,NTERM);
	free_ivector(listap,1,NTERM);
    
    free_matrix(covarn,1,NTERM,1,NTERM);
	free_vector(sigman,1,nppm);
	free_vector(yn,1,nppm);
	free_vector(an,1,NTERM);
	free_ivector(listan,1,NTERM);
    
    return;
}
