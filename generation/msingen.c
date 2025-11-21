// To compile into a mex file in Matlab, use
// mex -D__MATLAB msingen.c fft.c

/* Call Library source file */
#ifdef __LABVIEW
#include "extcode.h"
#endif

#ifdef __MATLAB
#include "mex.h"
#endif

#include "fft.h"
#include <math.h>

#define TWO_PI (6.28318530717959)
#define ITNO 200
#define EPSILON 2.220446049250313e-016
/* EPSILON is the distance from 1.0 to the next double precision floating point number.
   See Matlab help for "eps" function for more info. */

#ifdef __LABVIEW
_declspec(dllexport) long msingen(unsigned long numFFTpoints, unsigned long numFreq, 
	double multisinex[]);
#endif
void vectorCopy(double *vec1, double *vec2, unsigned long n);
double interpolate_clx(double crx);
void getTimeFunction(double *cx, unsigned long numFreq, double *multisinex, unsigned long numFFTpoints, 
					 double *crx);
double sign(double x);

#ifdef __LABVIEW
_declspec(dllexport) long msingen(unsigned long numFFTpoints, unsigned long numFreq, 
	double multisinex[])
#endif
#ifdef __MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
#endif
{
#ifdef __MATLAB
	unsigned long numFFTpoints = (unsigned long)mxGetScalar(prhs[0]);
	unsigned long numFreq = (unsigned long)mxGetScalar(prhs[1]);
#endif

	double *ampv;               /* Amplitude vector */
	double *cx;                 /* Vector of complex amplitudes;
						           Storage of complex numbers: cx is a vector of length 2*numFreq;
						           then cx[2*i] = Re(cx[i]) and cx[2*i+1] = Im(cx[i]); */
	double crx, crxold, crxopt; /* crest factors */
	double *cxopt;              /* optimum vector or complex amplitudes - see description of cx */
	int nodecrease;             /* counter of nonimproving steps before increase clipping level */
	int crxconst;               /* counter of steps  with no change in crx*/
	int iopt;                   /* number of the optimal interation cycle */
	double clx;                 /* clipping */
	const double clxyuw1=0.1, clxyuw2=0.1;   /* clx and cly update weightings */
	unsigned long numClip;      /* number of time function values to clip */
	double min5absValues[5];    /* smallest 5 absolute values in time function */

	unsigned long i,j; int k, kk; /* Loop counters */
	double phase; /* random phase angle while setting up initial vector of complex amplitudes */
	double magfft, refft, imfft;

// In LabVIEW, this array gets allocated in the LabVIEW part
#ifdef __MATLAB
	double *multisinex;
	plhs[0] = mxCreateDoubleMatrix(1, numFFTpoints, mxREAL);
	multisinex = mxGetPr(plhs[0]);
#endif

	/* DEBUG STUFF
	double dbgmsinex[2048], dbgcx[64];
	float dbgvar1, dbgvar2;
	int idbg;
	char dbgstr[80];
	FILE *dbgfile1; */

	/* Setup Frequency and desired amplitude vectors */
	/* ampv = (double *)malloc(numFreq*sizeof(double)); */
#ifdef __LABVIEW
	ampv = (double *)AZNewPtr(numFreq*sizeof(double));
#endif
#ifdef __MATLAB
	ampv = (double *)malloc(numFreq*sizeof(double));
#endif
	for (i=0; i<numFreq; i++) ampv[i] = 1.0;

	/* Create initial vector is complex amplitudes */
	/* cx = (double *)malloc(2*numFreq*sizeof(double));
         cxopt = (double *)malloc(2*numFreq*sizeof(double)); */
#ifdef __LABVIEW
	cx = (double *)AZNewPtr(2*numFreq*sizeof(double));
	cxopt = (double *)AZNewPtr(2*numFreq*sizeof(double));
#endif
#ifdef __MATLAB
	cx = (double *)malloc(2*numFreq*sizeof(double));
	cxopt = (double *)malloc(2*numFreq*sizeof(double));
#endif
	for (i=0; i<numFreq; i++) {
		/* phase = TWO_PI*(double)rand()/((double)RAND_MAX + 1.0); --> This is the call in windows, but in labview, */
#ifdef __LABVIEW        
		RandomGen(&phase); phase*=TWO_PI;
#endif
#ifdef __MATLAB
        phase = TWO_PI*(double)rand()/((double)RAND_MAX + 1.0);
#endif
		cx[2*i] = ampv[i]*cos(phase);
		cx[2*i+1] = ampv[i]*sin(phase);
	}

	/* DEBUG STUFF
	dbgfile1 = fopen("cxreal.dat","r");
	for (i=0; i<numFreq; i++) {
		fgets(dbgstr,80,dbgfile1);
		sscanf(dbgstr,"%f",&dbgvar1);
		cx[2*i] = dbgvar1;
		fgets(dbgstr,80,dbgfile1);
		sscanf(dbgstr,"%f",&dbgvar2);
		cx[2*i+1] = dbgvar2;
	}
	fclose(dbgfile1);
	for (idbg=0; idbg<64; idbg++) dbgcx[idbg] = cx[idbg]; */

	crx = (double)numFFTpoints;
	nodecrease = 0;
	crxconst = 0;
	iopt = 0;
	/* Optimization Iterations */
	for (i=0; i<=ITNO; i++) {
		
		if (i>=1) { /* No clipping for i = 0 */
			if (i==1) clx = interpolate_clx(crx); /* Set clipping level */
			
			if (crx < crxold) { /* Improvement found */
				clx = clx - clxyuw2*(1.0-clx);
				clx = (clx > 0.2) ? clx : 0.2;
			}
			else {
				clx=clx+clxyuw1*(1.0-clx);
			} /* Improvement found */

			/* Clipping */
			numClip = 0;
			for (j=0; j<numFFTpoints; j++) {
				if (fabs(multisinex[j])>clx*crx) 
					numClip++;
			}
			if (numClip > numFFTpoints - 5) {
				for (k=0; k<5; k++) min5absValues[k] = fabs(multisinex[k]);
				for (j=5; j<numFFTpoints; j++) {
					for (k=0; k<5; k++) {
						if (fabs(multisinex[j]) < min5absValues[k]) {
							for (kk=4; kk>k; kk--)
								min5absValues[kk] = min5absValues[kk-1];
							min5absValues[k] = fabs(multisinex[j]);
							break;
						}
					}
				}
				clx = min5absValues[4]/crx;
			}
			for (j=0; j<numFFTpoints; j++) {
				if (fabs(multisinex[j])>clx*crx) 
					multisinex[j] = clx*crx*sign(multisinex[j]);
			}
			/* End clipping */

			/* DEBUG STUFF
			for (idbg=0; idbg<2048; idbg++) dbgmsinex[idbg] = multisinex[idbg]; */

			/* Find New Complex Amplitudes */
			realft(multisinex, numFFTpoints, 1);
			for (j=0; j<numFreq; j++) {
				refft = multisinex[2*(j+1)];
				imfft = multisinex[2*(j+1)+1];
				if (refft*refft + imfft*imfft > 0) {
					magfft = sqrt(refft*refft + imfft*imfft);
					cx[2*j] = ampv[j]*refft/magfft;
					cx[2*j+1] = -ampv[j]*imfft/magfft;
				}
			}
			/* DEBUG STUFF
			for (idbg=0; idbg<64; idbg++) dbgcx[idbg] = cx[idbg]; */
		} /* if i>= 1*/

		crxold = crx;
		getTimeFunction(cx, numFreq, multisinex, numFFTpoints, &crx);
		/* DEBUG STUFF
		for (idbg=0; idbg<2048; idbg++) dbgmsinex[idbg] = multisinex[idbg]; */

		/* Check for optimum */
		if (i==0) {
			vectorCopy(cx, cxopt, 2*numFreq);
			crxopt = crx;
		}
		else {
			if (crx < crxopt - (1.01e-4)) {
				vectorCopy(cx, cxopt, 2*numFreq);
				crxopt = crx;
				nodecrease = 0;
				iopt = i;
			}
			else {
				nodecrease ++;
				if ( fabs(crx-crxold)<10*EPSILON )
					crxconst++;
				else
					crxconst=0;
			}
		}

		if ( ((i!=0) && (clx>=0.99)) || (crxconst>10) ) break;

	} /* End of Optimization Iterations */
	
	/* DEBUG STUFF 
	for (idbg=0; idbg<64; idbg++) dbgcx[idbg] = cxopt[idbg]; */
	getTimeFunction(cxopt, numFreq, multisinex, numFFTpoints, &crx);
	for (i=0; i<numFFTpoints; i++) multisinex[i]/=crx;

	/* Clean up */
#ifdef __LABVIEW
    AZDisposePtr(ampv);
	AZDisposePtr(cx);
	AZDisposePtr(cxopt);
#endif
#ifdef __MATLAB
    free(ampv);
    free(cx);
    free(cxopt);
#endif
	return 0;
}

/**************************************************************************************************************/

void vectorCopy(double *vec1, double *vec2, unsigned long n)
{
	unsigned long i;

	for (i=0; i<n; i++) vec2[i] = vec1[i];
}

/**************************************************************************************************************/

double interpolate_clx(double crx)
{
	double a[] = {0,1.6,2,3.5,4.5,1e6};
	double b[] = {.9,.9,.8,.3,.2,.2};
	int i=0;

	while (crx>a[i]) i++;

	if (i==0)
		return b[i];
	else
		return b[i-1] + (b[i] - b[i-1])/(a[i] - a[i-1])*(crx - a[i-1]);
}

/**************************************************************************************************************/

void getTimeFunction(double *cx, unsigned long numFreq, double *multisinex, unsigned long numFFTpoints, 
					 double *crx)
{
	double rms;                 /* RMS value of multisine function */
	double maxms, minms;        /* max and min values in the multisine */

	unsigned long j; /* loop counter */

	/* DEBUG STUFF
	double dbgmsinex[2048];
	int idbg; 
	FILE *dbgfile2; */

	for (j=0; j<numFFTpoints; j++) multisinex[j] = 0;
	for (j=0; j<numFreq; j++) {
		multisinex[2*(j+1)] = cx[2*j]; 	/* +1 because there are no DC */
		multisinex[2*(j+1)+1] = -cx[2*j+1]; /* components in the multisine */
	}
	realft(multisinex, numFFTpoints, -1);
	for (j=0; j<numFFTpoints; j++) multisinex[j] = multisinex[j]*2.0; /* Because ifft requires multiplying by 2/N 
																	     See NR documentation */

	/* DEBUG STUFF
	dbgfile2 = fopen("tempfft.dat","w");
	for (idbg=0; idbg<2048; idbg++) {
		dbgmsinex[idbg] = multisinex[idbg];
		fprintf(dbgfile2,"%f\n",(float)multisinex[idbg]);
	}
	fclose(dbgfile2);*/

	rms = 0.0;
	for (j=0; j<numFreq; j++) rms += cx[2*j]*cx[2*j] + cx[2*j+1]*cx[2*j+1];
	rms = sqrt(2.0*rms); /* factor 2 to account for complex conj Fourier amplitudes in real transform */
	
	multisinex[0] = multisinex[0]/rms; maxms = multisinex[0]; minms = maxms;
	for (j=1; j<numFFTpoints; j++) {
		multisinex[j] = multisinex[j]/rms;
		if (multisinex[j]>maxms) maxms = multisinex[j];
		if (multisinex[j]<minms) minms = multisinex[j];
	}
	*crx = (fabs(maxms)>fabs(minms)) ? fabs(maxms) : fabs(minms);
}

/**************************************************************************************************************/

double sign(double x)
{
	if (x<0)
		return -1.0;
	else if (x>0)
		return 1.0;
	else
		return 0.0;
}
