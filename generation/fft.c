#include <math.h>
#include "fft.h"

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* The following functions are taken from Numerical Recipies in C */

void four1(double *data, unsigned long nn, int isign)
/* Replaces data[0..2*nn-1] by its discrete Fourier transform, if isign is input as 1; or replaces
   data[0..2*nn-1] by nn times its inverse discrete Fourier transform, if isign is input as -1.
   data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
   be an integer power of 2 (this is not checked for!).*/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine.*/
		if (j > i) {
			SWAP(data[j-1],data[i-1]); /* Exchange the two complex numbers. */
			SWAP(data[j],data[i]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	/* Here begins the Danielson-Lanczos section of the routine.*/
	mmax=2;
	while (n > mmax) { /* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; /* This is the Danielson-Lanczos formula: */
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence.*/
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void realft(double *data, unsigned long n, int isign)
/* Calculates the Fourier transform of a set of n real-valued data points. Replaces this data (which
   is stored in array data[1..n]) by the positive frequency half of its complex Fourier transform.
   The real-valued first and last components of the complex transform are returned as elements
   data[0] and data[1], respectively. n must be a power of 2. This routine also calculates the
   inverse transform of a complex data array if it is the transform of real data. (Result in this case
   must be multiplied by 2/n.) */
{
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;
	theta=3.141592653589793/(double) (n>>1); /* Initialize the recurrence. */
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1); /* The forward transform is here. */
	} else {
		c2=0.5; /* Otherwise set up for an inverse transform. */
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) { /* Case i=1 done separately below. */
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1-1]+data[i3-1]); /* The two separate transforms are separated out of data. */
		h1i=c1*(data[i2-1]-data[i4-1]);
		h2r = -c2*(data[i2-1]+data[i4-1]);
		h2i=c2*(data[i1-1]-data[i3-1]);
		data[i1-1]=h1r+wr*h2r-wi*h2i; /* Here they are recombined to form the true 
				                       transform of the original real data. */
		data[i2-1]=h1i+wr*h2i+wi*h2r;
		data[i3-1]=h1r-wr*h2r+wi*h2i;
		data[i4-1] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr; /* The recurrence. */
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[0] = (h1r=data[0])+data[1]; /* Squeeze the first and last data together
			                                to get them all within the original array. */
		data[1] = h1r-data[1];
	} 
	else {
		data[0]=c1*((h1r=data[0])+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1(data,n>>1,-1); /* This is the inverse transform for the case isign=-1.*/
	}
}