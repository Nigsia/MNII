#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double 	f		(double x, double y);
double* nabla_f	(double x, double y);
double g( double x, double y, double x0, double y0, double delta );
double* nabla_g ( double x, double y, double x0, double y0 ); 
double* comp_zero();
double* newton	(double x0, double y0, double delta);


int main( void )
{
	double *zero;

	zero=comp_zero();
	printf("%.16lf, %.16lf\n", zero[0], zero[1]);
	zero=newton(zero[0], zero[1], 0.01);
	printf("%.16lf, %.16lf\n", zero[0], zero[1]);

	return 0;	
}

double f( double x, double y )
{
	return (3*x*x+3*y*y-1)*(x*x+y*y-5)*(x*x+y*y-3*x+2)+1;
} 

double* nabla_f ( double x, double y )
{
	double *nf = (double*)malloc(2*sizeof(double));
	nf[0]=(6*x)*(x*x+y*y-5)*(x*x+y*y-3*x+2) + /*
	*/	  (3*x*x+3*y*y-1)*(2*x)*(x*x+y*y-3*x+2) + /*
	*/    (3*x*x+3*y*y-1)*(x*x+y*y-5)*(2*x-3);
	nf[1]=(6*y)*(x*x+y*y-5)*(x*x+y*y-3*x+2) + /*
	*/	  (3*x*x+3*y*y-1)*(2*y)*(x*x+y*y-3*x+2) + /*
	*/	  (3*x*x+3*y*y-1)*(x*x+y*y-5)*(2*y);
	return nf;
} 

double g( double x, double y, double x0, double y0, double delta )
{
	return pow(x-x0, 2) + pow(y-y0, 2) - pow(delta, 2);
}

double* nabla_g ( double x, double y, double x0, double y0 ) 
{
	double *ng = (double*)malloc(2*sizeof(double));
	ng[0] = 2*(x-x0);
	ng[1] = 2*(y-y0);
	return ng;
}

double* comp_zero()
{
	double y = 3.0;
	double *n_f;
	double *p = (double*)malloc(2*sizeof(double));

	while(fabs(f(0,y))>1e-10)
	{
		n_f = nabla_f(0,y);
		y = y - f(0,y)/n_f[1];
	}
	p[0] = 0.0;
	p[1] = y;

	free(n_f);

	return p;
}

double* newton ( double x0, double y0, double delta )
{
	double *xn, *h, *fxn;
	double **J;
	double dJ; 
	int i;

	J = (double**)malloc(2*sizeof(double*));
	for(i=0; i<2; i++)
		J[i] = (double*)malloc(2*sizeof(double));

	xn = (double*)malloc(2*sizeof(double));
	h = (double*)malloc(2*sizeof(double));
	fxn = (double*)malloc(2*sizeof(double));

	xn[0] = 0; 
	xn[1] = 2.5;

	i = 0;
	while(fabs(f(xn[0], xn[1]))>1e-10)
	{
		printf("------- Iter %d\n\n", i++);
		J[0] = nabla_f(xn[0],xn[1]);
		J[1] = nabla_g(xn[0],xn[1],x0,y0);

		printf("%+.4e %+.4e\n%+.4e %+.4e\n\n", J[0][0], J[0][1], J[1][0], J[1][1]);

		dJ = (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

		printf("%.10lf\n\n", dJ);

		fxn[0] = f(xn[0],xn[1]);
		fxn[1] = g(xn[0],xn[1], x0, y0, delta);

		/* Hem de resoldre el sistema lineal Df*h = -b */
		h[0] = ( (fxn[0] * J[1][1]) - (J[0][1] * fxn[1]) ) / dJ;
		h[1] = ( (J[0][0] * fxn[1]) - (fxn[0] * J[1][0]) ) / dJ;

		printf("%+.4e %+.4e\n\n", h[0], h[1]);

		xn[0] = xn[0] - h[0];
		xn[1] = xn[1] - h[1];

		printf("%+.4e %+.4e\n\n", xn[0], xn[1]);

	}

	return xn;

}