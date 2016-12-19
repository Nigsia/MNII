#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1e-10

double 	f			( double x, double y );
double* nabla_f		( double x, double y );
double 	g 			( double x, double y, double x0, double y0, double delta );
double* nabla_g 	( double x, double y, double x0, double y0 ); 
double* comp_zero	();
double* newton		( double x, double y, double x0, double y0, double delta );

int main( void )
{
	double *zeroInic, *zero;
	double delta = 0.01;
	double normDir;
	double *dir, *dir1;
	double aux;
	FILE *file;

	file=fopen("test", "w");

	zeroInic = comp_zero(); // Començem aquí
	fprintf(file,"%.14lf, %.14lf\n", zeroInic[0], zeroInic[1]);


	zero = (double*)malloc(2*sizeof(double));
	dir1 = (double*)malloc(2*sizeof(double));

	zero[0] = zeroInic[0];
	zero[1] = zeroInic[1];

	dir1[0] = 1;
	dir1[1] = 1;

	do
	{
		// obtenim una nova direcció
		dir = nabla_f( zero[0], zero[1] ); 
		// Projectem al tangent
		aux=dir[0];
		dir[0]=dir[1];
		dir[1]=-aux;
		// Normalitzar
		normDir = sqrt(dir[0]*dir[0] + dir[1]*dir[1]);
		dir[0] /= normDir;
		dir[1] /= normDir;
		
		if( dir[0]*dir1[0] + dir[1]*dir1[1] < 0 )
			zero = newton( zero[0]-dir[0]*delta, zero[1]-dir[1]*delta, zero[0], zero[1], delta );
		else
			zero = newton( zero[0]+dir[0]*delta, zero[1]+dir[1]*delta, zero[0], zero[1], delta );
		fprintf(file,"%.12lf, %.12lf\n", zero[0], zero[1]);

		dir1[0] = dir[0];
		dir1[1] = dir[1];

	}while( sqrt(pow(zeroInic[0]-zero[0],2)+pow(zeroInic[1]-zero[1],2)) > delta/2 );
	
	fclose(file);

	free(zeroInic);
	free(zero);
	free(dir);
	free(dir1);

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

	while(fabs(f(0,y))>TOL)
	{
		n_f = nabla_f(0,y);
		y = y - f(0,y)/n_f[1];
	}
	p[0] = 0.0;
	p[1] = y;

	free(n_f);

	return p;
}

double* newton ( double x, double y, double x0, double y0, double delta )
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

	xn[0] = x; 
	xn[1] = y;

	while(fabs(f(xn[0], xn[1]))>=TOL)
	{
		J[0] = nabla_f(xn[0],xn[1]);
		J[1] = nabla_g(xn[0],xn[1],x0,y0);

		dJ = (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

		if(fabs(dJ) < TOL )
		{
			printf("%+.4lf %+.4lf\n%+.4lf %+.4lf\n", J[0][0], J[0][1], J[1][0], J[1][1]);
			printf("El determinant és zero.\n");
			break;
		}

		fxn[0] = f(xn[0], xn[1]);
		fxn[1] = g(xn[0], xn[1], x0, y0, delta);

		h[0] = ( (fxn[0] * J[1][1]) - (J[0][1] * fxn[1]) ) / dJ;
		h[1] = ( (J[0][0] * fxn[1]) - (fxn[0] * J[1][0]) ) / dJ;

		xn[0] = xn[0] - h[0];
		xn[1] = xn[1] - h[1];
	}

	free(h);
	free(fxn);
	for(i=0;i<2;i++)
		free(J[i]);
	free(J);

	return xn;

}