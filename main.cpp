#include <stdio.h>
#include <math.h>
#include "chevalier.hpp"

int main(int argc, char **argv)
{
	Chevalier C;
	int		n_r = 1000;
	double	r_min = 1.0e-2;
	double	r_max = 1000.0*pow(10.,1.0); //r_max in pc
	double	r;

	double	M;
	double	u_star, u;
	double 	rho_star, rho;
	double	P_star, P;
  double  T;

	FILE *fp;
	char fname[200];

	r = 1.0e-2;
	printf("r %e P0 %e rho0 %e u0 %e T0 %e\n",r,C.Pressure(r)/C.kb_cgs,C.Density(r)/C.mp_in_grams,C.WindVelocity(r),C.Temperature(r));
	printf("r %e Mach %e u %e\n",r,C.MachNumber(r),C.WindVelocity(r));
	printf("x %e rho_star %e P_star %e u_star %e\n\n",r/C.R,C.rho_star(r),C.P_star(r),C.u_star(r)/(r/C.R));

	r = 300.;
	printf("r %e P0 %e rho0 %e u0 %e T0 %e\n",r,C.Pressure(r),C.Density(r),C.WindVelocity(r),C.Temperature(r));
	printf("r %e Mach %e u %e\n",r,C.MachNumber(r),C.WindVelocity(r));
	printf("x %e rho_star %e P_star %e u_star %e\n\n",r/C.R,C.rho_star(r),C.P_star(r),C.u_star(r)/(r/C.R));

	r = 500.;
	printf("r %e P0 %e rho0 %e u0 %e T0 %e\n",r,C.Pressure(r),C.Density(r),C.WindVelocity(r),C.Temperature(r));
	printf("r %e Mach %e u %e\n",r,C.MachNumber(r),C.WindVelocity(r));
	printf("x %e rho_star %e P_star %e u_star %e\n\n",r/C.R,C.rho_star(r),C.P_star(r),C.u_star(r)/(r/C.R));

	r = 700.;
	printf("r %e P0 %e rho0 %e u0 %e T0 %e\n",r,C.Pressure(r),C.Density(r),C.WindVelocity(r),C.Temperature(r));
	printf("r %e Mach %e u %e\n",r,C.MachNumber(r),C.WindVelocity(r));
	printf("x %e rho_star %e P_star %e u_star %e\n\n",r/C.R,C.rho_star(r),C.P_star(r),C.u_star(r)/(r/C.R));

	r = 1000.;
	printf("r %e P0 %e rho0 %e u0 %e T0 %e\n",r,C.Pressure(r)/C.kb_cgs,C.Density(r)/C.mp_in_grams,C.WindVelocity(r),C.Temperature(r));
	printf("r %e Mach %e u %e\n",r,C.MachNumber(r),C.WindVelocity(r));
	printf("x %e rho_star %e P_star %e u_star %e\n\n",r/C.R,C.rho_star(r),C.P_star(r),C.u_star(r)/(r/C.R));

	r = 10000.;
	printf("r %e P0 %e rho0 %e u0 %e T0 %e\n",r,C.Pressure(r),C.Density(r),C.WindVelocity(r),C.Temperature(r));
	printf("r %e Mach %e u %e\n",r,C.MachNumber(r),C.WindVelocity(r));
	printf("x %e rho_star %e P_star %e u_star %e\n",r/C.R,C.rho_star(r)*pow(r/C.R,2),C.P_star(r)*pow(r/C.R,10./3.),C.u_star(r));


	sprintf(fname,"chevalier.txt");

	fp = fopen(fname,"w");

	for(int i=0;i<n_r;i++)
	{
		r = (r_max-r_min)*((double) i)/((double) (n_r-1)) + r_min;
		M = C.MachNumber(r);
		u_star = C.u_star(r);
		P_star = C.P_star(r);
		rho_star = C.rho_star(r);
		fprintf(fp,"%e\t%e\t%e\t%e\t%e\n",r/C.R,M,u_star,P_star,rho_star);
	}
	fclose(fp);

	sprintf(fname,"chevalier.dimensionfull.txt");

	fp = fopen(fname,"w");

	for(int i=0;i<n_r;i++)
	{
		r = (r_max-r_min)*((double) i)/((double) (n_r-1)) + r_min;
		M = C.MachNumber(r);
		u = C.WindVelocity(r);
		P = C.Pressure(r);
		rho = C.Density(r);
    T = C.Temperature(r);
		fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\n",r,M,u,P,rho,T);
	}
	fclose(fp);
	//end
	return 0;
}
