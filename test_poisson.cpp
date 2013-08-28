#include "slu_ddefs.h"
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <vector>
#include <string>
#include "lu_solver.h"
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;
#define EPS0 8.85e-12
#define PI 3.141592651
#define CTOM 1.758820088e11
#define EC 1.60217657e-19
#define EM 9.10938291e-31
#define AMTOEM 1822.887473
#define CL 2.99e8
/*
This file illustrates how to use a SuperLU solver for elliptical geometry.

*/
int main(int argc, char *argv[])
	{
	//GRID SIZE AND PARAMETERS...SIS100
	//PIPE TRANSVERSE DIMENSIONS
	double Rx = 0.02;
	double Ry = 0.02;
	//NUMBER OF CELLS+1 IN X AND Y DIRECTIONS
	int Nx=200;
	int Ny=200;
	//TOTAL NUMBER OF GRID POINTS
	int NxNy=Nx*Ny;

	//GRID DIMENSIONS
	double dx=2.0*Rx/(Nx-1);
	double dy=2.0*Ry/(Ny-1);

	//GENERATING DENSITY GRID
	double *rho_e=(double*) malloc(NxNy*sizeof(double));
	//TRANSVERSE ELECTRIC FIELDS
	double *Ex=(double*) malloc(NxNy*sizeof(double));
	double *Ey=(double*) malloc(NxNy*sizeof(double));
	//SETTING GRID VALUES TO ZERO
	set_zero(rho_e,NxNy);
	//DENSITY IS GIVEN AS A CHARGE PER CELL (PHYSICAL CHARGE DENSITY IS 1)
	for(int i=0;i<NxNy;i++)
		{
		rho_e[i]=1.0*dx*dy;
		}
	set_zero(Ex,NxNy);
	set_zero(Ey,NxNy);
	//CUT CELL PARAMETERS FOR THESE INDICES (EDGE LENGTH IN X AND Y RELATIVE TO CORRESPONDING dX OR dY)
	vector <double> alphax;
	vector <double> alphay;
	//STARTING THE BUNCH TRANSFER
	FILE *fd=fopen(string("potential.dat").c_str(),"w");
	FILE *fdx=fopen(string("Ex.dat").c_str(),"w");
	FILE *fdy=fopen(string("Ey.dat").c_str(),"w");

	int bordlength=2*Nx+2*Ny-4;
	//INDICES LYING INSIDE THE PHYSICAL DOMAIN
	vector <int> iindex;
	//INDICES INSIDE THE DOMAIN CLOSE TO BOUNDARY (i+-1 || i+-Ny ARE OUTSIDE)
	vector <int> bindex;
	//A STRUCTURE WITH ALL THE PARAMETERS NEEDED FOR SUPERLU MANIPULATIONS
	superlu_params Par;
	//PREPARE THE MATRIX AND ITS LU DECOMPOSITION AND SAVING THE LENGTH OF CUTCELLS
	preparation_of_matrix_superlu(Rx, Ry, Nx, Ny, dx, dy, Par, alphax, alphay, iindex, bindex);
	//SOLVING FOR THE POTENTIAL
	calc_potential_superlu(rho_e, &Par, Nx, Ny, dx, dy, iindex);
	//SOLVING FOR ELECTRIC FIELD
	//COEFFICIENT 1.0 MEANS 
	calc_field_cut_cell(rho_e, Ex, Ey, iindex, bindex, alphax, alphay,  Nx, Ny, dx, dy, 1.0);
	//SAVING CALCULATED VALUES
	for(int j=0;j<Ny;j++)
		{
		for(int i=0; i<Nx;i++)
			{
			fprintf(fd,"%e ",rho_e[j+i*Ny]);
			fprintf(fdx,"%e ",Ey[j+i*Ny]);
			fprintf(fdy,"%e ",Ex[j+i*Ny]);
			}
		fprintf(fd,"\n");
		fprintf(fdx,"\n");
		fprintf(fdy,"\n");
		}		
	fclose(fd);
	fclose(fdx);
	fclose(fdy);

	return 0;
	}


