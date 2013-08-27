#include "slu_ddefs.h"
#include <math.h>
#include <time.h>
#include <ctime>
#include <vector>
#include "supermatrix.h"
#include "lu_solver.h"

#define EPS0 8.85e-12
#define PI 3.1415926535897932384626433832795 //

using namespace std;


//SETTING GRID VALUES TO ZERO
void set_zero(double *rho, int NxNy)
	{
	for(int i=0;i<NxNy;++i)
		{
		rho[i]=0.0;
		}
	}

double derivative(double V1, double V2, double V3, double dL, int num)
	{
	//num=-1 - left point; num=0 - middle point; num=1 - right point
	return ((V1-V2*2+V3)*num+(V3-V1)/2.0)/dL;
	}

//PREPARATION OF LU DECOMPOSITION FOR POISSON SOLVER USING SUPERLU ROUTINES
void preparation_of_matrix_superlu(double Rx, double Ry,  long int Nx, long int Ny, double dx, double dy,
struct superlu_params &Par,
 vector <double> &alphax, vector <double> &alphay,
vector <int> &iindex,
vector <int> &bindex,
//ANGLE BETWEEN THE X AND Y AXIS AND THE WALL NEAR THE CUT CELL
vector <double> &tanx, vector <double> &tany)
	{
double dxdy=dx*dy;
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif
//     Defaults 
    Par.lwork = 0;
    Par.nrhs  = 1;
    Par.equil = YES;
    //Par.IterRefine = ::DOUBLE;	
    Par.u     = 0.1;
    Par.trans = NOTRANS;
    set_default_options(&Par.options);

    Par.options.Equil = Par.equil;
    Par.options.DiagPivotThresh = Par.u;
   // Par.options.IterRefine = SLU_DOUBLE;
    Par.options.Trans = Par.trans;
    Par.options.ColPerm=MMD_AT_PLUS_A;
    Par.options.SymmetricMode = YES;


    if ( Par.lwork > 0 ) {
	Par.work = SUPERLU_MALLOC(Par.lwork);
	if ( !Par.work ) {
	    ABORT("DLINSOLX: cannot allocate work[]");
	}
    }
printf("Creating arrays\n");
vector <double> Aa;
vector <int> Arows;
vector <int> Acolumn;
vector <int> incolumn;
vector <int> ix;
vector <int> iy;
//TEMPORARY VARIABLE FOR THE CUT-CELL INFORMATION
vector <double> alphatempx;
vector <double> alphatempy;
double Dx=dx*(Nx-1)/2.0;
double Dy=dy*(Ny-1)/2.0;
//HORIZONTAL INDEXING
Acolumn.push_back(0);

for(int i=0;i<Nx;++i)
	{
	int num=0;
	for(int j=0;j<Ny;++j)
		{
		double xx=(dx*i-Dx);
		double yy=(dy*j-Dy);
		double x=fabs(xx);
		double y=fabs(yy);
		//TAKING ONLY THE POINT INSIDE THE PIPE
		if(pow(x/Rx,2)+pow(y/Ry,2)<=1.0)
			{
			//ADDING INDICES TO TO VECTORS
			ix.push_back(i);
			iy.push_back(j);
			iindex.push_back(i*Ny+j);
			num++;
			//IF SHIFT ALONG X OR ALONG Y PER ONE CELL GETS US OUTSIDE THE PIPE
			if( (pow((x+dx)/Rx,2)+pow((y)/Ry,2)>1.0) || (pow((x)/Rx,2)+pow((y+dy)/Ry,2)>1.0) )
					{
					//FINDING THE HORIZONTAL CUT CELL EDGE LENGTH
					if(pow((x+dx)/Rx,2)+pow((y)/Ry,2)>=1.0)
						{
						double alphatemp=(Rx*sqrt(1.0-y*y/Ry/Ry)-x)/dx;
						/*if(alphatemp==1.0)
							{
							alphatemp=0.0;
							}*/
						alphax.push_back( alphatemp);
						//HERE THE  yy IS USED BECAUSE THE SIGN MATTERS
						tanx.push_back(-Rx*y/Ry/Ry/sqrt(1.0-y*y/Ry/Ry));
						//alphax.push_back(1.0);
						}
					else
						{
						alphax.push_back(1.0);
						}
					//FINDING THE VERTTICAL CUT CELL EDGE LENGTH
					if(pow((x)/Rx,2)+pow((y+dy)/Ry,2)>=1.0)
						{
						double alphatemp=(Ry*sqrt(1.0-x*x/Rx/Rx)-y)/dy;
						/*if(alphatemp==1.0)
							{
							alphatemp=0.0;
							}*/
						alphay.push_back( alphatemp);
						//alphay.push_back(1.0);
						//HERE THE  xx IS USED BECAUSE THE SIGN MATTERS
						tany.push_back(-Ry*x/Rx/Rx/sqrt(1.0-x*x/Rx/Rx));
						}
					else
						{
						alphay.push_back(1.0);
						}

					if(tanx.size()>tany.size())
						{
						//JUST SOME VALUE THAT WILL NOT BE USED
						tany.push_back(1.0);
						}
					else if(tany.size()>tanx.size())
						{
						//JUST SOME VALUE THAT WILL NOT BE USED
						tanx.push_back(1.0);
						}
					
					alphatempx.push_back(alphax.back());
					alphatempy.push_back(alphay.back());
					//bx.push_back(i);
					//by.push_back(j);	
					bindex.push_back(i*Ny+j);						
					}
				else
					{
					alphax.push_back(1.0);
					alphay.push_back(1.0);
					}
			}
		}
	if(num>0)
		{
		incolumn.push_back(num);
		}
	}

int totsize=ix.size();
int ixmin=ix.at(0);
double dxdx=dx*dx;
double dydy=dy*dy;
//CONSTRUCTING MATRIX X
//vector <double> matrix(totsize*totsize,0.0);
for(int i=0;i<totsize;++i)
	{
	int index=ix[i]-ixmin;
	//UPPER DIAGONAL
	if(index>0)
		{
		int deltam=(incolumn.at(index)+incolumn.at(index-1))/2;
		//IF SHIFT BY deltam PRESERVES THE ROW NUMBER
		if(i-deltam>=0 && iy[i]==iy.at(i-deltam))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i-deltam);
			//matrix.at(i*totsize+i-deltam)=Aa.back()/EPS0;
			}
		}
	//UP DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i>0 && ix[i]==ix.at(i-1))
		{
		//printf("Before up\n");
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i-1);
		//matrix.at(i*totsize+i-1)=Aa.back()/EPS0/dydy;
		}
	//MAIN DIAGONAL
		//Aa.push_back(2.0*EPS0*(1.0/dxdx/alphax[i]+1.0/dydy/alphay[i]));
		//Aa.push_back(EPS0*(2.0/dxdx/alphax[i]+2.0/dydy/alphay[i]) );
		Aa.push_back(EPS0*dxdy*((1.0+alphax[i])/alphax[i]/dxdx+(1.0+alphay[i])/alphay[i]/dydy) );
		Arows.push_back(i);
		//matrix.at(i*totsize+i)=Aa.back()/EPS0;
	//LOW DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i+1<totsize && ix[i]==ix.at(i+1))
		{
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i+1);
		//matrix.at(i*totsize+i+1)=Aa.back()/EPS0;
		}
	//LOWER DIAGONAL
	if(index+1<incolumn.size())
		{
		int deltap=(incolumn.at(index)+incolumn.at(index+1))/2;
		//IF SHIFT BY deltap PRESERVES THE ROW NUMBER
		if(i+deltap<totsize && iy[i]==iy.at(i+deltap))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i+deltap);
			//matrix.at(i*totsize+i+deltap)=Aa.back()/EPS0;
			}
		}
	Acolumn.push_back(Aa.size());
	}
	/*if(totsize<=10)
	{
	FILE *fd=fopen("matrice.dat","w");
	for(int i=0;i<totsize;i++)
		{
		for(int j=0;j<totsize;j++)
			{
			fprintf(fd,"%e ",matrix[i*totsize+j]);
			}
		fprintf(fd, "\n");
		}
	fclose(fd);
	}*/
	printf("Rewriting vectors to dynamic arrays\n");
	double *a=new double[Aa.size()];
	memcpy(a,&Aa[0],Aa.size()*sizeof(double));
	int *asub=new int[Arows.size()];
	memcpy(asub,&Arows[0],Arows.size()*sizeof(int));
	int *xa=new int[Acolumn.size()];
	memcpy(xa,&Acolumn[0],Acolumn.size()*sizeof(int));

Par.nnz=Aa.size();
Par.rhsb=new double[totsize];
Par.rhsx=new double[totsize];
for(int j=0;j<totsize;++j)
	{
	Par.rhsb[j]=1.0;
	Par.rhsx[j]=1.0;
	}
    Par.m=totsize;
    Par.n=totsize;
	printf("Creating matrices in SuperLU format\n");
    dCreate_CompCol_Matrix(&(Par.A), Par.m, Par.n, Par.nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.B), Par.m, Par.nrhs, Par.rhsb, Par.m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.X), Par.m, Par.nrhs, Par.rhsx, Par.m, SLU_DN, SLU_D, SLU_GE);
    Par.xact = doubleMalloc(Par.n * Par.nrhs);
    Par.ldx = Par.n;
    //dGenXtrue(n, nrhs, xact, ldx);
    //dFillRHS(trans, nrhs, xact, ldx, &A, &B);

	printf("Other SuperLU variables\n");
    if ( !(Par.etree = intMalloc(Par.n)) ) ABORT("Malloc fails for etree[].");
    if ( !(Par.perm_r = intMalloc(Par.m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(Par.perm_c = intMalloc(Par.n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(Par.R = (double *) SUPERLU_MALLOC(Par.A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(Par.C = (double *) SUPERLU_MALLOC(Par.A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(Par.ferr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(Par.berr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");
   //ONLY PERFORM THE LU DECOMPOSITION 
    Par.B.ncol = Par.nrhs;  //Indicate not to solve the system 
	printf("Finally LU decomposition\n");
    StatInit(&Par.stat);
    dgssvx(&(Par.options), &(Par.A), Par.perm_c, Par.perm_r, Par.etree, Par.equed, Par.R, Par.C,
           &(Par.L), &(Par.U), Par.work, Par.lwork, &(Par.B), &(Par.X), &(Par.rpg), &(Par.rcond), Par.ferr, Par.berr,
           &(Par.mem_usage), &(Par.stat), &(Par.info));

    StatFree(&(Par.stat));
    Par.options.Fact = FACTORED; // Indicate the factored form of A is supplied.
    Par.B.ncol = Par.nrhs;  // Set the number of right-hand side
	printf("Freeing stuff\n");
	delete [] xa;
	delete [] asub;
	delete [] a;
	printf("ix_size %d alpha_size %d\n",ix.size(),alphax.size());
	alphax=alphatempx;
	alphay=alphatempy;
	printf("New alphax_size %d\n",alphax.size());
	printf("Size of tangent vectors %d %d\n",tanx.size(),tany.size());
	}


//FOR RECTANGULAR GEOMETRY
void preparation_of_matrix_superlu_rectangular(long int Nx, long int Ny, double dx, double dy,
struct superlu_params &Par,
 vector <double> &alphax, vector <double> &alphay,
vector <int> &iindex,
vector <int> &bindex,
//ANGLE BETWEEN THE X AND Y AXIS AND THE WALL NEAR THE CUT CELL
vector <double> &tanx, vector <double> &tany)
	{
double dxdy=dx*dy;
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif
//     Defaults 
    Par.lwork = 0;
    Par.nrhs  = 1;
    Par.equil = YES;
    //Par.IterRefine = ::DOUBLE;	
    Par.u     = 0.1;
    Par.trans = NOTRANS;
    set_default_options(&Par.options);

    Par.options.Equil = Par.equil;
    Par.options.DiagPivotThresh = Par.u;
   // Par.options.IterRefine = SLU_DOUBLE;
    Par.options.Trans = Par.trans;
    Par.options.ColPerm=MMD_AT_PLUS_A;
    Par.options.SymmetricMode = YES;


    if ( Par.lwork > 0 ) {
	Par.work = SUPERLU_MALLOC(Par.lwork);
	if ( !Par.work ) {
	    ABORT("DLINSOLX: cannot allocate work[]");
	}
    }
printf("Creating arrays\n");
vector <double> Aa;
vector <int> Arows;
vector <int> Acolumn;
vector <int> incolumn;
vector <int> ix;
vector <int> iy;
//TEMPORARY VARIABLE FOR THE CUT-CELL INFORMATION
vector <double> alphatempx;
vector <double> alphatempy;
double Dx=dx*(Nx-1)/2.0;
double Dy=dy*(Ny-1)/2.0;
//HORIZONTAL INDEXING
Acolumn.push_back(0);

for(int i=1;i<Nx-1;++i)
	{
	int num=0;
	for(int j=1;j<Ny-1;++j)
		{
		double xx=(dx*i-Dx);
		double yy=(dy*j-Dy);
		double x=fabs(xx);
		double y=fabs(yy);
		//TAKING ONLY THE POINT INSIDE THE PIPE

			//ADDING INDICES TO TO VECTORS
			ix.push_back(i);
			iy.push_back(j);
			iindex.push_back(i*Ny+j);
			num++;
			//IF SHIFT ALONG X OR ALONG Y PER ONE CELL GETS US OUTSIDE THE PIPE
			{
			alphax.push_back(1.0);
			alphay.push_back(1.0);
			}
			
		}
	if(num>0)
		{
		incolumn.push_back(num);
		}
	}

int totsize=ix.size();
int ixmin=ix.at(0);
double dxdx=dx*dx;
double dydy=dy*dy;
//CONSTRUCTING MATRIX X
//vector <double> matrix(totsize*totsize,0.0);
for(int i=0;i<totsize;++i)
	{
	int index=ix[i]-ixmin;
	//UPPER DIAGONAL
	if(index>0)
		{
		int deltam=(incolumn.at(index)+incolumn.at(index-1))/2;
		//IF SHIFT BY deltam PRESERVES THE ROW NUMBER
		if(i-deltam>=0 && iy[i]==iy.at(i-deltam))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i-deltam);
			//matrix.at(i*totsize+i-deltam)=Aa.back()/EPS0;
			}
		}
	//UP DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i>0 && ix[i]==ix.at(i-1))
		{
		//printf("Before up\n");
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i-1);
		//matrix.at(i*totsize+i-1)=Aa.back()/EPS0/dydy;
		}
	//MAIN DIAGONAL
		//Aa.push_back(2.0*EPS0*(1.0/dxdx/alphax[i]+1.0/dydy/alphay[i]));
		//Aa.push_back(EPS0*(2.0/dxdx/alphax[i]+2.0/dydy/alphay[i]) );
		Aa.push_back(EPS0*dxdy*((1.0+alphax[i])/alphax[i]/dxdx+(1.0+alphay[i])/alphay[i]/dydy) );
		Arows.push_back(i);
		//matrix.at(i*totsize+i)=Aa.back()/EPS0;
	//LOW DIAGONAL//IF SHIFT BY ONE ELEMENT PRESERVES THE COLUMN
	if(i+1<totsize && ix[i]==ix.at(i+1))
		{
		Aa.push_back(-EPS0*dxdy/dydy/*/(1.0+alphay[i])*/);
		Arows.push_back(i+1);
		//matrix.at(i*totsize+i+1)=Aa.back()/EPS0;
		}
	//LOWER DIAGONAL
	if(index+1<incolumn.size())
		{
		int deltap=(incolumn.at(index)+incolumn.at(index+1))/2;
		//IF SHIFT BY deltap PRESERVES THE ROW NUMBER
		if(i+deltap<totsize && iy[i]==iy.at(i+deltap))
			{
			Aa.push_back(-EPS0*dxdy/dxdx/*/(1.0+alphax[i])*/);
			Arows.push_back(i+deltap);
			//matrix.at(i*totsize+i+deltap)=Aa.back()/EPS0;
			}
		}
	Acolumn.push_back(Aa.size());
	}
	/*if(totsize<=10)
	{
	FILE *fd=fopen("matrice.dat","w");
	for(int i=0;i<totsize;i++)
		{
		for(int j=0;j<totsize;j++)
			{
			fprintf(fd,"%e ",matrix[i*totsize+j]);
			}
		fprintf(fd, "\n");
		}
	fclose(fd);
	}*/
	printf("Rewriting vectors to dynamic arrays\n");
	double *a=new double[Aa.size()];
	memcpy(a,&Aa[0],Aa.size()*sizeof(double));
	int *asub=new int[Arows.size()];
	memcpy(asub,&Arows[0],Arows.size()*sizeof(int));
	int *xa=new int[Acolumn.size()];
	memcpy(xa,&Acolumn[0],Acolumn.size()*sizeof(int));

Par.nnz=Aa.size();
Par.rhsb=new double[totsize];
Par.rhsx=new double[totsize];
for(int j=0;j<totsize;++j)
	{
	Par.rhsb[j]=1.0;
	Par.rhsx[j]=1.0;
	}
    Par.m=totsize;
    Par.n=totsize;
	printf("Creating matrices in SuperLU format\n");
    dCreate_CompCol_Matrix(&(Par.A), Par.m, Par.n, Par.nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.B), Par.m, Par.nrhs, Par.rhsb, Par.m, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&(Par.X), Par.m, Par.nrhs, Par.rhsx, Par.m, SLU_DN, SLU_D, SLU_GE);
    Par.xact = doubleMalloc(Par.n * Par.nrhs);
    Par.ldx = Par.n;
    //dGenXtrue(n, nrhs, xact, ldx);
    //dFillRHS(trans, nrhs, xact, ldx, &A, &B);

	printf("Other SuperLU variables\n");
    if ( !(Par.etree = intMalloc(Par.n)) ) ABORT("Malloc fails for etree[].");
    if ( !(Par.perm_r = intMalloc(Par.m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(Par.perm_c = intMalloc(Par.n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(Par.R = (double *) SUPERLU_MALLOC(Par.A.nrow * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(Par.C = (double *) SUPERLU_MALLOC(Par.A.ncol * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(Par.ferr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(Par.berr = (double *) SUPERLU_MALLOC(Par.nrhs * sizeof(double))) ) 
        ABORT("SUPERLU_MALLOC fails for berr[].");
   //ONLY PERFORM THE LU DECOMPOSITION 
    Par.B.ncol = Par.nrhs;  //Indicate not to solve the system 
	printf("Finally LU decomposition\n");
    StatInit(&Par.stat);
    dgssvx(&(Par.options), &(Par.A), Par.perm_c, Par.perm_r, Par.etree, Par.equed, Par.R, Par.C,
           &(Par.L), &(Par.U), Par.work, Par.lwork, &(Par.B), &(Par.X), &(Par.rpg), &(Par.rcond), Par.ferr, Par.berr,
           &(Par.mem_usage), &(Par.stat), &(Par.info));

    StatFree(&(Par.stat));
    Par.options.Fact = FACTORED; // Indicate the factored form of A is supplied.
    Par.B.ncol = Par.nrhs;  // Set the number of right-hand side
	printf("Freeing stuff\n");
	delete [] xa;
	delete [] asub;
	delete [] a;
	printf("ix_size %d alpha_size %d\n",ix.size(),alphax.size());
	alphax=alphatempx;
	alphay=alphatempy;
	printf("New alphax_size %d\n",alphax.size());
	printf("Size of tangent vectors %d %d\n",tanx.size(),tany.size());
	}

void calc_potential_superlu(double *rho, struct superlu_params *Par, long int Nx, long int Ny, double dx, double dy,
vector <int> iindex)
	{
	int totsize=iindex.size();
	for(int i=0;i<totsize;++i)
		{
		Par->rhsb[i]=rho[iindex[i]];
		Par->rhsx[i]=rho[iindex[i]];
		}
	dCreate_Dense_Matrix(&Par->B, Par->m, Par->nrhs, Par->rhsb, Par->m, SLU_DN, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&Par->X, Par->m, Par->nrhs, Par->rhsx, Par->m, SLU_DN, SLU_D, SLU_GE);
    StatInit(&(Par->stat));
    dgssvx(&(Par->options), &(Par->A), Par->perm_c, Par->perm_r, Par->etree, Par->equed, Par->R, Par->C,
           &(Par->L), &(Par->U), Par->work, Par->lwork, &(Par->B), &(Par->X), &(Par->rpg), &(Par->rcond), Par->ferr, Par->berr,
           &(Par->mem_usage), &(Par->stat), &(Par->info));
    StatFree(&(Par->stat));
	double *data=(double*) ((DNformat*) Par->X.Store)->nzval;
	
	set_zero(rho,Nx*Ny);
	for(int i=0;i<totsize;++i)
		{
		rho[iindex[i]]=data[i];
		}
	//free(data);
	}


//ELECTRIC FIELD SOLUTION. GRADIENT OF CALCULATED POTENTIAL
//ctomdt IS 1.0 IF YOU WANT TO CALCULATE THE FIELD, IF DIRECTLY THE ACCELERATION THEN ELECTRON CHARGE TO MASS
void calc_field(double *phi, double *Ex, double *Ey, int Nx, int Ny, double dx, double dy, double ctomdt)
	{
	int flagx=0;
	int flagy=0;
	int globind=0;
	for(int i=0;i<Nx;++i)
		{
		for(int j=0;j<Ny;++j)
			{
			//GLOBAL INDEX
			globind=i+Nx*j;
			//DERIVATIVE FLAG FOR X DIRECTION
			if(i==0)
				{
				flagx=-1;
				}
			else if(i==Nx-1)
				{
				flagx=1;
				}	
			else
				{
				flagx=0;
				}
			//DERIVATIVE FLAG FOR Y DIRECTION
			if(j==0)
				{
				flagy=-1;
				}
			else if(j==Ny-1)
				{
				flagy=1;
				}	
			else
				{
				flagy=0;
				}
			//TAKING DERIVATIVE OF POTENTIAL
			//Ex[globind]=ctomdt*derivative(phi[globind-1-flagx],phi[globind-flagx],phi[globind+1-flagx],dx,flagx);
			//Ey[globind]=ctomdt*derivative(phi[globind-(1+flagy)*Nx],phi[globind-flagy*Nx],phi[globind+(1-flagy)*Nx],dy,flagy);
			Ex[globind]=ctomdt*derivative(phi[globind-(1+flagx)*Ny],phi[globind-flagx*Ny],phi[globind+(1-flagx)*Ny],dx,flagx);
			Ey[globind]=ctomdt*derivative(phi[globind-1-flagy],phi[globind-flagy],phi[globind+1-flagy],dy,flagy);
			}
		}
	};
//ELECTRIC FIELD SOLUTION. GRADIENT OF CALCULATED POTENTIAL. TAKES CUT CELL INFORMATION TO CALCULATE THE CORRECT FIELD AT THE BOUNDARY AND EXTRAPOLATE IT ONE CELL AWAY FROM THE BOUNDARY

void calc_field_cut_cell(
//FIELD ARRAYS
double *phi, double *Ex, double *Ey,
//VECTOR OF INTERNAL GRID COORDINATES
vector <int> iindex,
//VECTOR OF POINTS NEAR THE BOUNDARY
vector <int> bindex,
//VECTORS OF RELATIVE CUT-CELL EDGES
vector <double> alphax, vector <double> alphay,
//SIZE AND STEP SIZE OF THE GRID
 int Nx, int Ny, double dx, double dy, 
double ctomdt)
	{
	double Dx=dx*(Nx-1)/2;
	double Dy=dy*(Ny-1)/2;
	int flagx=0;
	int flagy=0;
	int globind=0;
	int totsize=iindex.size();
	double dxdx=dx*dx;
	double dydy=dy*dy;
	double constx=ctomdt/2.0/dx;
	double consty=ctomdt/2.0/dy;
	for(int i=0;i<totsize;++i)
		{
		globind=iindex[i];
		//INSIDE OF PIPE
		if(phi[globind+Ny]!=0.0 && phi[globind-Ny]!=0.0)
			{
			Ex[globind]=constx*(phi[globind+Ny]-phi[globind-Ny]);
			}
		if(phi[globind+1]!=0.0 && phi[globind-1]!=0.0)
			{
			Ey[globind]=consty*(phi[globind+1]-phi[globind-1]);
			}
		}

	int bordersize=bindex.size();
	//FIRST TIME WE CALCULATE THE SIMPLEST DERIVATIVES
	for(int i=0;i<bordersize;++i)
		{
		globind=bindex[i];
		int by=globind%Ny;
		int bx=globind/Ny;
		//AT THE RIGHT SIDE OF A PIPE
		if(phi[globind+Ny]==0.0)
			{
			double E1=(phi[globind]-phi[globind-Ny])/dx;
			double E2=-phi[globind]/alphax[i]/dx;
			if(alphax[i]==1.0)
				{			
				printf("%e %e %e\n",alphax[i],(dx*bx-Dx)/0.0609, (dy*by-Dy)/0.02425);
				}
			Ex[globind]=ctomdt*(E1+E2*alphax[i])/(1.0+alphax[i]);
			//FIELD AT THE BOUNDARY
			//double E0x=2.0*E2-Ex[globind];
			Ex[globind+Ny]=2.0*Ex[globind]-Ex[globind-Ny];
			}
		//AT THE LEFT SIDE OF A PIPE
		if(phi[globind-Ny]==0.0)
			{
			double E1=(phi[globind+Ny]-phi[globind])/dx;
			double E2=phi[globind]/alphax[i]/dx;
			Ex[globind]=ctomdt*(E1+E2*alphax[i])/(1.0+alphax[i]);
			Ex[globind-Ny]=2.0*Ex[globind]-Ex[globind+Ny];
			}
		//AT THE UPPER SIDE OF A PIPE
		if(phi[globind+1]==0.0)
			{
			double E1=(phi[globind]-phi[globind-1])/dy;
			double E2=-phi[globind]/alphay[i]/dy;
			Ey[globind]=ctomdt*(E1+E2*alphay[i])/(1.0+alphay[i]);
			Ey[globind+1]=2.0*Ey[globind]-Ey[globind-1];
			}
		//AT THE LOWER SIDE OF A PIPE
		if(phi[globind-1]==0.0)
			{
			double E1=(phi[globind+1]-phi[globind])/dy;
			double E2=phi[globind]/alphay[i]/dy;
			Ey[globind]=ctomdt*(E1+E2*alphay[i])/(1.0+alphay[i]);
			Ey[globind-1]=2.0*Ey[globind]-Ey[globind+1];
			}
		}
	//SECOND TIME WE INTERPOLATE TRANSVERSE DERIVATIVES (NOT USING ANGLES AT A TIME)
	for(int i=0;i<bordersize;++i)
		{
		globind=bindex[i];
		//int by=globind%Ny;
		//int bx=globind/Ny;
		if(phi[globind+Ny]==0.0 && Ey[globind+Ny]==0.0)
			{
			//WE INTERPOLATE OUTSIDE THE GRID
			Ey[globind+Ny]=2.0*Ey[globind]-Ey[globind-Ny];
			}
		if(phi[globind+1]==0.0 && Ex[globind+1]==0.0)
			{

			Ex[globind+1]=2.0*Ex[globind]-Ex[globind-1];
			}
		if(phi[globind-Ny]==0.0 && Ey[globind-Ny]==0.0)
			{
			//WE CALCULATE THE Ex ON THE BORDER
			Ey[globind-Ny]=2.0*Ey[globind]-Ey[globind+Ny];
			}
		if(phi[globind-1]==0.0 && Ex[globind-1]==0.0)
			{

			Ex[globind-1]=2.0*Ex[globind]-Ex[globind+1];
			}
		}

	};







