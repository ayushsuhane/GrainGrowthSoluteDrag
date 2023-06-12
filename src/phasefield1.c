#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<omp.h>
#include"Global.h"
#include"functions.h"
//#include "gsl/gsl_matrix_double.h"
//#include "gsl/gsl_linalg.h"
//#include "gsl/gsl_blas.h"



void pfsolver(int t)
{

	extern int nx, ny, nx, dim, grid2, grid3, ngrainsmax;
	extern float *grid, *newgrid, *curvature;
	extern float *dphidtold;
	extern int *flaggrid, *newflaggrid;
	extern float **orientation, **rotationmatrix;
	extern int bc, N, M;
	extern float dt, mob, eps, omega;

	extern float c_length, c_time;
	extern float threshold, p_threshold;
	extern float E0, T;

	extern float  mu_anisotropy, sigma_anisotropy, sigma_e0anisotropy;

	extern int flag_anisotropy, flag_e0anisotropy;

	extern int trackid, flag_curvature;
	extern int t_output;
	extern int NUM_THREADS;





	int left, right, top, down, up, bottom;
	int iter, counter;
	float sumphi2=0, sumphi;
	float misorientation = 45;
	
	int i=0, j=0, k=0, l=0, id, m;
	float grad;
	int index, count, gind;
	float velocity, velocity_prev, dphidt, dv, vel;
	float gradphi[3];
	float constantrate=0.5, cr;


	float dphidtpast, error, mag_gradphi, dummy;
	//gradphi = malloc(3 * sizeof(float));
	for(int i=0; i<3; i++) gradphi[i]=0;

	if(flag_curvature && (t%t_output)==0) {for(i=0; i<N; i++) curvature[i] = 0.0;}

	//printf("First : Flaggrid and Grid!!");
	//for(l=0; l<ngrainsmax; l++) printf("index : %d, flaggrid: %d phi: %lf\n", pos2ind(trackid, 0, 0), flaggrid[l*grid3 + pos2ind(trackid, 0, 0)], grid[l*grid3 + pos2ind(trackid, 0, 0)]);
	if(t==0) updatefornextstep();

	//printf("t : %d\n",t);
	
	/*Loop over all the gridpoints*/
	#pragma omp parallel for collapse(3) num_threads(NUM_THREADS) schedule(auto) private(left, right, top, bottom, up, down, gradphi, count, sumphi2, sumphi, index, gind, id, l, grad, misorientation, iter, dphidt, velocity, velocity_prev, dv, vel, cr, m, error, dphidtpast, mag_gradphi, counter, dummy) shared(t, flaggrid, newflaggrid, grid, newgrid, curvature, orientation, rotationmatrix, dphidtold)
	for(i=0; i<nx; i++)
	{
	  for(j=0; j<ny; j++)
	  {
	  	for(k=0; k<nz; k++)
	  	{
	  	  //printinfo_gridpoint(i, j, k);
	  	  count=0; sumphi2=0;
	  	  sumphi=0;
	      index = k*grid2 + j*nx + i;
	      //Only solve when more than 2 active grains are present
	  	  //if(i==trackid && j==0 && k==0){  for(l=0; l<ngrainsmax;l++) printf("flaggrid: %d\n", flaggrid[l*grid3+index]);}
	  	  for(l=0; l<ngrainsmax;l++){if(flaggrid[index + l*grid3]==-1) break; count+=1;}
	  	  if(t>2) {if(count<2) continue;}
	  	  //printf("index: %d, count : %d\n", index, count);

	  	  //printf("t = %d, i : %d, j:%d, k:%d index : %d, Nthread : %d count : %d\n", t, i, j, k, index, omp_get_thread_num(), count);
	      
	  	  // Evaluate the sum of squares of phi
	  	  for(l=0; l<ngrainsmax; l++){ if(flaggrid[index + l*grid3]==-1) break; sumphi2 += grid[index + l*grid3]*grid[index + l*grid3];sumphi += grid[index + l*grid3];}	

	  	  left = -1; right = -1;  top=-1;  bottom=-1; up=-1; down=-1; //set everything to -1
	  	  // If no flux condition, only solve for internal nodes - also considers 2d and 3d system
	  	  if(bc==0 || bc==2){id = boundarycondition(i, j, k, dim);if(!id) continue;}
          
          // Assign the neighbouring index considering boundary conditions
	  	  left = pos2ind(i-1, j, k); right = pos2ind(i+1, j, k);
	  	  if(dim>1)
	  	  {
	  	  	top  = pos2ind(i, j+1, k);bottom = pos2ind(i, j-1, k);
	  	  }
	  	  if(dim>2)
	  	  {
	  	  	up = pos2ind(i, j, k+1);down =  pos2ind(i, j, k-1);
	  	  }

	  	  /***Misorientation calculations *****/
	  	  /***Only calculate when anisotropic flag is activated ****/
	  	  if((flag_anisotropy||flag_e0anisotropy) && t>2)
	  	  {
	  	  	int maxindex=0, maxindex1=1;
	        float maxphi=0.0, maxphi1=0.0;
	        int grain1 = maxindex, grain2 = maxindex1;
	        for(l=0; l<ngrainsmax;l++)
	        {
		      if(grid[index+l*grid3]>maxphi) {maxindex=index+l*grid3; maxphi=grid[maxindex]; grain1 = flaggrid[maxindex];}
	        }
	        for(l=0; l<ngrainsmax;l++)
	        {
		      if(grid[index+l*grid3]>maxphi1 && ((index+l*grid3)!=maxindex)) {maxindex1=index+l*grid3; maxphi1=grid[maxindex1]; grain2=flaggrid[maxindex1];}
	        }
	        /*
	        for(l=0; l<ngrainsmax; l++) printf("l : %d, Grainnum : %d, phi : %lf\n", l, flaggrid[index+l*grid3], grid[index+l*grid3]);
	        printf("i : %d, j: %d, k: %d, grain1 : %d, grain2 : %d\n", i, j, k, grain1, grain2);
	        */
	        //printf("i : %d, j: %d, k: %d, grain1 : %d, grain2 : %d, orientation[grain1] : %lf, orientation[grain2] : %lf\n", i, j, k, grain1, grain2, orientation[grain1], orientation[grain2]);
	        misorientation = evaluatemisorientation(grain1, grain2);
	        //if(i==200 && j<ny/2) printf("%d %d, Grain Info: %d, %d, Misorientation : %lf, phi: %le\n", i, j, grain1, grain2, misorientation, grid[maxindex1]);
	        //misorientation = fabs(orientation[grain1] - orientation[grain2]);
	  	  }
	  	  else misorientation = mu_anisotropy;

	  	  /********/

	  	  
	  	  //if(i==trackid && j==0 && k==0){  for(l=0; l<ngrainsmax;l++) printf("flaggrid: %d\n", flaggrid[l*grid3+index]);}
	  	  //for(l=0; l<ngrainsmax; l++) {printf("index : %d l:%d flag:%d ", index, l, flaggrid[l*grid3 + index]);}; printf("\n");

	  	  for(l=0; l<ngrainsmax; l++)
	  	  {

	  	  	gind = index + l*grid3;
	  	  	if(flaggrid[gind]==-1) { break;}
	  	  	//if(i==trackid && j==0 && k==0){ for(l=0; l<ngrainsmax;l++) printf("flaggrid: %d\n", flaggrid[l*grid3+index]);}
	  	  	// While forming grad for the next timestep take care to only keep all the grains with non-zero gradients only elimintate phi with 0 at the end
	  	  	grad = calculatelaplacian(index, left, right, top, bottom, up, down, dim, flaggrid[gind]); //Grain number at the end  
	  	  	if(flag_curvature && (t%(t_output)==0) && t>2) {/*printf("i : %d, j:%d k:%d phi: %lf\n", i, j, k, grid[gind]);*/curvature[gind]=calculatecurvature(grid[gind], grad);}


	  	  	newgrid[gind] = grid[gind] + dt*(anisotropic_mobility(mob, misorientation, flag_anisotropy))*(eps*eps*grad - omega*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 3*grid[gind]*(sumphi2-grid[gind]*grid[gind])));
	  	  	//newgrid[gind] = grid[gind] + dt*(mob)*(eps*eps*grad - omega*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 3*grid[gind]*(sumphi2)));
	  	  	
	  	  	cr = constantrate;
	  	  	dphidt = (newgrid[gind]-grid[gind])/dt;

	  	  	   
	  	  	if(flag_sd && (grid[gind]>0.001) && (grid[gind]<(1-0.001)) && (t>100))
	  	  	{
	  	  	  counter=0;
	  	  	  error = 1;
	  	  	  iter=0;
	  	  	  for(m=0; m<3; m++) gradphi[m]=0;
	  	  	  
	  	  	  gradientvector(index, left, right, top, bottom, up, down, dim, flaggrid[gind], gradphi);
	  	  	  mag_gradphi = sqrt(gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1]);
	  	  	  //dphidtpastt = 0;
	  	  	  dummy = dphidtold[gind];

	  	  	  //dphidtpast = dphidt;
	  	  	  //printf("dphidt : %lf\n", dphidt);
	  	  	  //if(i==(int)(nx/2)) printf("Before entering SD loop : \n t : %d, i : %d,  j: %d, l: %d dphidt : %lf, grid : %le, newgrid : %le\n", t, i, j, l, dphidt, grid[gind], newgrid[gind]);
	  	  	  cr = 0.5;
	  	  	  while((error > 0.05) && (iter<800))
	  	  	  {
	  	  	  	
	  	  	  	//dphidtpast = dphidtold[gind];
	  	  	  	dphidtpast = dphidt;
	  	  	  	//if(counter==0) dphidtpast = dphidtold[gind];
	  	  	  	if(counter==0) dphidtpast = dummy;
	  	  	  	velocity = findvelocity(gradphi, dphidtpast, dim);
	  	  	  	if(fabs(mag_gradphi)<1e-5) {velocity=0; break;}
	  	  	  	
	  	  	  	
	  	  	  	dphidt = anisotropic_mobility(mob, misorientation, flag_anisotropy)*(eps*eps*grad - omega*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 3*grid[gind]*(sumphi2-grid[gind]*grid[gind])) + 3*dG_drag(anisotropic_segregation(E0, misorientation, flag_e0anisotropy), c_init, velocity*c_length/c_time, T)*grid[gind]*(sumphi-grid[gind]));
	  	  	  	error = fabs((dphidtpast-dphidt)/dphidtpast);
	  	  	  	dphidt = dphidtpast + cr*(dphidt - dphidtpast);	
	  	  	  	//dphidtold[gind] = dphidt;
	  	  	  	counter+=1;
				//if(i==(int)(nx/2)) printf("index : %d, Velocity : %lf dphidt : %lf, dphidtpast : %lf, iter : %d grad: %le mag_gradphi : %le, error = %le\n",index, velocity, dphidt, dphidtpast, iter, grad, sqrt(gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1]), error);
	  	  	  	/*if((dphidt*dphidtpast<0))*/ {if(iter%40==0 && iter<400) {velocity=0;dummy/=10;counter=0;cr/=10;}/*dphidt=0.0; break;*/}
	  	  	  	if(iter==400) {velocity=0;dummy*=1e4;counter=0;cr=0.5;}
	  	  	  	{if(iter%40==0 && iter>400) {velocity=0;dummy*=10;counter=0;cr/=2;}}
	  	  	  	iter += 1;

				//if(iter==500) {dphidt = fmax(dphidt, dphidtpast);}
				if(iter==400) printf("THIS IS NOT WORKING YAAAR, counter: %d error : %le, i: %d, j:%d, l:%d t:%d\n", counter,error, i, j, l, t);
				//if(i==8 && j==89) printf("t : %d, i : %d, j: %d, l : %d, index : %d, Velocity : %lf dphidt : %lf, dphidtpast : %lf, iter : %d grad: %le mag_gradphi : %le, error = %le, grid : %le\n",t, i, j, l, index, velocity, dphidt, dphidtpast, iter, grad, sqrt(gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1]), error, grid[gind]);
						//printf("iter = %d, dphidt = %le, dphidtpast = %le, error = %le\n",iter, dphidt, dphidtpast, error);
				//}
	  	  	  }

	  	  	  //if(iter>1) printf("index : %d, iter : %d, dphidt : %le, newgrid : %le\n", gind, iter, dphidt, newgrid[gind]);
	  	  	  newgrid[gind] = grid[gind] + dt*dphidt;
	  	  	  if(newgrid[gind]<0.0) newgrid[gind]=0.0;
	  	  	  else if(newgrid[gind]>1.0) newgrid[gind]=1.0;

/*	  	  	
	  	  	  iter = 0;
	  	  	  dphidt = (newgrid[gind]-grid[gind])/dt;
	  	  	  for(m=0; m<3; m++) gradphi[m]=0;
	  	  	  gradphi = gradientvector(index, right, top, up, dim, flaggrid[gind]);
	  	  	  velocity = findvelocity(gradphi, dphidt, dim);
	  	  	  dv = velocity; velocity_prev = velocity;
	  	  	  

	  	  	
	  	  	 if(fabs(velocity)>1e-14 && grid[gind]>0.001)
	  	  	 {

	  	  	  	while((abs(dv/velocity_prev) > 0.01) && iter<501)
	  	  	    { 

	  	  	  	  newgrid[gind] = grid[gind] + dt*(anisotropic_mobility(mob, misorientation, flag_anisotropy))*(eps*eps*grad - 
	  	  	  	  omega*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 3*grid[gind]*(sumphi2-grid[gind]*grid[gind]))
	  	  	  	  + 3*dG_drag2(anisotropic_segregation(E0, misorientation, flag_e0anisotropy), c_init, velocity*c_length/c_time, T)*grid[gind]*(sumphi-grid[gind]));
	  	  	  	  
	  	  	  	  if(newgrid[gind]<0.0) newgrid[gind]=0.0;
	  	  		  else if(newgrid[gind]>1.0) newgrid[gind]=1.0;

	  	  	  	  dphidt = (newgrid[gind]-grid[gind])/dt;
	  	  	  	  //gradphi = gradientvector(index, right, top, up, dim, flaggrid[gind]);
	  	  	  	  vel = findvelocity(gradphi, dphidt, dim);
	  	  	  	  velocity = velocity_prev + cr*(vel-velocity_prev);
	  	  	  	  dv = velocity-velocity_prev; 
	  	  	  	  if(iter==500) {printf("THIS IS NOT WORKING AT ALL\n"); printf("Index : %d, velocity : %le, previous_velocity : %le, dv = %le phi : %le grad : %le Eseg : %le ActualVelocity : %le cr : %le\n", index, velocity, velocity_prev, dv, grid[gind], grad, anisotropic_segregation(E0, misorientation, flag_e0anisotropy), velocity*c_length/c_time, cr);}
	  	  	  	  
	  	  	  	  //Put velocity 0 and move on if the sign of velocity changes
	  	  	  	  iter += 1;

	  	  	  	  if((((iter)%(50))==0)) 
	  	  	  	  	{cr=cr*0.5;} 
	  	  	  	  velocity_prev=velocity;  
	  	  	  	  
	  	  	    }
	  	  	  }
*/
	  	  	}
	  	  	dphidtold[gind]=dphidt;
	  	  	if(newgrid[gind]<0.0) newgrid[gind]=0.0;
	  	  	else if(newgrid[gind]>1.0) newgrid[gind]=1.0;

	  	  	  //if(i==nx/2 && t>100) printf("iter : %d velocity : %e\n", iter, velocity*c_length/c_time);
	  	  	//}
	  	  	//if(((i==trackid) || (i==trackid-1) || (i==trackid+1)) && j==0 && k==0) printf("index : %d, l: %d, grainnum : %d, grad: %lf, val: %lf, newgrid:%lf\n", index, l, flaggrid[gind], grad, dt*mob*(eps*eps*grad - omega*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 3*grid[gind]*(sumphi2-grid[gind]*grid[gind]))), newgrid[gind]);
	  	  	newflaggrid[gind] = flaggrid[gind];
	  	  	//printinfo_gridpoint(i, j, k);  	  	
	  	  }
	  	}
	  }
	}
	updatefornextstep();
	//free(gradphi);
	return;
}

int pos2ind(int i, int j, int k)
{
  extern int nx, ny, nz, grid2;
  extern int bc;
  if(bc!=0)
  {
  	if(i<0) i=i+nx;
    if(j<0) j=j+ny;
    if(k<0) k=k+nz;
    if(i>nx-1) i=i-nx;
    if(j>ny-1) j=j-ny;
    if(k>nz-1) k=k-nz;
  }
  return k*grid2 + j*nx + i;
}

int boundarycondition(int i, int j, int k, int dim)
{
	extern int nx, ny, nz;
	if(dim==1 && (i<1 || i>nx-2))  return 0;
	if(dim==2 && ((i<1 || i>nx-2) || (j<1 || j>ny-2))) return 0;
	if(dim==3 && ((i<1 || i>nx-2) || (j<1 || j>ny-2) || (k<1 || k>nz-2))) return 0;
	return 1;
}

float calculatecurvature(float phi, float lap)
{
	//extern float dx, c_length;
	extern float eps, omega;
	extern float p_threshold;
	if(phi<p_threshold) return 0;
	float width = 2*eps/(sqrt(omega/2));
	float curvilinear_der = (4/width)*phi*(1.0-phi);
	//if(curvilinear_der < threshold) return 0;
	float curvilinear_lap = (32/(width*width))*(phi)*(1.0-phi)*(0.5-phi);
	//printf("lap : %le, curvil_lap: %le, curvil_der : %le, Curvature : %lf\n", lap, curvilinear_lap, curvilinear_der, (lap - curvilinear_lap)/(curvilinear_der));
	return (lap - curvilinear_lap)/(curvilinear_der);
}

float calculatelaplacian(int c, int l, int r, int t, int b, int u, int d, int dim, int grainnum)
{
	extern int nx, ny, nz, grid2;
	extern float *grid;
	extern int *flaggrid;
	extern float dx;
	extern int trackid;

	float cp, lp, rp, tp, bp, up, dp;// Corresponding phi values for neighbouring indices
	float lap=0;

	cp = searchphi(c, grainnum);
	lp = searchphi(l, grainnum);
	rp = searchphi(r, grainnum);

	//if(c==trackid) printf("Grain num : %d, l: %d, lp: %lf, c: %d, cp: %lf, r : %d, rp: %lf\n", grainnum, l, lp, c, cp, r, rp);

	lap = (lp+rp-2*cp)/(dx*dx);
	if(dim>1){if(t!=-1) tp = searchphi(t, grainnum); if(b!=-1) bp = searchphi(b, grainnum);lap += (tp+bp-2*cp)/(dx*dx);}
	if(dim>2){if(u!=-1) up = searchphi(u, grainnum); if(d!=-1) dp = searchphi(d, grainnum); lap += (up+dp-2*cp)/(dx*dx);}
	return lap;
}

float calculategradient(int l, int r, int grainnum)
{
	extern float dx;
	
	float lp, rp;// Corresponding phi values for neighbouring indices, gradient is rp-lp
	float g=0;
	lp = searchphi(l, grainnum);
	rp = searchphi(r, grainnum);
	g = (rp-lp)/(dx);
	return g;
}

float searchphi(int index, int grainnum)
{
	// Search the value of phi corresponding to a given grain number in grid
	extern float *grid;
	extern int *flaggrid;
	extern int ngrainsmax, grid3;
	int l, gind;

	for(l=0; l<ngrainsmax; l++)
	{
		gind  = index + l*grid3;
		if(flaggrid[gind]==-1) break;
		if(flaggrid[gind]==grainnum) return grid[gind];
	}
	return 0;
}

void updatefornextstep()
{
	// Form grid and  flaggrid from newgrid and newflaggrid respectively

	extern float *grid, *newgrid;
	extern int *flaggrid, *newflaggrid;
	extern int nx, ny, nz, N, ngrainsmax;

	extern int trackid;
	extern int *activegrain;
	extern int NUM_THREADS;

	int i, j, k, l, index, gind, count, grainnum;
	float phival;
	float phi[ngrainsmax];
	int flag[ngrainsmax]; //dummy to hold the new values for a grid

	//activegrain = malloc(ngrainsmax * sizeof(int));
	//for(l=0; l<ngrainsmax;l++) activegrain[l]=-1;
	//printf("UPDATING!!\n");
	//#pragma omp parallel for collapse(3) num_threads(NUM_THREADS) schedule(auto) private(count, index, gind,  l, activegrain,  grainnum, phival, phi, flag) shared(flaggrid, grid, newgrid, newflaggrid)
	// For pragma to work have to fix active grain, and move it to form local variable rathre than global variable
	
	for(i=0; i<nx; i++)
	{
	  for(j=0; j<ny; j++)
	  {
	  	for(k=0; k<nz; k++)
	  	{ 
	  	  //printf("%d %d %d\n", i, j, k);
	  	  if(bc==0 || bc==2)
	  	  {
	  	  	if(dim==1) {if(i==0 || i==(nx-1)) continue;}
	  	  	if(dim>1) {if(j==0 || j==(ny-1)) continue;}
	  	  	if(dim>2) {if(k==0 || k==(nz-1)) continue;}
	  	  }
	  	  //First create a list of active grains from all the neighbouring index by using a threshold function	
	  	  count=0;
	  	  index = k*grid2 + j*nx + i;
	  	  
	  	  
	  	  //printf("i: %d, j: %d, k:%d\n", i, j, k);
	  	  
	  	  //if(i==trackid && j==0 && k==0) {printf("Before Assembling grains\n"); for(l=0; l<ngrainsmax; l++) printf("index : %d, l: %d, flaggrain: %d\n", index, l, newflaggrid[l*grid3 + index]);}
	  	  //activegrain = assembleactivegrains(i, j, k, activegrain);//only include non-zero phi(greater than a threshold) from surrounding grains
	  	  assembleactivegrains(i, j, k);
	  	  //if(i==trackid && j==0 && k==0) {printf("Assembling grains\n"); for(l=0; l<ngrainsmax; l++) printf("index : %d, l: %d, activegrain: %d\n", index, l, activegrain[l]);}
	  	  for(l=0; l<ngrainsmax; l++)
	  	  {
	  	  	if(activegrain[l]==-1) break;
	  	  	grainnum = activegrain[l];
	  	  	phival = checkgradient(index, grainnum); //returns the value of phi for the gridpoint if gradient exists; use newgrid here//returns the value if grain exists or returns 0
	  	    phi[count]=phival; flag[count]=grainnum;count+=1;
	  	  }
	  	  
	  	  for(l=count; l<ngrainsmax; l++) {phi[l]=0.0; flag[l]=-1;}
	  	  for(l=0; l<ngrainsmax; l++) {grid[index + l*grid3] = phi[l]; flaggrid[index + l*grid3]=flag[l];}
	  	  
	  	  //if(i==trackid && j==0 && k==0) {printf("Updated!!\n");for(l=0; l<ngrainsmax; l++) printf("Index : %d, l: %d Activegrain: %d, Phi: %lf \n", index, l, flaggrid[l*grid3 + index], grid[l*grid3 + index]);}
	  	  //for(l=0; l<ngrainsmax; l++) {grid[index + l*grid3] = newgrid[index + l*grid3]; flaggrid[index + l*grid3]=newflaggrid[index + l*grid3];} 
	  	}
	  }
	}
	//for(l=0; l<ngrainsmax; l++) printf("Before : l: %d, Flaggrid: %d phi : %lf\n", l, flaggrid[l*grid3 + pos2ind(trackid-4, 0, 0)], grid[l*grid3 + pos2ind(trackid-4, 0, 0)]);
	
	if(bc==0)
	{
	 if(dim==1){for(l=0; l<ngrainsmax; l++) {grid[l*grid3 + pos2ind(0, 0, 0)]=grid[l*grid3 + pos2ind(1, 0, 0)]; grid[pos2ind((nx-1), 0, 0) + l*grid3] = grid[pos2ind((nx-2), 0, 0) + l*grid3]; flaggrid[pos2ind(0, 0, 0) + l*grid3]=flaggrid[pos2ind(1, 0, 0) + l*grid3]; flaggrid[pos2ind(nx-1, 0, 0) + l*grid3] = flaggrid[pos2ind((nx-2), 0, 0) + l*grid3]; }}
	 
	 if(dim==2)
	 {
	 
	   //printf("in BC\n");
	   //j==0
	   for(i=0; i<nx; i++){j=0;k=0;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i, j+1, k)];flaggrid[gind] = flaggrid[l*grid3 + pos2ind(i, j+1, k)];}}
	   //j==ny-1
	   for(i=0; i<nx; i++){j=ny-1;k=0;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i, j-1, k)];flaggrid[gind] = flaggrid[l*grid3 + pos2ind(i, j-1, k)];}}
	   //i==0
	   for(j=0; j<ny; j++){i=0;k=0;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i+1, j, k)];flaggrid[gind] = flaggrid[l*grid3 + pos2ind(i+1, j, k)];}}
	   //i==nx-1
	   for(j=0; j<ny; j++){i=nx-1;k=0;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i-1, j, k)];flaggrid[gind] = flaggrid[l*grid3 + pos2ind(i-1, j, k)];}}
	   
	 }

	 if(dim==3)
	 {
	   //k=0
	   for(i=0; i<nx; i++){for(j=0; j<ny; j++){k=0;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i, j, k+1)];flaggrid[gind] = grid[l*grid3 + pos2ind(i, j, k+1)];}}}
	   //k=nz-1	
	   for(i=0; i<nx; i++){for(j=0; j<ny; j++){k=nz-1;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i, j, k-1)];flaggrid[gind] = grid[l*grid3 + pos2ind(i, j, k-1)];}}}
	   //j=0
	   for(i=0; i<nx; i++){for(k=0; k<nz; k++){j=0;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i, j+1, k)];flaggrid[gind] = grid[l*grid3 + pos2ind(i, j+1, k)];}}}
	   //j=ny-1
	   for(i=0; i<nx; i++){for(k=0; k<nz; k++){j=ny-1;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i, j-1, k)];flaggrid[gind] = grid[l*grid3 + pos2ind(i, j-1, k)];}}}
	   //i=0
	   for(j=0; j<ny; j++){for(k=0; k<nz; k++){i=0;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i+1, j, k)];flaggrid[gind] = grid[l*grid3 + pos2ind(i+1, j, k)];}}}
	   //i=nx-1
	   for(j=0; j<ny; j++){for(k=0; k<nz; k++){i=nx-1;for(l=0; l<ngrainsmax; l++){gind = l*grid3 + pos2ind(i, j, k);grid[gind] = grid[l*grid3 + pos2ind(i-1, j, k)];flaggrid[gind] = grid[l*grid3 + pos2ind(i-1, j, k)];}}}
	 
	  }
	}
	if(bc==2)
	{
	  if(dim==1){for(l=0; l<ngrainsmax; l++) {grid[l*grid3 + pos2ind(0, 0, 0)]=grid[l*grid3 + pos2ind(1, 0, 0)]; grid[pos2ind((nx-1), 0, 0) + l*grid3] = grid[pos2ind((nx-2), 0, 0) + l*grid3]; flaggrid[pos2ind(0, 0, 0) + l*grid3]=flaggrid[pos2ind(1, 0, 0) + l*grid3]; flaggrid[pos2ind(nx-1, 0, 0) + l*grid3] = flaggrid[pos2ind((nx-2), 0, 0) + l*grid3]; }}
	  if(dim==2)
	  {
	  	float p1, p2;
	  	//j==0
	  	for(i=0; i<nx; i++){j=0;k=0;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i, j+1, k)]; p1=searchphi(pos2ind(i, j+1, k), grainnum); p2=searchphi(pos2ind(i, j+2, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}
	     //j==ny-1
	  	for(i=0; i<nx; i++){j=ny-1;k=0;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i, j-1, k)]; p1=searchphi(pos2ind(i, j-1, k), grainnum); p2=searchphi(pos2ind(i, j-2, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}	   
	     //i==0
	  	for(j=0; j<ny; j++){i=0;k=0;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i+1, j, k)]; p1=searchphi(pos2ind(i+1, j, k), grainnum); p2=searchphi(pos2ind(i+2, j, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}
	     //i==nx-1
	  	for(j=0; j<ny; j++){i=nx-1;k=0;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i-1, j, k)]; p1=searchphi(pos2ind(i-1, j, k), grainnum); p2=searchphi(pos2ind(i-2, j, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}
	  }
	  if(dim==3)
	  {
	  	int p1, p2;
	  	//k=0
	    for(i=0; i<nx; i++){for(j=0; j<ny; j++){k=0;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i, j, k+1)]; p1=searchphi(pos2ind(i, j, k+1), grainnum); p2=searchphi(pos2ind(i, j, k+2), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}}
	    //k=nz-1
	    for(i=0; i<nx; i++){for(j=0; j<ny; j++){k=nz-1;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i, j, k-1)]; p1=searchphi(pos2ind(i, j, k-1), grainnum); p2=searchphi(pos2ind(i, j, k-2), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}}
	   	//j=0
	    for(i=0; i<nx; i++){for(k=0; k<nz; k++){j=0;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i, j+1, k)]; p1=searchphi(pos2ind(i, j+1, k), grainnum); p2=searchphi(pos2ind(i, j+2, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}}
	    //j=ny-1
	    for(i=0; i<nx; i++){for(k=0; k<nz; k++){j=ny-1;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i, j-1, k)]; p1=searchphi(pos2ind(i, j-1, k), grainnum); p2=searchphi(pos2ind(i, j-2, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}}
	    //i=0
	    for(j=0; j<ny; j++){for(k=0; k<nz; k++){i=0;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i+1, j, k)]; p1=searchphi(pos2ind(i+1, j, k), grainnum); p2=searchphi(pos2ind(i+2, j, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}}
	    //i=nx-1
	    for(j=0; j<ny; j++){for(k=0; k<nz; k++){i=nx-1;for(l=0; l<ngrainsmax; l++){index = k*grid2+j*nx+i; gind = l*grid3 + index; grainnum = flaggrid[l*grid3 + pos2ind(i-1, j, k)]; p1=searchphi(pos2ind(i-1, j, k), grainnum); p2=searchphi(pos2ind(i-2, j, k), grainnum); grid[gind] = 2*p1-p2; flaggrid[gind] = grainnum;}}}
	  }
	}
	//printf("done BC\n");
	return;
}

void assembleactivegrains(int i, int j, int k)
{
	extern int ngrainsmax, nx, ny, nz, grid2, grid3;
	extern float *newgrid;
	extern int *newflaggrid;
	extern float threshold, p_threshold;
	extern int *activegrain;

	int index = k*grid2 + j*nx + i;
	int count=0, newindex, b=0, l=0, a=0, gind, id=0;

	//if(i==trackid && j==0 && k==0) {printf("Assembling\n"); for(l=0; l<ngrainsmax; l++) printf("index : %d, l: %d, flaggrid: %d phi: %lf\n", index, l, flaggrid[l*grid3 + pos2ind(i, j, k)], newgrid[l*grid3 + pos2ind(i, j, k)]);}

	for(l=0; l<ngrainsmax; l++) activegrain[l]=-1;// Initialize the active list as -1

	for(l=0; l<ngrainsmax; l++) 
	{
	  id=1;
	  gind = l*grid3 + index; 
	  if(newflaggrid[gind]==-1) break;
	  if(newgrid[gind]>p_threshold)
	  {
	  	for(b=0; b<count+1; b++) 
	  	{
	  	  if(activegrain[b]==newflaggrid[gind]) {id=0; break;}
	  	}
	  	if(id && (count<ngrainsmax-1)){activegrain[count] = newflaggrid[gind];count+=1;}
	  }
	  //if(i==trackid && j==0 && k==0) printf(" first During assembling: index : %d, l: %d, count : %d, flaggrid: %d\n", pos2ind(i, j, k), l, count, activegrain[count-1]);
	}
	// Loop for all neighbouring elements in x
	//if(i==trackid-4 && j==0 && k==0){for(l=0; l<count; l++) printf("index : %d, count:%d , activegrain: %d\n", index, l, activegrain[l]);}
	for(a=-1;a<2; a++)
	{
	  if(a==0) continue;
 	  newindex = pos2ind(i+a, j, k); 
	  for(l=0; l<ngrainsmax; l++)
	  {
	  	gind = newindex + l*grid3;
		if(newflaggrid[gind]==-1)break;
		//if(i==trackid && j==0 && k==0){printf("index : %d, l:%d , gind : %d, phi: %lf flaggrid: %d\n", index, l, gind, newgrid[gind], newflaggrid[gind]);}
	  	
		// First check which grains have more than a threshold value, then check if they are already included in the list
		if(newgrid[gind]>p_threshold) 
		{
		  id = 1;
		  for(b=0; b<count+1; b++) 
		  {	
		  	if(newflaggrid[gind]==activegrain[b])
		  	{
		  	  id=0; 
		  	  break;
		  	}
		  }
		  if(id && (count<ngrainsmax-1)) {activegrain[count] = newflaggrid[gind];count+=1;}
		
		}
		//if(i==trackid && j==0 && k==0) printf(" During assembling: a : %d,  index : %d, l: %d, count : %d, flaggrid: %d\n", a, pos2ind(i, j, k), l, count, activegrain[count-1]);
	   }
	}
	//if(i==trackid && j==0 && k==0){for(l=0; l<count; l++) printf("index : %d, count:%d , activegrain: %d\n", index, l, activegrain[l]);}

	if(dim>1)
	{
	  for(a=-1;a<2; a++)
	  { 
	  	if(a==0) continue; 
		newindex = pos2ind(i, j+a, k);
	  	for(l=0; l<ngrainsmax; l++)
		{
		  gind  = newindex + l*grid3;
		  if(newflaggrid[gind]==-1)break;
		  if(newgrid[gind]>p_threshold) 
		  {
		  	id=1;
		  	for(b=0; b<count+1; b++) { if(count==0) continue; if(activegrain[b]==newflaggrid[gind]) {id=0; break;}}
		  	if(id && (count<ngrainsmax-1)) {activegrain[count] = newflaggrid[gind];count+=1;}
		  }
		   
		}
	  }
	}
	if(dim>2)
	{
	  for(a=-1;a<2; a++)
	  { 
	  	if(a==0) continue;
		newindex = pos2ind(i, j, k+a);
		for(l=0; l<ngrainsmax; l++)
		{
		  gind  = newindex + l*grid3;
		  if(newflaggrid[gind]==-1)break;
		  if(newgrid[gind]>p_threshold)
		  {
		  	id=1;
		  	for(b=0; b<count+1; b++) {if(count==0) continue;id=1;if(activegrain[b]==newflaggrid[gind]) {id=0; break;}}
		  	if(id && (count<ngrainsmax-1)) {activegrain[count] = newflaggrid[gind];count+=1;}
		  }
		   
		}
	  }
	}
	//return activegrain;
	return;
}

float checkgradient(int index,int grainnum)
{
	extern float *newgrid;
	extern int *newflaggrid;
	extern int grid3;
	int l, gind;
	for(l=0; l<ngrainsmax; l++)
	{
		gind = index + l*grid3;
		if(newflaggrid[gind]==-1) break;
		if(newflaggrid[gind]==grainnum) return newgrid[gind];
	}
	return 0.0; 



}

void gradientvector(int c, int l, int r, int t, int b, int u, int d, int dim, int grainnum, float *gradphi)
{
	//static float gv[3];
	for(int i=0; i<3; i++) gradphi[i]=0.0;
	gradphi[0] = calculategradient(l, r, grainnum)/2;
	if(dim==2) gradphi[1] = calculategradient(b, t, grainnum)/2;
	if(dim==3) gradphi[2] = calculategradient(d, u, grainnum)/2;
	//return gradphi;

}

float findvelocity(float *gradvector, float dphidt, int dim)
{
	float mag;
	mag = gradvector[0]*gradvector[0];
	if(dim==2) mag += gradvector[1]*gradvector[1];
	if(dim==2) mag += gradvector[2]*gradvector[2];
	
	return -(dphidt/sqrt(mag));
}

float dG_drag(float E0, float c0, float v, float T)
{
	extern float Vm, c_energy;
	return ((adrag(E0, T)/Vm)*v*c0/(1.0 + bdrag(E0, T)*v*v))*(1/c_energy);
}

float adrag(float E0, float T)
{
	extern float dalpha, R, c_diffusivity, sd_delta;
	return ((sd_delta)/(dalpha*E0*c_diffusivity))*(R*T)*(R*T)*(sinh(E0/(R*T)) - E0/(R*T));
}
float bdrag(float E0, float T)
{
	extern float R, dalpha, c_diffusivity, sd_delta;
	return ((sd_delta/2)*adrag(E0, T)*(R*T)/(2*dalpha*E0*E0*c_diffusivity));
}

float dG_drag2(float a, float b, float v, float T)
{
	extern float Vm, c_energy;
	//printf("a: %lf, b : %lf\n", a, b);
	return (a*v/(1.0 + b*v*v))*(1/c_energy);
}
float dG_dragfit(float a, float b, float v, float T)
{
	extern float Vm, c_energy;
	extern float sd_delta, dalpha, c_diffusivity;
	float vnd = (v*sd_delta)/(dalpha*c_diffusivity);
	return (a*vnd/(1.0 + b*vnd*vnd));
}
float iw()
{
	extern float eps, omega, c_length;
	return 2*eps/sqrt(omega/2)*c_length;
}

void printinfo_gridpoint(int i, int j, int k)
{
	int i_print, j_print, k_print;
	i_print = 201; j_print=1; k_print=1;
	if(i!=i_print-1 || j!=j_print-1 || k!=k_print-1) return;
	extern int ngrainsmax, grid2, grid3;
	int l=0;
	int index = k*grid2 + j*nx + i;
	int gind;  
	printf("Position : i : %d, j : %d, k : %d ThreadID : %d\n",i,j,k, omp_get_thread_num());
	for(l=0;l<ngrainsmax;l++)
	{
	  gind = index + l*grid3;
	  if(flaggrid[gind]==-1) break;
	  printf("Phi Value : %lf, Grain number : %d\n", grid[gind], flaggrid[gind]);	
	}
	return;
	
}

float anisotropic_mobility(float mob, float misorientation, int flag)
{
	extern float mu_anisotropy, sigma_anisotropy;
	float value=0;
	if(flag==0) return mob; //Isotropic case
	value = mob*exp(-((misorientation - mu_anisotropy)*(misorientation - mu_anisotropy))/(2*sigma_anisotropy*sigma_anisotropy));
	//if(fabs(misorientation-mu_anisotropy)>1e-2) printf("Mmax : %e, M : %e, misorientation : %lf, mean : %lf, std : %lf\n", mob, value, misorientation, mu_anisotropy, sigma_anisotropy);
	return value;
}



float anisotropic_segregation(float Emax, float misorientation, int flag)
{
	extern float mu_anisotropy, sigma_e0anisotropy;
	float value=0;
	if(flag==0) return Emax; //Isotropic case
	value = Emax*exp(-((misorientation - mu_anisotropy)*(misorientation - mu_anisotropy))/(2*sigma_e0anisotropy*sigma_e0anisotropy));
	//if(fabs(misorientation-mu_anisotropy)>1e-2) printf("Emax : %e, E : %e, misorientation : %lf, mean : %lf, std : %lf\n", Emax, value, misorientation, mu_anisotropy, sigma_e0anisotropy);
	
	//printf("MO : %e, M : %e, misorientation : %lf, mean : %lf, std : %lf\n", mob, value, misorientation, mu_anisotropy, sigma_anisotropy);
	return value;
}

float evaluatemisorientation(int grain1, int grain2)
{
	extern float **orientation, **fullmisorientation;
	extern int flag_anisotropy;
	float value;

	if(flag_anisotropy==2) value = fabs(orientation[grain1][0]-orientation[grain2][0]);

	else value = fullmisorientation[grain1][grain2];

	return value;
}

