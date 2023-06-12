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



void pfsolver(int t){
  extern int nx, ny, nx, dim, grid2, grid3, ngrainsmax;
  extern float *grid, *newgrid, *curvature;
  extern float *dphidtold;
  extern int *flaggrid, *newflaggrid;
  extern float **orientation, **rotationmatrix;
  extern int bc, N, M;
  extern float dt, mob, eps, omega;
  extern float gamma_max, gamma_min;
  extern float eps_max, eps_min;
  extern float c_length, c_time;
  extern float threshold, p_threshold;
  extern float E0, T;
  extern float beta, beta_scaled;
  extern float  mu_anisotropy, sigma_anisotropy, sigma_e0anisotropy;
  extern int flag_anisotropy, flag_e0anisotropy, flag_gbeanisotropy;
  extern float gbe_ratio;
  extern int trackid, flag_curvature;
  extern int t_output;
  extern int NUM_THREADS;
  extern float iter_E0;
  extern float expected_velocity, pvelocity;
  extern int vel_count;
  extern float m0, qm;
  extern int flag_tt;
  extern float sd_delta;

  int solution;
  int left, right, top, down, up, bottom;
  int iter, counter;
  int iter_first=200, iter_later=30, iter_sim;
  float sumphi2=0, sumphi;
  float misorientation = 45;
  float beta_eff=0;
	
  int i=0, j=0, k=0, l=0, id, m;
  float grad;
  int index, count, gind;
  float velocity, velocity_prev, dphidt, dv, vel;
  float gradphi[3];
  float constantrate=0.5, cr;
  float cr_first=0.5, cr_later=0.1;
  float Eeff_iter=E0, Meff_iter=mob;
  float eps_iter = eps, omega_iter = omega;
  float gamma_iter;
  float sum_gamma=0;
  float delta_theta;
  float dphidtpast, error, mag_gradphi, dummy;
  int l2;
  float totaldragcomponent, mis2;
  for(int i=0; i<3; i++) gradphi[i]=0;
  if(flag_curvature && (t%t_output)==0) {for(i=0; i<N; i++) curvature[i] = 0.0;}
  // Preparing the stencils
  if(t==0) updatefornextstep();

  // initializing constants mobility and dt
  evaluate_constants();
  T = calc_temperature();
  pvelocity = 0;
  vel_count = 0;
	
/*Loop over all the gridpoints*/
  #pragma omp parallel for collapse(3) num_threads(NUM_THREADS) schedule(auto) private(left, right, top, bottom, up, down, gradphi, count, sumphi2, sumphi, beta_eff, index, gind, id, l, l2,grad, misorientation, iter, dphidt, velocity, velocity_prev, dv, vel, cr, m, error, dphidtpast, mag_gradphi, counter, dummy, iter_sim, solution, Meff_iter, Eeff_iter, eps_iter, omega_iter, gamma_iter, delta_theta, sum_gamma, totaldragcomponent) shared(t, flaggrid, newflaggrid, grid, newgrid, curvature, orientation, rotationmatrix, dphidtold, pvelocity, vel_count, beta_scaled)
  for(i=0; i<nx; i++){
    for(j=0; j<ny; j++){
      for(k=0; k<nz; k++){
        // initializers at each grid point
        Meff_iter = mob;Eeff_iter = E0;eps_iter = eps;omega_iter = omega;gamma_iter = gamma_max;
        count=0; sumphi2=0; sum_gamma = 0;sumphi=0;beta_eff = 0;
        
        //Only solve if more than 2 active grains
	index = k*grid2 + j*nx + i; //index in the 3D grid if present, grid2 = nx*ny, grid3 =nx*ny*nz - For 2D grid2 = grid3, k=0
	for(l=0; l<ngrainsmax;l++){
	  if(flaggrid[index + l*grid3]==-1) break; 
	  count+=1;
	  gind = index + l*grid3; sumphi2 += grid[gind]*grid[gind];
	  sumphi += grid[gind];
	  if(dim==1) beta_eff += beta_table[flaggrid[gind]]*grid[gind]*grid[gind]; 
	}
	if(dim==1) beta_eff /= sumphi2;
	
	if(t>2) {if(count<2) continue;}
	left = -1; right = -1;  top=-1;  bottom=-1; up=-1; down=-1; //set everything to -1
	// Boundary conditions
	// If no flux condition, only solve for internal nodes - also considers 2d and 3d system
	if(bc==0 || bc==2){id = boundarycondition(i, j, k, dim);if(!id) continue;}
          
        // Assign the neighbouring index considering boundary conditions
	left = pos2ind(i-1, j, k); right = pos2ind(i+1, j, k);
	if(dim>1) {top  = pos2ind(i, j+1, k);bottom = pos2ind(i, j-1, k);}
	if(dim>2) {up = pos2ind(i, j, k+1);down =  pos2ind(i, j, k-1);}
	  	
	/***Misorientation and effective parameters calculations *****/
	// M = sum_i sum_j!=i M_ij eta_i^2 eta_j^2/ sum_i sum_j!=i eta_i^2 eta_j^2
	/***Only calculate when anisotropic flag is activated ****/
	if((flag_anisotropy||flag_e0anisotropy||flag_gbeanisotropy) && t>2){
	  //Initializers
	  Meff_iter = 0; Eeff_iter=0; eps_iter=0; gamma_iter=0;
	  float fac=0, sumsq2=0;
	  float maxphi=0.0, maxphi1=0.0;
	  int maxindex=0, maxindex1=1;
	  int grain1=maxindex, grain2=maxindex1;
	  for(l=0; l<ngrainsmax;l++){
	    gind = index+l*grid3;
	    if(grid[gind]<p_threshold) continue;
	    for(l2=0; l2<ngrainsmax; l2++){
	      if((grid[index+l2*grid3]<p_threshold)||(l2==l)) continue;
	      fac = grid[index+l2*grid3]*grid[index+l2*grid3]*grid[gind]*grid[gind]; sumsq2+=fac;
	      delta_theta = evaluatemisorientation(flaggrid[gind], flaggrid[index+l2*grid3]);
	      Meff_iter += anisotropic_mobility(mob, delta_theta, flag_anisotropy)*fac; 
	      Eeff_iter += anisotropic_segregation(E0, delta_theta, flag_e0anisotropy)*fac;
	      eps_iter += anisotropic_eps(eps_max, delta_theta, flag_gbeanisotropy)*fac;        					
	    }
	    if(grid[gind]>maxphi) {maxindex=gind; maxphi=grid[maxindex]; grain1 = flaggrid[maxindex];}
	  }
	  for(l=0; l<ngrainsmax;l++){
	    gind = index+l*grid3;
	    if(grid[gind]>maxphi1 && (gind!=maxindex)) {maxindex1=gind; maxphi1=grid[maxindex1]; grain2=flaggrid[maxindex1];}
	  }
	  if(sumsq2<1e-6){
	    misorientation = evaluatemisorientation(grain1, grain2);
	    Meff_iter = anisotropic_mobility(mob, misorientation, flag_anisotropy); Eeff_iter = anisotropic_segregation(E0, misorientation, flag_e0anisotropy);
	    eps_iter = anisotropic_eps(eps, misorientation, flag_gbeanisotropy);
	  }
	  else{
	    Meff_iter = Meff_iter/sumsq2; Eeff_iter = Eeff_iter/sumsq2;eps_iter = eps_iter/sumsq2;
	  }
	  misorientation = evaluatemisorientation(grain1, grain2); // misorientation between grains with maximum indices
	}
	else misorientation = mu_anisotropy;

	
	//Phase field equation for each order parameter eta_i
	//d(eta_i)/d(t) = L(eps^2(grad eta_i)^2 - m(eta_i^3 - eta_i + 2*eta_i*sum_j eta_j^2 gamma_ij) + 3*eta_i*sum_j eta_j Gsd_ij)
	for(l=0; l<ngrainsmax; l++){
	  gind = index + l*grid3;
	  if(flaggrid[gind]==-1) break;
	  // While forming grad for the next timestep take care to only keep all the grains with non-zero gradients only elimintate phi with 0 at the end
	  grad = calculatelaplacian(index, left, right, top, bottom, up, down, dim, flaggrid[gind]); //Grain number at the end  
	  if(flag_curvature && (t%(t_output)==0) && t>2) {/*printf("i : %d, j:%d k:%d phi: %lf\n", i, j, k, grid[gind]);*/curvature[gind]=calculatecurvature(grid[gind], grad);}
	  
	  //Calculation of sum_j eta_j^2 gamma_ij
	  sum_gamma = 0;
	  for(l2=0; l2<ngrainsmax; l2++){
	    if((grid[index+l2*grid3]<p_threshold)||(l2==l)) continue; 
	    if(flaggrid[gind]==-1) break;
	    if((flag_anisotropy||flag_e0anisotropy||flag_gbeanisotropy) && t>2) delta_theta = evaluatemisorientation(flaggrid[gind], flaggrid[index+l2*grid3]);
	    else delta_theta = mu_anisotropy; 
	    sum_gamma+=anisotropic_gamma(gamma_max, delta_theta, flag_gbeanisotropy)*grid[index + l2*grid3]*grid[index + l2*grid3];
	  }
	  if(sum_gamma<1e-6) sum_gamma=anisotropic_gamma(gamma_max, misorientation, flag_gbeanisotropy)*(sumphi2 - grid[gind]*grid[gind]);
	  //printf("%d, %lf, %lf, %lf %lf\n",flaggrid[gind], anisotropic_gamma(gamma_max, misorientation, flag_gbeanisotropy), sum_gamma, grid[gind]*grid[gind], sumphi2);
	  //for(l2=0;l2<ngrainsmax;l2++){printf("%d, %lf    ", flaggrid[index + l2*grid3], grid[index +  l2*grid3] );}
	  //printf("\n");
	  // If dim =1 and artificial driving force is active or else normal grain growth
	  if(dim==1){
	    newgrid[gind] = grid[gind] + dt*(Meff_iter)*(eps_iter*eps_iter*grad - omega_iter*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 2*sum_gamma*grid[gind]) + (2/0.95)*grid[gind]*(beta_table[flaggrid[gind]] - beta_eff)/sumphi2);
	    } 
	  else newgrid[gind] = grid[gind] + dt*(Meff_iter)*(eps_iter*eps_iter*grad - omega_iter*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 2*sum_gamma*grid[gind]));
	 
	 
	  // if solute drag is active
	  cr = constantrate;
	  dphidt = (newgrid[gind]-grid[gind])/dt; // from normal grain growth solution
	  
	  //Can only solve equation at the interface since velocity is computed only at the interface
	  if(flag_sd && (grid[gind]>0.01) && (grid[gind]<(1-0.01)) && (t>100)){
	    //Initializers
	    solution=0;counter=0;error = 1;iter=0;
	    for(m=0; m<3; m++) gradphi[m]=0;
	    //fill the gradient vectors
	    gradientvector(index, left, right, top, bottom, up, down, dim, flaggrid[gind], gradphi);
	    mag_gradphi = sqrt(gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1]);
	    dummy = dphidtold[gind];
	    if(t<500) {iter_sim = iter_first; cr=cr_first;}
	    else {iter_sim = iter_later; cr=cr_later;}
	    
	    //Iteration loop
	    while((error > 0.1) && (iter<iter_sim)){
	      //printf("Velocity: %lf\n",findvelocity(gradphi, dphidtpast, dim));
	      
	      dphidtpast = dphidt;
	      if(counter==0) dphidtpast = dummy;
	      velocity = findvelocity(gradphi, dphidtpast, dim);
	      totaldragcomponent=0; 
	      if(fabs(mag_gradphi)<1e-5) {velocity=0; break;}
	      int dim_factor=1;
	      if(dim>1) dim_factor=0;
	      // If reading from a lookup table
	      if(flag_sd==2){
	        totaldragcomponent = dG_drag_lookup(velocity*c_length/c_time)*(sumphi-grid[gind]);
	      }
	      else{
	      /*
	        for(l2=0; l2<ngrainsmax; l2++){
	          if((grid[index+l2*grid3]<p_threshold)||(l2==l)) continue; 
	          if(flaggrid[index + l2*grid3]==-1) break;
	          if((flag_anisotropy||flag_e0anisotropy||flag_gbeanisotropy) && t>2) misorientation = evaluatemisorientation(flaggrid[gind], flaggrid[index+l2*grid3]);
	          else misorientation = mu_anisotropy; 
	          totaldragcomponent+=dG_drag(anisotropic_segregation(E0, misorientation, flag_e0anisotropy), c_init, velocity*c_length/c_time, T, delta_fn(misorientation, flag_e0anisotropy))*grid[index + l2*grid3];
	        }
	      */ 
	        totaldragcomponent = dG_drag(Eeff_iter, c_init, velocity*c_length/c_time, T, delta_fn(misorientation, flag_e0anisotropy))*(sumphi - grid[gind]);
	      }
	      //dphidt = Meff_iter*(eps_iter*eps_iter*grad - omega_iter*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 2*sum_gamma*grid[gind]) + 3*totaldragcomponent*grid[gind] + dim_factor*(2/0.95)*grid[gind]*(beta_table[flaggrid[gind]] - beta_eff)/sumphi2);	  	  	
	      dphidt = Meff_iter*(eps_iter*eps_iter*grad - omega_iter*(grid[gind]*grid[gind]*grid[gind] - grid[gind] + 2*sum_gamma*grid[gind]) + 3*totaldragcomponent*grid[gind]);	  	  	
	      // check the error for loop
	      error = fabs((dphidtpast-dphidt)/dphidtpast);
	      dphidt = dphidtpast + cr*(dphidt - dphidtpast);	
	      counter+=1;
	      if(dphidt*dphidtpast<0) {counter=0; cr/=2;}
	      {if(iter%30==0 && iter<iter_sim) {velocity=0;dummy/=2;counter=0;}}
	      if(fabs(dphidt)<1e-15) {dphidt=0;break;}
	      iter += 1;
	  	  	  	
	      //if(iter==51) printf("THIS IS NOT WORKING YAAAR, time : %d, counter: %d error : %le, i: %d, j:%d, l:%d t:%d cr :%le dphidt: %le dphidtpast: %le\n", t, counter,error, i, j, l, t, cr, dphidt, dphidtpast);
	    }
	    
	    if(velocity>1e-15){
	      vel_count +=1;
	      pvelocity+=velocity*c_length/c_time;
	    }
	    newgrid[gind] = grid[gind] + dt*dphidt;
	    if(newgrid[gind]<0.0) newgrid[gind]=0.0;
	    else if(newgrid[gind]>1.0) newgrid[gind]=1.0;
	  }
	  dphidtold[gind]=dphidt;
	  
	  if(newgrid[gind]<0.0) newgrid[gind]=0.0;
	  else if(newgrid[gind]>1.0) newgrid[gind]=1.0;
	  newflaggrid[gind] = flaggrid[gind];
	  	  	
	}
      }
    }
  }
  pvelocity /= vel_count;
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
	float curvilinear_lap = (32/(width*width))*(phi)*(1.0-phi)*(0.5-phi);
	//float curvilinear_der = -(8/width)*phi*(1.0-phi);
	//float curvilinear_lap = (32*2/(width*width))*(phi)*(1.0-phi)*(0.5-phi);
	
	return fabs((lap - curvilinear_lap)/(curvilinear_der));
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
	  	  	if(dim>0) {if(i==0 || i==(nx-1)) continue;}
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

float dG_drag(float E0, float c0, float v, float T, float delta)
{
	extern float threshold;
	if(fabs(E0) < threshold) return 0;
	extern float Vm, c_energy;
	return ((adrag(E0, T, delta))*v*c0/(1.0 + bdrag(E0, T, delta)*v*v))*(1/c_energy);
}

float adrag(float E0, float T, float delta)
{
	extern float dalpha, R, c_diffusivity, sd_delta;
	//return ((sd_delta)/(dalpha*E0*Vm*c_diffusivity))*(R*T)*(R*T)*(sinh(E0/(R*T)) - E0/(R*T));
	return ((delta)/(dalpha*E0*Vm*c_diffusivity))*(R*T)*(R*T)*(sinh(E0/(R*T)) - E0/(R*T));
}
float bdrag(float E0, float T, float delta)
{
	extern float R, dalpha, c_diffusivity, sd_delta;
	//return (Vm*(sd_delta/2)*adrag(E0, T)*(R*T)/(2*dalpha*E0*E0*c_diffusivity));
	return (Vm*(delta/2)*adrag(E0, T, delta)*(R*T)/(2*dalpha*E0*E0*c_diffusivity));
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

float dG_drag_lookup(float v)
{
	//input is actual velocity
	extern float Vm, c_energy;
	extern float **sdinfo;
	extern float sdmin, sdmax;
	extern int sdcount;
	extern int flag_sd;
	if(flag_sd!=2) {printf("Set flag SD = 2 for lookup table calculations\n"); exit(0);}
    

	// lookup table should have nondimensional velocity (vb/D) and dG in J/mol
	float vnd = fabs(v*2.88e-10/(dalpha*c_diffusivity));
	if((vnd)<1e-8) return 0;
	if((vnd)>1e+4) return 0;
	
	float sign;
	if(vnd*v<0) sign=-1;
	if(vnd*v>=0) sign=1;
	
	int l_index = (int)((log10(vnd) - sdmin)*(sdcount-1)/(sdmax-sdmin));
	int h_index = l_index+1;
	//printf("vnd : %e, h_index : %d l_index : %d\n", vnd, h_index, l_index);
	float l_dg = sdinfo[l_index][1];
	float h_dg = sdinfo[h_index][1];

	float dg_Jmol = h_dg + ((h_dg-l_dg)/(sdinfo[h_index][0]-sdinfo[l_index][0]))*(vnd - sdinfo[h_index][0]);
	
	//printf("v->actual : %e, v->nd : %e, hi-index : %d, high-dG : %e dg_Jmol : %e\n", v, vnd, h_index, h_dg, dg_Jmol);
	float dg_inJm3 = dg_Jmol/Vm;
	return sign*dg_inJm3/c_energy;
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
	extern float mu_anisotropy, sigma_anisotropy, mob_ratio;
	extern float mu_anisotropy_array[4], mob_ratio_array[4];
	float value=0;
	if(flag==0) return mob; //Isotropic case

	// For the cutoff, here mu_aisotropy refers to the cutoff angle
	if(flag==3)
	{
	  //printf("Mobility ratio : %lf\n", mob_ratio);
	  if(mob_ratio<1)
	  {
	  	//printf("r<1");
	  	if(misorientation<mu_anisotropy) value=mob;
	    else if(misorientation>=mu_anisotropy) value=mob_ratio*mob;
	  }
	  else
	  {
	  	//printf("r>1");
	  	if(misorientation<mu_anisotropy) value=(1.0/mob_ratio)*mob;
	    else if(misorientation>=mu_anisotropy) value=mob;
	  }
	  
	}
	else if(flag==4)
	{
		if((misorientation>=0.0) && (misorientation<mu_anisotropy_array[0])) {value = mob*(1.0/mob_ratio_array[0]);}
		if((misorientation>mu_anisotropy_array[0]) && (misorientation<mu_anisotropy_array[1])) {value = mob*(1.0/mob_ratio_array[1]);}
		if((misorientation>mu_anisotropy_array[1]) && (misorientation<mu_anisotropy_array[2])) {value = mob*(1.0/mob_ratio_array[2]);}
		if((misorientation>mu_anisotropy_array[2]) && (misorientation<mu_anisotropy_array[3])) {value = mob*(1.0/mob_ratio_array[3]);}
	}
	// Remove the top portion to accomodate for gaussian change in mobility
	else value = mob*exp(-((misorientation - mu_anisotropy)*(misorientation - mu_anisotropy))/(2*sigma_anisotropy*sigma_anisotropy));
	//if(fabs(misorientation-mu_anisotropy)>1e-2) printf("Mmax : %e, M : %e, misorientation : %lf, mean : %lf, std : %lf\n", mob, value, misorientation, mu_anisotropy, sigma_anisotropy);
	return value;
}



float anisotropic_segregation(float Emax, float misorientation, int flag)
{
	extern float mu_e0anisotropy, sigma_e0anisotropy;
	extern float mu_e0anisotropy_array[4], E0_array[4];
	
	float value=0;
	if(flag==0) return Emax; //Isotropic case

	if(flag==1)
	{
	  // For the cutoff, here mu_aisotropy refers to the cutoff angle
	  value = Emax*exp(-((misorientation - mu_e0anisotropy)*(misorientation - mu_e0anisotropy))/(2*sigma_e0anisotropy*sigma_e0anisotropy));
	}
	if(flag==2)
	{
	  if(misorientation<mu_e0anisotropy) value = Emax;
	  else if(misorientation>=mu_e0anisotropy) value = (1+dEpct/100)*Emax;	
	}
	if(flag==3)
	{
	  if(misorientation>=mu_e0anisotropy) value = Emax;
	  else if(misorientation<mu_e0anisotropy) value = (1+dEpct/100)*Emax;	
	}
	if(flag==4)
	{
		if((misorientation>=0.0) && (misorientation<mu_e0anisotropy_array[0])) {value = E0_array[0];}
		if((misorientation>mu_e0anisotropy_array[0]) && (misorientation<mu_e0anisotropy_array[1])) {value = E0_array[1];}
		if((misorientation>mu_e0anisotropy_array[1]) && (misorientation<mu_e0anisotropy_array[2])) {value = E0_array[2];}
		if((misorientation>mu_e0anisotropy_array[2]) && (misorientation<mu_e0anisotropy_array[3])) {value = E0_array[3];}
		//printf("Misorientation: %f, E0: %e\n", misorientation, value);
	}

	//if(fabs(misorientation-mu_anisotropy)>1e-2) printf("Emax : %e, E : %e, misorientation : %lf, mean : %lf, std : %lf\n", Emax, value, misorientation, mu_anisotropy, sigma_e0anisotropy);
	
	//printf("MO : %e, M : %e, misorientation : %lf, mean : %lf, std : %lf\n", mob, value, misorientation, mu_anisotropy, sigma_anisotropy);
	return value;
}
float delta_fn(float misorientation, int flag)
{
	extern float sd_delta, sd_delta_array[4];
	extern float mu_e0anisotropy_array[4];
	float value=1e-10;
	//return sd_delta;
	
	if(flag!=4) return sd_delta;
	else
	{
		if((misorientation>=0.00) && (misorientation<mu_e0anisotropy_array[0])) {value = sd_delta_array[0];}
		if((misorientation>=mu_e0anisotropy_array[0]) && (misorientation<mu_e0anisotropy_array[1])) {value = sd_delta_array[1];}
		if((misorientation>=mu_e0anisotropy_array[1]) && (misorientation<mu_e0anisotropy_array[2])) {value =sd_delta_array[2];}
		if((misorientation>=mu_e0anisotropy_array[2]) && (misorientation<mu_e0anisotropy_array[3])) {value =sd_delta_array[3];}
	}
	//if(misorientation<0) printf("Misorientation: %f, delta: %le\n", misorientation, value);
	return value;
	
}

float evaluatemisorientation(int grain1, int grain2)
{
	extern float **orientation, **fullmisorientation;
	extern int flag_anisotropy;
	extern int ngrainsmax, ngrainsor;
	float value;
	//printf("Grain 1: %d, Grain 2: %d, orientation id 1: %d, orientation id 2: %d\n", grain1, grain2, grain1%ngrainsor, grain2%ngrainsor);
	int i=0, j=0;
	if(flag_anisotropy==2) value = fabs(orientation[grain1%ngrainsor][0]-orientation[grain2%ngrainsor][0]);

	else value = fullmisorientation[grain1%ngrainsor][grain2%ngrainsor];

	return value;
}

void evaluate_constants()
{
	extern float k_dt, dx, mob, eps;
	extern float m0, qm, T, d0, qd, R, c_diffusivity;
	extern float dalpha;
	extern float sim_time;
	extern float mob_ratio;
	extern int flag_tt;
	extern float mu_anisotropy;

	float m1;
	T = calc_temperature();

	mob = mob_nd(mob_int(m0, qm, T));
	m1 = mob;
  	dalpha = d0*exp(-qd/(R*T))/c_diffusivity;
    
    //printf("Mobility constant = %le Phase field mobility = %le Interface Mobility : %le\n", c_mobility, mob, mob_int(m0, qm, T));
  	//printf("Diffusion constant = %le\n", dalpha);
  	//printf("time : %e, Temperature = %le\n", sim_time, T);
    if(flag_anisotropy==3)
  	{
  	  float m1 = anisotropic_mobility(mob, mu_anisotropy+1, flag_anisotropy);	
  	  float m2 = anisotropic_mobility(mob, mu_anisotropy-1, flag_anisotropy);	
  		dt = fmin(k_dt*dx*dx/(m2*eps*eps), k_dt*dx*dx/(m1*eps*eps));	
  	  //dt = fmin(k_dt*dx*dx/(m2), k_dt*dx*dx/(m1));	
  	  dt = fmin(dt, k_dt*dx*dx/dalpha);  	  
  	} 
  	else
  	{
  	  dt = fmin(k_dt*dx*dx/(mob*eps*eps), k_dt*dx*dx/dalpha);
  	}
  	//printf("dt = %le\n", dt);
  	
  	return;
}

float calc_temperature()
{
	extern float sim_time;
	extern int flag_tt;
	extern float T;
	extern float **ttinfo;
	extern int ttcount;

	float updatedT;
	if(flag_tt==0) return T;
	
	int index=0;

	if(sim_time<1e-3*threshold) return ttinfo[0][1];
	for(int i=0; i<ttcount; i++)
	{
		if(sim_time<ttinfo[i][2]) {index=i; break;} 
	}
	updatedT = ttinfo[index][1] + ((ttinfo[index][3]-ttinfo[index][1])/(ttinfo[index][2]-ttinfo[index][0]))*(sim_time - ttinfo[index][0]); 
	//printf("updatedT : %le\n", updatedT);
	return updatedT;
}

float effective_velocity(float E0, float T)
{
	extern float m0, qm;
	extern float d0, qd, dx;
	extern float gbenergy;
	extern int flag_geometry;
	extern int grid3;
	extern float sd_delta;
	
	float average_radius;
	float effective_mobility;
	float v;

	if(flag_geometry==2) average_radius= dx*sqrt(evaluatevolumefraction()*grid3/PI);
	else average_radius = dx*evaluateaverageradius();
	effective_mobility = pow((1/mob_int(m0, qm, T) + adrag(E0, T, sd_delta)/Vm), -1);
	return effective_mobility*gbenergy/average_radius;
}

float anisotropic_eps(float e, float misorientation, int flag)
{
	extern float mu_anisotropy, sigma_anisotropy, gbe_ratio;
	extern float eps_min,eps_max;
	float value=e;
	if(flag==1){
	  /*
	  if(misorientation<mu_anisotropy) value=sqrt(1.0/gbe_ratio)*e;
	  else if(misorientation>=mu_anisotropy) value=e;
	  */
	  if(misorientation<mu_anisotropy) value=eps_min;
	  else if(misorientation>=mu_anisotropy) value=eps_max;
	  
	}
	return value;
}

float anisotropic_omega(float w, float misorientation, int flag)
{
	// Make sure the gbe_ratio is > 1
	// otherwise the timestep requirement of FDM maynot be satisfied
	extern float mu_anisotropy, sigma_anisotropy, gbe_ratio;
	float value=w;
	if(flag==1){
	  if(misorientation<mu_anisotropy) value=(1.0/gbe_ratio)*w;
	  else if(misorientation>=mu_anisotropy) value=w;
	}
	return value;
}

float anisotropic_gamma(float gm, float misorientation, int flag){
	extern float mu_anisotropy, sigma_anisotropy;
	extern float gamma_max, gamma_min;
	float value = gm;
	if(flag==1){
	  if(misorientation<mu_anisotropy) value=gamma_min;
	  else value=gamma_max; 	
	}
	return value;
}
