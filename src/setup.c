#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<omp.h>
#include"Global.h"
#include"functions.h"
#include <gsl/gsl_rng.h>


void simulationsetup()
{
	int success;
	extern int flag_geometry;
	extern int readfromrestart;
	extern int numbergrains;
	extern FILE *fs;
	extern float real_totaltime, sim_time, dt;
	extern float c_time, time;

	readparamsfromfile();
	nondimensionalize();

	if(readfromrestart) flag_geometry=0;
	success = initializegeometry(flag_geometry);
	if(success) printf("Geometry Initialized successfully\n\n");
	else
	{
		printf("Error!!! Look in the code in setup.c\n");
		exit(0);
	}

	time = (real_totaltime-sim_time)/c_time;
  timesteps = (int) (time/dt);
  printf("Estimated timesteps : %d\n", timesteps);

  if(readfromrestart && (sim_time>0)) fs = fopen("statistics.txt", "a"); // if restarting from a non-zero time, append else create new file
	else fs = fopen("statistics.txt", "w");

	numbergrains = evaluatenumberofgrains();
  printf("Evaluated number of grains!");
	return ;
}

void readparamsfromfile()
{
  char *inputfilename = "inp/input.txt";
  FILE *fp = fopen(inputfilename, "r");

  extern int dim, nx, ny, nz;
  extern float m0, qm, d0, qd;
  extern float eps, omega, e0, T;
  extern float gbenergy, gbpoints; 
  extern float c_init, vf, dalpha, dbeta;
  extern float R, Tref;
  extern float time;
  extern float beta, vf;
  extern float sd_delta;
  extern float iter_E0;
  extern float mob_ratio, dEpct;
  extern float mob_ratio_array[4];
  extern float Vm;

  extern float real_dx, real_totaltime;
  extern float eps, omega;
  extern float E0;
  extern float E0_array[4]; //initializing incase the special case is activated
  extern float mu_e0anisotropy_array[4];
  extern float mu_anisotropy_array[4];
  extern float sd_delta_array[4];
  extern int bc;
  extern int ngrainsmax, ngrainsinit, ngrainsor; 
  //ngrainsmax is the list of maimum active grains
  //ngrainsinit is the initial number of grains
  //ngrainsor is the limit of orientation generation and should be smaller than 5000 

  extern int flag_geometry, flag_curvature, flag_anisotropy, flag_e0anisotropy;
  extern int flag_sd, flag_tt;
  extern int readfromrestart, writetorestart;
  extern int NUM_THREADS;
  extern int t_output, t_print, t_createcheckpoint;
  extern int graincontrol;
  extern float threshold, p_threshold;
  extern float initaverageradius;
  extern float mu_anisotropy, sigma_anisotropy, mu_e0anisotropy, sigma_e0anisotropy;
  extern int randomseed;
  extern int flag_gbeanisotropy;
  extern float gbe_ratio; //Only supports two types of grain boundaries
  
  NUM_THREADS = 1;
  E0 = 0;
  flag_sd = 0;
  flag_tt = 0;
  threshold = 1e-7; //for all other measures
  p_threshold=1e-7; // Threshold for storing phi grains
  flag_e0anisotropy=0;
  flag_anisotropy=0;
  flag_gbeanisotropy=0;
  mu_anisotropy = 45; //Typically the maxima
  mu_e0anisotropy = 45; //Typically the maxima
  sigma_anisotropy = 120; // High value such that leads to constant value
  sigma_e0anisotropy = 120; // High value such that leads to constant value
  sd_delta = -1e-4;
  dEpct=0;
  gbe_ratio = 1;

  printf("%s : Reading Parameters\n", inputfilename);
  char tmpstr[20], tmpval[20];
  char tempbuffer[100];
  while(!feof(fp))
  {
	  if(fgets(tempbuffer, 100, fp))
	  {
	    sscanf(tempbuffer, "%20s : %20[^;]", tmpstr, tmpval);
	    tmpstr[strcspn(tmpstr, "\r\n")] = 0;
	    tmpval[strcspn(tmpval, "\r\n")] = 0;

	    if(!strcmp(tmpstr, "dim")) dim = atof(tmpval);
	    if(!strcmp(tmpstr, "nx")) nx = atof(tmpval);
	    if(!strcmp(tmpstr, "ny")) ny = atof(tmpval);
	    if(!strcmp(tmpstr, "nz")) nz = atof(tmpval);
	    if(!strcmp(tmpstr, "dx")) real_dx = atof(tmpval);
	  
	    if(!strcmp(tmpstr, "m0")) m0 = atof(tmpval);
	    if(!strcmp(tmpstr, "qm")) qm = atof(tmpval);
	    if(!strcmp(tmpstr, "d0")) d0 = atof(tmpval);
	    if(!strcmp(tmpstr, "qd")) qd = atof(tmpval);
	    
      if(!strcmp(tmpstr, "gbenergy")) gbenergy = atof(tmpval);
      if(!strcmp(tmpstr, "gbpoints")) gbpoints = atof(tmpval);
     
	    //if(!strcmp(tmpstr, "eps")) eps = atof(tmpval);
	    //if(!strcmp(tmpstr, "omega")) omega = atof(tmpval);
	   
	    if(!strcmp(tmpstr, "totaltime")) real_totaltime = atof(tmpval);

	  //if(!strcmp(tmpstr, "tsteps")) timesteps = atoi(tmpval);
	  
	    if(!strcmp(tmpstr, "cinit")) c_init = atof(tmpval);
	  
	  
	    if(!strcmp(tmpstr, "beta")) beta = atof(tmpval); //in J/mol
	    if(!strcmp(tmpstr, "temperature")) T = atof(tmpval);
	    if(!strcmp(tmpstr, "bc")) bc = atoi(tmpval);
	    if(!strcmp(tmpstr, "ngrainsmax")) ngrainsmax = atoi(tmpval); // Maximum active grains
	    if(!strcmp(tmpstr, "vf")) vf = atof(tmpval);
	    if(!strcmp(tmpstr, "flag_geometry")) flag_geometry = atoi(tmpval);
	    if(!strcmp(tmpstr, "flag_curvature")) flag_curvature = atoi(tmpval);
      if(!strcmp(tmpstr, "flag_anisotropy")) flag_anisotropy = atoi(tmpval);
      if(!strcmp(tmpstr, "flag_e0anisotropy")) flag_e0anisotropy = atoi(tmpval);
      if(!strcmp(tmpstr, "flag_gbeanisotropy")) flag_gbeanisotropy = atoi(tmpval);
      if(!strcmp(tmpstr, "mu_anisotropy"))
      {
      	if(flag_anisotropy==4) sscanf(tempbuffer, "%20s : %f %f %f %f", tmpstr, &mu_anisotropy_array[0], &mu_anisotropy_array[1], &mu_anisotropy_array[2], &mu_anisotropy_array[3]);
      	else{mu_anisotropy = atof(tmpval);}
      } 
      
      if(!strcmp(tmpstr, "mu_e0anisotropy"))
      {
      	if(flag_e0anisotropy==4) sscanf(tempbuffer, "%20s : %f %f %f %f", tmpstr, &mu_e0anisotropy_array[0], &mu_e0anisotropy_array[1], &mu_e0anisotropy_array[2], &mu_e0anisotropy_array[3]);
      	else{mu_e0anisotropy = atof(tmpval);}
      }

      if(!strcmp(tmpstr, "E0"))
      {
      	if(flag_e0anisotropy==4) sscanf(tempbuffer, "%20s : %f %f %f %f", tmpstr, &E0_array[0], &E0_array[1], &E0_array[2], &E0_array[3]);
      	else{E0 = atof(tmpval);}
      }
      
      
      if(!strcmp(tmpstr, "flag_sd")) flag_sd = atoi(tmpval);
      if(!strcmp(tmpstr, "flag_tt")) flag_tt = atoi(tmpval);
	  
	    if(!strcmp(tmpstr, "readfromrestart")) readfromrestart = atoi(tmpval);
	    if(!strcmp(tmpstr, "writetorestart")) writetorestart = atoi(tmpval);
	    if(!strcmp(tmpstr, "kdt")) k_dt = atof(tmpval);
	    if(!strcmp(tmpstr, "ngrainsinit")) ngrainsinit = atoi(tmpval);
	    if(!strcmp(tmpstr, "ngrainsor")) ngrainsor = atoi(tmpval);
	    if(!strcmp(tmpstr, "NUM_THREADS")) NUM_THREADS = atoi(tmpval);
      if(!strcmp(tmpstr, "toutput")) t_output = atoi(tmpval);
      if(!strcmp(tmpstr, "tprint")) t_print = atoi(tmpval);
      if(!strcmp(tmpstr, "tcheckpoint")) t_createcheckpoint = atoi(tmpval);
      
      if(!strcmp(tmpstr, "ngrainscontrol")) graincontrol = atoi(tmpval);
      if(!strcmp(tmpstr, "averageradius")) initaverageradius = atof(tmpval);
      if(!strcmp(tmpstr, "mob_ratio"))
      {
      	if(flag_anisotropy==4) sscanf(tempbuffer, "%20s : %f %f %f %f", tmpstr, &mob_ratio_array[0], &mob_ratio_array[1], &mob_ratio_array[2], &mob_ratio_array[3]);
      	else{mob_ratio = atof(tmpval);}
      } 
      if(!strcmp(tmpstr, "gbe_ratio")) gbe_ratio = atof(tmpval); //gbenergy/gbe_ratio = gb_min
      
	    if(!strcmp(tmpstr, "sigma_anisotropy")) sigma_anisotropy = atof(tmpval);
      //if(!strcmp(tmpstr, "mu_e0anisotropy")) mu_e0anisotropy = atof(tmpval);
      if(!strcmp(tmpstr, "sigma_e0anisotropy")) sigma_e0anisotropy = atof(tmpval);
      if(!strcmp(tmpstr, "dEpct")) dEpct = atof(tmpval);
      
      if(!strcmp(tmpstr, "seed")) randomseed = atoi(tmpval);
      if(!strcmp(tmpstr, "delta")) 
      {
      	if(flag_e0anisotropy==4) sscanf(tempbuffer, "%20s : %f %f %f %f", tmpstr, &sd_delta_array[0], &sd_delta_array[1], &sd_delta_array[2], &sd_delta_array[3]);
      	else{sd_delta = atof(tmpval);}
      	//sd_delta = atof(tmpval);
      }
    	if(!strcmp(tmpstr, "Vm")) Vm = atof(tmpval);
    
    
	  }
  }

  /*
  flag_anisotropy - mobility
  0 : No anisotropy  
  1 : Gauss anisotropy with full misorientation eeuler angle
  2 : Gauss anisotropy with tilt misorientation
  3 : mobility ratio with euler angle
  4 : Multiple mobility - 4 types
  */

  /*
  flag_e0anisotropy - segregation
  0 : No anisotropy
  1 : Gauss anisotropy
  2 : Cutoff anisotropy with cutoff euler angle theta<thetaC = Emin, theta>thetaC =Emax
  3 : Cutoff anisotropy with cutoff euler angle theta>thetaC = Emin, theta<thetaC =Emax 
  4 : special case of 4 anisotropic properties given by an array of 4 theta c (mu_anisotropy, and corresponding 4 E0 for the grain boundaries
  */
  //if(flag_e0anisotropy==4) sd_delta = sd_delta_array[0];
  if(flag_e0anisotropy==4)
  {
  	for(int i=0; i<4; i++) E0_array[i] *=1000.0;
  	E0 = E0_array[0];
  	iter_E0 = E0;
  	sd_delta = sd_delta_array[0];
  }
  else {E0 *= 1000;iter_E0 = E0;}
  if((abs(E0)>0.0001 && flag_sd!=2) || flag_e0anisotropy) {printf("!!!!!!!SOLUTE DRAG ACTIVATED!!!!!!!\n");flag_sd=1;}
  
  printf("E0 : %lf\n", E0);
  //printf("mob_ratio : %lf\n", mob_ratio);
  // E0 is defined in kJ/mol in input file, so convert in J/mol
  fclose(fp);
  if(flag_anisotropy) printf("Mobility Anisotropy activated!\n");
  if(flag_anisotropy==1) {printf("Gauss anisotropy with full misorientation space\n");printf("Mean: %lf, STD: %lf\n", mu_anisotropy, sigma_anisotropy);}
  if(flag_anisotropy==2) {printf("Gauss anisotropy with tilt angle only space\n");printf("Mean: %lf, STD: %lf\n", mu_anisotropy, sigma_anisotropy);}
  if(flag_anisotropy==3) {printf("Mobility ratio with euler angle and cutoff");}
  if(flag_anisotropy==4){for(int i=0; i<4; i++) printf("%f %f\n",mu_anisotropy_array[i], mob_ratio_array[i]); printf("\n");}
  if(flag_e0anisotropy) printf("Segregation anisotropy activated!\n");
  if(flag_e0anisotropy==1) {printf("Gauss anisotropy with full misorientation\n");printf("Mean: %lf, STD: %lf\n", mu_e0anisotropy, sigma_e0anisotropy);}
  if(flag_e0anisotropy==2||flag_e0anisotropy==3){printf("dE anisotropy with cutoff\n");printf("cutoff : %lf, below cutoff : %lf, above cutoff : %lf\n", mu_e0anisotropy, anisotropic_segregation(E0, mu_e0anisotropy-1, flag_e0anisotropy), anisotropic_segregation(E0, mu_e0anisotropy+1, flag_e0anisotropy));}
  if(flag_e0anisotropy==4){for(int i=0; i<4; i++) printf("%f %f %e\n",E0_array[i], mu_e0anisotropy_array[i], sd_delta_array[i]);}
  if(flag_e0anisotropy){printf("TEST SEG\n");printf("E(10): %f, E(45): %f, E(60): %f\n", anisotropic_segregation(E0, 10, flag_e0anisotropy), anisotropic_segregation(E0, 45, flag_e0anisotropy), anisotropic_segregation(E0, 60, flag_e0anisotropy));}
  if(flag_gbeanisotropy){printf("GB energy anisotropy activated\n GBE_Max: %lf, GBE_Min: %lf\n", gbenergy, gbenergy/gbe_ratio);}
  if(flag_sd==2){printf("Reading dG-v relationship from file\n"); readSDinfo();}
  if(flag_tt==1){printf("Reading time-temperature relationship\n"); readTTinfo();}
  return;
}


void nondimensionalize()
{
  extern float c_length, c_energy, c_time, c_mobility, c_diffusivity;
  extern float R, T, Vm;
  extern float gbe_ratio;
  extern float eps, omega, dx, Tref;
  extern float eps_min, eps_max, gamma_min, gamma_max;
  extern float gbenergy, gbpoints;
  extern float real_dx, real_totaltime;

  extern float time, dt, k_dt, sim_time;
  extern int timesteps;

  extern float m0, qm;
  extern float d0, qd;
  extern float dalpha;
  extern float mob;
  extern float sd_delta;
  extern float beta, beta_scaled;
  extern float beta_table[2];

  extern int timesteps;

  float real_eps, real_omega;
  float gbe_init = 0.5*gbenergy*(1+1/gbe_ratio);
  float gamma_init = 1.5;
  float eps_init, omega_init;

  R = 8.314;
  //Vm = 7.1e-6;
  //Vm = 10.21e-6;
  //Tref = 1073.0;
  Tref = T;
  c_diffusivity = d0*exp(-qd/(R*T));
  //c_diffusivity = 1e-15;
  //c_diffusivity = 1;

  c_length = real_dx;

  dx = real_dx/c_length;
  printf("dx : %lf\n", dx);
  c_energy = R*Tref/Vm;
  //c_energy = 1.0;
  c_time = c_length*c_length/(c_diffusivity);
  c_mobility = 1/(c_energy*c_time);

  beta_scaled = beta/(Vm*c_energy);
  beta_table[0] = 0;
  beta_table[1] = beta_scaled;

  printf("C_energy : %e\n", c_energy);
  printf("C_time : %e\n", c_time);
  
  
  
  //real_eps = sqrt(3*gbpoints*c_length*gbenergy)/2; //#Units - sqrt(J/m)
  //real_omega = 6*gbenergy/(gbpoints*c_length); //#Units - J/m^3
  
  real_eps = sqrt(3*gbpoints*c_length*gbe_init)/2; //#Units - sqrt(J/m)
  real_omega = 6*gbe_init/(gbpoints*c_length); //#Units - J/m^3
 // if(flag_gbeanisotropy==0){eps_min = eps_init; gamma_min=gamma_init;eps_max = eps_init; gamma_max=gamma_init;}
  
  
  eps_init = real_eps/(c_length*sqrt(c_energy));
  omega_init = real_omega/(c_energy);
  printf("gamma_init: %lf\n",1.0/(103.397*pow(1/6.0, 6) - 165.393*pow(1/6.0, 5) + 105.3469*pow(1/6.0, 4) - 44.55661*pow(1/6.0, 3) + 24.7348*pow(1/6.0, 2) -11.25718*1/6.0 + 1.999642));
  if(flag_gbeanisotropy==0){eps_min = eps_init; gamma_min=gamma_init;eps_max = eps_init; gamma_max=gamma_init;}
  float G1 = (gbenergy/gbe_ratio)/(6.0*gbe_init);
  printf("sigma_min: %lf, G1: %lf\n", gbenergy/gbe_ratio, G1);
  gamma_min = 1.0/(103.397*pow(G1, 6) - 165.393*pow(G1, 5) + 105.3469*pow(G1, 4) - 44.55661*pow(G1, 3) + 24.7348*pow(G1, 2) -11.25718*G1 + 1.999642);
  double f0c_min = pow((-0.072966*pow((1.0/gamma_min), 5) + 0.35784*pow((1.0/gamma_min), 4) - 0.68325*pow((1.0/gamma_min), 3) + 0.63578*pow((1.0/gamma_min), 2) - 0.48566*((1.0/gamma_min)) + 0.53703), 2.0);
  double f0c_init = pow((-0.072966*pow((1.0/gamma_init), 5) + 0.35784*pow((1.0/gamma_init), 4) - 0.68325*pow((1.0/gamma_init), 3) + 0.63578*pow((1.0/gamma_init), 2) - 0.48566*(1.0/gamma_init) + 0.53703), 2);
  eps_min = sqrt(6.0*f0c_min*gbe_init*gbpoints*c_length)/(c_length*sqrt(c_energy));
  //eps_min = eps_init*sqrt(f0c_min/f0c_init);
  
  float G2 = (gbenergy)/(6*gbe_init);
  gamma_max = 1/(103.397*pow(G2, 6) - 165.393*pow(G2, 5) + 105.3469*pow(G2, 4) - 44.55661*pow(G2, 3) + 24.7348*pow(G2, 2) -11.25718*G2 + 1.999642);
  double f0c_max = pow((-0.072966*pow(1.0/gamma_max, 5) + 0.35784*pow(1.0/gamma_max, 4) - 0.68325*pow(1.0/gamma_max, 3) + 0.63578*pow(1.0/gamma_max, 2) - 0.48566*(1.0/gamma_max) + 0.53703), 2);
  //eps_max = eps_init*sqrt(f0c_max/f0c_init);
  eps_max = sqrt(6*f0c_max*gbe_init*gbpoints*c_length)/(c_length*sqrt(c_energy));
  
  eps = eps_max;
  omega = omega_init;
  
  printf("Kappa: %le, m: %le\n", real_eps*real_eps, real_omega);
  printf("FOC: init: %lf, max: %lf, min: %lf\n", f0c_init, f0c_max, f0c_min);
/*
  eps_max = eps;
  gamma_max = 1.5;
  if(flag_gbeanisotropy==0){eps_min = eps; gamma_min=1.5;}
  float G = 1.0/(6*gbe_ratio);
  // from Minar 2022
  gamma_min = 1/(103.397*pow(G, 6) - 165.393*pow(G, 5) + 105.3469*pow(G, 4) - 44.55661*pow(G, 3) + 24.7348*pow(G, 2) -11.25718*G + 1.999642);
  float f0c_min = pow(0.072966*pow(1.0/gamma_min, 5) + 0.35784*pow(1.0/gamma_min, 4) - 0.68325*pow(1.0/gamma_min, 3) + 0.63578*pow(1.0/gamma_min, 2) - 0.48566*(1.0/gamma_min) + 0.53703, 2);
  float f0c_max = pow(0.072966*pow(1.0/gamma_max, 5) + 0.35784*pow(1.0/gamma_max, 4) - 0.68325*pow(1.0/gamma_max, 3) + 0.63578*pow(1.0/gamma_max, 2) - 0.48566*(1.0/gamma_max) + 0.53703, 2);
  eps_min = eps_max*(sqrt(f0c_min/f0c_max));
 */
  printf("EPS: init: %lf, max: %lf, min: %lf\n", sqrt(6*f0c_init*gbe_init*gbpoints*c_length)/(c_length*sqrt(c_energy)), eps_max, eps_min);
  printf("gamma: init: %lf, max: %lf, min: %lf\n", gamma_init, gamma_max, gamma_min);
  
  printf("Parameter : Actual Values Used (ND) : Values (Real) \n");
  printf("epsilon : %e : %e\n", eps, eps*c_length*sqrt(c_energy));
  printf("omega : %e : %e\n", omega, omega*c_energy);
	printf("Dx : %e : %e\n", dx, dx*c_length);
  printf("Total time : %le\n", real_totaltime);

  //printf("time : %le", time);

  evaluate_constants();  
  if(sd_delta<=0) sd_delta = 2*eps/sqrt(omega/2)*c_length;
  
  //dt = fmin(k_dt*dx*dx/dalpha, k_dt*dx*dx/(mob*eps*eps));

  printf("timescale: %le dt : %le : %le\n", c_time, dt, dt*c_time);
  
  if(flag_anisotropy) printf("Cutoff: %le, M(before cutoff): %le, M(after cutoff): %le\n", mu_anisotropy, anisotropic_mobility(mob, mu_anisotropy-1, flag_anisotropy), anisotropic_mobility(mob, mu_anisotropy+1, flag_anisotropy));
  printf("Mobility : %e : %e c_mobility : %e\n", mob, mob*c_mobility, c_mobility);
  printf("beta (in J/mol) : %le, beta_scaled:%le\n", beta, beta_scaled);
  printf("If in 1D , expected velocity : %le\n", beta*(mob_int(m0, qm, T)));
  if(flag_anisotropy) {printf("TEST MOB!!!!!\n");printf("M(10): %le, M(45): %le, M(60): %le\n", anisotropic_mobility(mob, 10, flag_anisotropy), anisotropic_mobility(mob, 45, flag_anisotropy), anisotropic_mobility(mob, 60, flag_anisotropy));}
  
  return;	
}

int initializegeometry(int flag)
{

  extern int nx, ny, nz;
  extern int dim;
  extern float *grid, *newgrid, *curvature;
  extern float *dphidtold;
  extern float **orientation, **rotationmatrix;
  extern int *flaggrid, *newflaggrid;
  extern int flag_anisotropy, flag_e0anisotropy, flag_gbeanisotropy;
  extern float *c, *newc;
  extern int N, M;
  extern int ngrainsmax, ngrainsor;;
  extern float vf, c_init;
  extern int trackid;
  extern int *activegrain;
  extern int ngrainsinit;
  extern float initaverageradius, vf;
  extern float sim_time;
  extern int NUM_THREADS;
  extern float gbenergy, gbe_ratio;
  extern float eps, omega;
  float thetamax = 90;
  float thetamin = 0;
  float thetawindow = thetamax - thetamin;
  int r=0;

  extern int grid2, grid3;

  grid2 = nx*ny;
  grid3 = nx*ny*nz;


  int index;

  N = grid3*ngrainsmax;
  M = grid3;

  /*Allocate*/
  grid = malloc(N * sizeof(float));
  newgrid = malloc(N * sizeof(float));
  dphidtold = malloc(N *sizeof(float));

  c = malloc(M * sizeof(float));
  newc = malloc(M * sizeof(float));

  flaggrid = malloc(N * sizeof(int));
  newflaggrid = malloc(N * sizeof(int));


  curvature=malloc(N * sizeof(float));
  activegrain = malloc(ngrainsmax * sizeof(int));

  uniformrandomnumbergenerator();
    

  /*Initialize*/
  for(int i=0; i<nx; i++)
  {
  	for(int j=0;j<ny;j++)
  	{
  	  for(int k=0; k<nz; k++)
  	  {
  	  	index = k*grid2 + j*nx + i;
  		  c[index] = c_init; newc[index] = c_init; curvature[index]=0.0;
  		  for(int l=0; l<ngrainsmax; l++)
  		  {
  			grid[l*grid3 + index] = 0.0; newgrid[l*grid3 + index] = 0.0;flaggrid[l*grid3 + index] = -1;newflaggrid[l*grid3 + index] = -1;dphidtold[l*grid3+index]=0.0;
  		  }
      }
  	}
  }

  /*Setup*/
  if(flag==0)
  {
  	// Read from previous file
  	char *readfile;
  	readfile = "output/restartfile.txt";
  	FILE *fr = fopen(readfile, "r");
  	int xx,  count;
  	fscanf(fr, "%e\n", &sim_time);
  	for(int i=0; i<nx; i++)
    {
  	  for(int j=0;j<ny;j++)
  	  {
  	    for(int k=0; k<nz; k++)
  	    {
  	      count = 1;
  	      index = k*grid2 + j*nx + i;
  	      fscanf(fr, "%d %d %d", &xx, &xx, &xx);
  	      fscanf(fr, " %d", &count);
  	      for(int l=0; l<count; l++)
  	      {
  	      	fscanf(fr, " %d %e", &flaggrid[l*grid3 + index], &grid[l*grid3 + index]);
  	      }
  	    }
  	  }
  	}
  	fclose(fr);
  	
  }
  else if(flag==1)
  {
  	trackid = (int)(nx/2)-4;
  	//Flat geometry
  	for(int i=0; i<nx; i++)
    {
  	  for(int j=0;j<ny;j++)
  	  {
  	    for(int k=0; k<nz; k++)
  	    {
  	      index = k*grid2 + j*nx + i;
  	      if(i<(int)(vf*nx)) 
  	      {
  	    	grid[0*grid3 + index] = 1.0;
  	    	flaggrid[0*grid3 + index] = 0;
  	      }
  	      else
  	      {
  	      	grid[0*grid3 + index] = 1.0;
  	    	flaggrid[0*grid3 + index] = 1;
  	      }
  	    }
  	  }
  	}

  }
  else if(flag==2)
  {

  	//circular geometry
  	if(dim == 1){ printf("Change the dimension > 2 for circular shrinkage"); exit(0);}
  	int rad, imid=(int)(nx/2), jmid=(int)(ny/2), kmid=(int)(nz/2);
  	if(dim==2) rad = sqrt(vf*grid3/(PI));
  	else if(dim==3) rad = pow(3*vf*grid3/(4*PI), 1/3.0);
  	trackid = (int)(imid - rad);
  	for(int i=0; i<nx; i++)
    {
  	  for(int j=0;j<ny;j++)
  	  {
  	    for(int k=0; k<nz; k++)
  	    {
  	      index = k*grid2 + j*nx + i;
  	      if((i-imid)*(i-imid) + (j-jmid)*(j-jmid) + (k-kmid)*(k-kmid) < rad*rad)
  	      {
  	    	  grid[0*grid3 + index] = 1.0;
  	    	  flaggrid[0*grid3 + index] = 0;
  	      }
  	      else
  	      {
  	      	grid[0*grid3 + index] = 1.0;
  	    	  flaggrid[0*grid3 + index] = 1;
  	      }
  	    }
  	  }
  	}
  }
  else if(flag==3)
  {
  	trackid = (int)(nx/2);
  	//three grains
  	if(dim == 1){ printf("Change the dimension > 2 for circular shrinkage"); exit(0);}
  	for(int i=0; i<nx; i++)
    {
  	  for(int j=0;j<ny;j++)
  	  {
  	    for(int k=0; k<nz; k++)
  	    {
  	      index = k*grid2 + j*nx + i;
  	      if(i<(int)(0.5*nx))
  	      {
  	      	if(j<(int)(0.5*nx))
  	      	{
  	      	  grid[0*grid3 + index] = 1.0;
  	    	    flaggrid[0*grid3 + index] = 0;
  	      	}
  	      	else
  	      	{
  	      	  grid[0*grid3 + index] = 1.0;
  	    	    flaggrid[0*grid3 + index] = 1;
  	      	}
  	      }
  	      else
  	      {
  	    	grid[0*grid3 + index] = 1.0;
  	    	flaggrid[0*grid3 + index] = 2;
  	      }
  	    }
  	  }
  	}
  	
  }
  else if(flag==4)
  {
  	trackid = (int)(nx/2);
  	//four grains
  	if(dim == 1){ printf("Change the dimension > 2 for circular shrinkage"); exit(0);}
  	for(int i=0; i<nx; i++)
    {
  	  for(int j=0;j<ny;j++)
  	  {
  	    for(int k=0; k<nz; k++)
  	    {
  	      index = k*grid2 + j*nx + i;
  	      if(i<(int)(0.5*nx))
  	      {
  	      	if(j<=(int)(0.5*nx))
  	      	{
  	      	  grid[0*grid3 + index] = 1.0;
  	    	    flaggrid[0*grid3 + index] = 0;
  	      	}
  	      	else
  	      	{
  	      	  grid[0*grid3 + index] = 1.0;
  	    	    flaggrid[0*grid3 + index] = 1;
  	      	}
  	      }
  	      else
  	      {
  	      	if(j<(int)(0.5*nx))
  	      	{
  	      	  grid[0*grid3 + index] = 1.0;
  	    	    flaggrid[0*grid3 + index] = 2;
  	      	}
  	      	else
  	      	{
  	      	  grid[0*grid3 + index] = 1.0;
  	    	    flaggrid[0*grid3 + index] = 3;
  	      	}
  	      }
  	    }
  	  }
  	}
  }
  else
  {
  	if(dim<2) {printf("Not applicable for 1 dimension, change the dimension to 2 or 3\n"); exit(0);}
    if(graincontrol!=0)
    {
      if(initaverageradius<4) {printf("Either change the grain control to 0 or increase the average radius");exit(0);}
      if(dim==2) {ngrainsinit = (int)((nx*ny)/(PI*initaverageradius*initaverageradius));}
      if(dim==3) {ngrainsinit = (int)((nx*ny*nz)/(4/3*PI*initaverageradius*initaverageradius*initaverageradius));}
    }
    if(ngrainsinit < ngrainsor) ngrainsor = ngrainsinit;
    printf("Number of initialized grains : %d Number of misorientation sample: %d\n", ngrainsinit, ngrainsor);
  	int xp[ngrainsinit], yp[ngrainsinit], zp[ngrainsinit]; //Centre of each grain
  	float xdis, ydis, zdis=0;
  	float dmin, d;
  	int a,  l, aindex;
  	for(int i=0;  i<ngrainsinit; i++)
  	{
  	  xp[i]=(gsl_rng_uniform(URN))*nx;
  	  yp[i]=(gsl_rng_uniform(URN))*ny;
  	  if(dim==3) zp[i]=(gsl_rng_uniform(URN))*nz;
  	  if(i==0) printf("xp : %d, yp: %d\n", xp[i], yp[i]);
  	}
  	// Voronoi Tessellation
//  	#pragma omp parallel for collapse(3) num_threads(NUM_THREADS) schedule(auto) private(index, dmin, a, aindex, xdis, ydis, d, l) shared(grid, flaggrid)
  	for(int i=0; i<nx; i++)
  	{
      for(int j=0;j<ny; j++)
  	  {
  	  	for(int k=0; k<nz; k++)
  	  	{
  	  		if(i%100==0 && j==0 && k==0) printf("i:%d\n", i);
  	  	  index=k*grid2 + j*nx + i;
  	  	  a=-1;
  	  	  aindex = 0;
  	  	  dmin = pow(((nx-1)*(nx-1) + (ny-1)*(ny-1)),0.5);
  		    if(dim==3) dmin = pow(((nx-1)*(nx-1) + (ny-1)*(ny-1) + (nz-1)*(nz-1)),0.5);
  	  	  for(l=0;l<ngrainsinit; l++)
  	  	  {
  	  	  	xdis = xp[l]-i;
  	  	  	ydis = yp[l]-j;
  	  	  	if(dim==3)zdis = zp[l]-k;
  	  	  	if(bc!=0)
  	  	  	{
  	  	  		if((xdis)>(int)(nx/2)) xdis -= nx;
  	  	  		if((xdis)<-(int)(nx/2)) xdis += nx;
  	  	  		if((ydis)>(int)(ny/2)) ydis -= ny;
  	  	  		if((ydis)<-(int)(ny/2)) ydis += ny;
  	  	  		if(dim==3)
  	  	  		{
  	  	  		  if((zdis)>(int)(nz/2)) zdis -= nz;
  	  	  		  if((xdis)<-(int)(nx/2)) zdis += nz;	
  	  	  		}
  	  	  	}
  	  	  	d = pow(((xdis)*(xdis) + (ydis)*(ydis)),0.5);
  			    if(dim==3) dmin = pow(((xdis)*(xdis) + (ydis)*(ydis) + (zdis)*(zdis)),0.5);
  			    if(d<dmin)
  			    {
  			      dmin = d;
  			      a = l;
  			      aindex = index;
  			      //printf("grain num: %d, index: %d\n", l, aindex);
  			    }
  	  	  }
  	  	  grid[aindex]=1.0;
  	  	  flaggrid[aindex]=a;
  	  	}

  	  }

  	}

    
  	//Voronoi tesselation
  	//return 0;

  }
  if(flag_anisotropy || flag_e0anisotropy || flag_gbeanisotropy)
  {
    // 3 Euler angles + 9 rotation matrix
    printf("Distributing Orientation\n");
    // Allocate space for each grain in the system to have an orientation variable
    orientation = (float **)malloc((ngrainsinit) * sizeof(float *)); 
    for(int i=0; i<ngrainsinit; i++) orientation[i] = (float *)malloc(3* sizeof(float));

    masterorientation = (float **)malloc((ngrainsor) * sizeof(float *)); 
    for(int i=0; i<ngrainsor; i++) masterorientation[i] = (float *)malloc(3* sizeof(float));

    // Create a master list with random sampling of orientation for ngrainsmax 
    for(int i=0; i<ngrainsor; i++)
    {
    	for(int j=0; j<3; j++)
    	{
    		masterorientation[i][j] = (gsl_rng_uniform(URN))*thetawindow + thetamin;
    	}
    }
    // Assign the orientation to each grain randomly from the master list
    for(int i=0; i<ngrainsinit; i++)
    {
    	//r = (int)((gsl_rng_uniform(URN))*ngrainsmax);
    	r = (i%ngrainsor);
    	orientation[i][0]=masterorientation[r][0];
    	orientation[i][1]=masterorientation[r][1];
    	orientation[i][2]=masterorientation[r][2];
    }	
    


    rotationmatrix=(float **)malloc((ngrainsor) * sizeof(float *)); // 9 rotation matrix 00, 01, 02, 10, 11, 12, 20, 21, 22
    for(int i=0; i<ngrainsor; i++) rotationmatrix[i] = (float *)malloc(9* sizeof(float));

    printf("Calculating Rotation Matrix and storing \n");
    float t1, t2, t3;
    for(int i=0; i<ngrainsor; i++)
    {
      t1=masterorientation[i][0]*RADIAN;
      t2=masterorientation[i][1]*RADIAN;
      t3=masterorientation[i][2]*RADIAN;
      
      rotationmatrix[i][0] = cos(t1)*cos(t2) - sin(t1)*sin(t2)*cos(t3);
      rotationmatrix[i][1] = sin(t1)*cos(t2) + cos(t1)*sin(t2)*cos(t3);
      rotationmatrix[i][2] = sin(t2)*sin(t3);
      rotationmatrix[i][3] = -cos(t1)*sin(t2) - sin(t1)*cos(t2)*cos(t3);
      rotationmatrix[i][4] = -sin(t1)*sin(t2) + cos(t1)*cos(t2)*cos(t3);
      rotationmatrix[i][5] = cos(t2)*sin(t3);
      rotationmatrix[i][6] = sin(t1)*sin(t3);
      rotationmatrix[i][7] = -cos(t1)*sin(t3);
      rotationmatrix[i][8] = cos(t3);
    }
	
    populatemisorientation();
    if(flag==3) {
      float eps_01, eps_12, eps_20;
      float gamma_01, gamma_12, gamma_20;
      float f0c_01, f0c_12, f0c_20;
      float gbe_init = 0.5*gbenergy*(1+1/gbe_ratio);
      eps_01 = anisotropic_eps(eps, evaluatemisorientation(0, 1), flag_gbeanisotropy); eps_12 = anisotropic_eps(eps, evaluatemisorientation(1, 2), flag_gbeanisotropy); eps_20 = anisotropic_eps(eps, evaluatemisorientation(2, 0), flag_gbeanisotropy);
      gamma_01  = anisotropic_gamma(gamma_max, evaluatemisorientation(0, 1), flag_gbeanisotropy); gamma_12  = anisotropic_gamma(gamma_max, evaluatemisorientation(1, 2), flag_gbeanisotropy); gamma_20  = anisotropic_gamma(gamma_max, evaluatemisorientation(2, 0), flag_gbeanisotropy);
      printf("GBE_init: %lf, Gamma list: %lf, %lf, %lf\n", gbe_init, gamma_01, gamma_12, gamma_20);
      printf("Misorientation list: %lf, %lf, %lf\n", evaluatemisorientation(0, 1), evaluatemisorientation(1, 2), evaluatemisorientation(2, 0));
      f0c_01 = pow((-0.072966*pow(1.0/gamma_01, 5) + 0.35784*pow(1.0/gamma_01, 4) - 0.68325*pow(1.0/gamma_01, 3) + 0.63578*pow(1.0/gamma_01, 2) - 0.48566*(1.0/gamma_01) + 0.53703), 2);
      f0c_12 = pow((-0.072966*pow(1.0/gamma_12, 5) + 0.35784*pow(1.0/gamma_12, 4) - 0.68325*pow(1.0/gamma_12, 3) + 0.63578*pow(1.0/gamma_12, 2) - 0.48566*(1.0/gamma_12) + 0.53703), 2);
      f0c_20 = pow((-0.072966*pow(1.0/gamma_20, 5) + 0.35784*pow(1.0/gamma_20, 4) - 0.68325*pow(1.0/gamma_20, 3) + 0.63578*pow(1.0/gamma_20, 2) - 0.48566*(1.0/gamma_20) + 0.53703), 2);
  
      float gbe_01 = (8.0)*(f0c_01)*gbe_init;
      float gbe_12 = (8.0)*(f0c_12)*gbe_init;
      float gbe_20 = (8.0)*(f0c_20)*gbe_init;
      
      printf("Mobility: Grain 12: %lf, Grain 23: %lf, Grain 31: %lf\n", 
    	anisotropic_mobility(mob, evaluatemisorientation(0, 1), flag_anisotropy), anisotropic_mobility(mob, evaluatemisorientation(1, 2), flag_anisotropy), 
    	anisotropic_mobility(mob, evaluatemisorientation(2, 0), flag_anisotropy));
    	//printf("Energy: Grain 12: %lf, Grain 23: %lf, Grain 31: %lf\n", 
    	//gbenergy*anisotropic_eps(eps, evaluatemisorientation(0, 1), flag_gbeanisotropy)*sqrt(anisotropic_omega(omega, evaluatemisorientation(0, 1), flag_gbeanisotropy))/(eps*sqrt(omega)),
    	//gbenergy*anisotropic_eps(eps, evaluatemisorientation(1, 2), flag_gbeanisotropy)*sqrt(anisotropic_omega(omega, evaluatemisorientation(1, 2), flag_gbeanisotropy))/(eps*sqrt(omega)),
    	//gbenergy*anisotropic_eps(eps, evaluatemisorientation(2, 0), flag_gbeanisotropy)*sqrt(anisotropic_omega(omega, evaluatemisorientation(2, 0), flag_gbeanisotropy))/(eps*sqrt(omega)));
    	printf("Energy: Grain 12: %lf, Grain 23: %lf, Grain 31: %lf\n", gbe_01, gbe_12, gbe_20);
    }
/*
    printf("Orientation list\n");
    FILE *fmis = fopen("orientationlist.txt", "w");
    for(int i=0; i<ngrainsinit; i++) 
    {
      printf("Grain num : %d, orientation : %lf %lf %lf\n", i, orientation[i][0], orientation[i][1], orientation[i][2]);
      fprintf(fmis, "%d %lf %lf %lf\n", i, orientation[i][0], orientation[i][1], orientation[i][2]);      
    }
    fclose(fmis);

    printf("Misorientation table\n");
    for(int i=0; i<ngrainsinit; i++)
    {
      for(int j=0; j<i; j++)
      {
        printf("%d %d : %lf\n", i, j, evaluatemisorientation(i, j));
      } 
    }
*/
  }

  
  
  for(int i=0; i<N; i++) {newgrid[i]=grid[i];newflaggrid[i]=flaggrid[i];}
  //for(int i=0; i<nx; i++) {for(int j=0; j<ny; j++) {for(int k=0; k<nz; k++) {index=k*grid2+j*nx+i; printf("index : %d, phi: %lf\n", index, grid[index+0*grid3]);}}}
  return 1;
}

float mob_int(float m0, float Q, float T)
{
	extern float R;
	return m0*exp(-Q/(R*T));
}

float mob_nd(float mob_int)
{
	extern float c_mobility, c_length, c_energy;
	extern float eps, omega;
	extern float sim_time, gbenergy, gbe_ratio;
	
	float gbe_init = 0.5*gbenergy*(1+1/gbe_ratio);
	
	float w = 6*gbe_init/(omega*c_energy);
	float se = gbe_init;
	
  	//float w = 2*eps/sqrt(omega/2);
	//float se = (eps/3)*sqrt(2*omega);
  	if(sim_time==0)
  	{
   	 printf("Surface Energy : %le\n", se);
   	 printf("Interface Width : %le\n", w);
  	}
  	
 
  	return (4*mob_int/(3*w))/c_mobility;
  	//return (mob_int/((3.0/4.0)*w*c_length))/c_mobility;
  	
	//return (mob_int/(3*eps/sqrt(omega/2)*c_length))/c_mobility;

	//return mob_int*0.1057/(c_length*c_mobility); /*the factor 0.1057 is corresponding to integral (dphi/dx)^2 dx - as observed from multiple simulations*/
}

void freememory()
{
	extern float *grid, *newgrid,*c, *newc;
  extern float *dphidtold;
	extern int *flaggrid, *newflaggrid;
	extern int *activegrain;
	extern float *curvature;
	extern FILE *fs;
  extern int flag_anisotropy;
  extern int ngrainsinit;
  extern gsl_rng *URN;
  extern flag_sd;
  gsl_rng_free(URN);

	free(grid);
	free(newgrid);
  free(dphidtold);
	free(c);
	free(newc);
	free(flaggrid);
	free(newflaggrid);
	free(activegrain);
	free(curvature);
  if(flag_anisotropy||flag_e0anisotropy) freemisorientation();
	if(flag_sd==2) free_sdinfo();
  fclose(fs);

	printf("Allocated memory is now free!!\n");
	return;
}

void uniformrandomnumbergenerator()
{
  extern int randomseed;
 	const gsl_rng_type * URNG;
	extern gsl_rng *URN;
	gsl_rng_env_setup();
  
	URNG = gsl_rng_default;
	URN = gsl_rng_alloc(URNG);
	gsl_rng_set(URN,randomseed);
//	gsl_rng_set(URN);
}


void readSDinfo()
{
  extern float **sdinfo;
  extern int sdcount;
  extern float sdmin, sdmax;
  
  char *sdfile;
  sdfile = "inp/sdinfo.txt";
  FILE *sd = fopen(sdfile, "r");
  fscanf(sd, "%d %f %f\n", &sdcount, &sdmin,&sdmax);

  sdinfo = (float **)malloc((sdcount) * sizeof(float *)); 
  for(int i=0; i<sdcount; i++) sdinfo[i] = (float *)malloc(2* sizeof(float));

  for(int i=0; i<sdcount; i++)
  {
    fscanf(sd, "%e %e\n", &sdinfo[i][0], &sdinfo[i][1]);
  }
  fclose(sd);
  printf("Datapoints : %d, min : %f, max: %f\n", sdcount, sdmin, sdmax);
  for(int i=0; i<sdcount;i++) printf("%d Velocity : %e, dg: %e\n", i, sdinfo[i][0], sdinfo[i][1]);
  return;
  
}

void free_sdinfo()
{
  extern float **sdinfo;
  extern int sdcount;
  for(int i=0; i<sdcount; i++) free(sdinfo[i]);
  return;
}

void readTTinfo()
{
  extern float threshold, real_totaltime;

  extern float **ttinfo;
  //ttinfo has ttstart-time ttstart-Temperature ttstop-time ttstop-Temperature 
  extern int ttcount;
  //ttcount has the number of steps used in the simulation
  
  char *ttfile;
  ttfile = "inp/ttinfo.txt";
  FILE *tt = fopen(ttfile, "r");
  fscanf(tt, "%d\n", &ttcount);

  ttinfo = (float **)malloc((ttcount) * sizeof(float *)); 
  for(int i=0; i<ttcount; i++) ttinfo[i] = (float *)malloc(4* sizeof(float));

  for(int i=0; i<ttcount; i++)
  {
    fscanf(tt, "%e %e %e %e\n", &ttinfo[i][0], &ttinfo[i][1], &ttinfo[i][2], &ttinfo[i][3]);
  }
  fclose(tt);
  printf("Datapoints : %d\n", ttcount);
  for(int i=0; i<ttcount;i++) printf("Step - %d , time-initial : %e,  time-final : %e, Temperature-initial : %e, Temperature-final: %e\n", i, ttinfo[i][0], ttinfo[i][2], ttinfo[i][1], ttinfo[i][3]);
  
  //Check if the time is consistent with the total time steps
  if(fabs(real_totaltime-ttinfo[ttcount-1][2]) > threshold) 
  {  
    printf("Inconsistent time temperature profile with total time in input file\n PLEASE keep the last time in ttinfo equal to the total time\n"); 
    printf("%e, %e, Different in time : %e\n", real_totaltime, ttinfo[ttcount-1][2], real_totaltime-ttinfo[ttcount-1][2]);
    exit(0);
  }  

  return;
  
}
