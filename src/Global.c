#include "Global.h"
#include "functions.h"

gsl_rng *URN;
gsl_rng *PDF;


int NUM_THREADS;

/************************************/
/************************************/
int dim, nx, ny, nz;
float m0, qm, d0, qd;
float eps, omega, E0, T; 
float c_init, vf, dalpha;
float R, Vm, Tref;
float gbenergy, gbpoints;
float beta;
float mob;
float averagecurvature;
float sd_delta;
float aver_radius;
float mob_ratio, gbe_ratio;
float maxratio;

float real_dx, real_omega, real_eps, real_totaltime, sim_time;
float eps_min, eps_max, gamma_min, gamma_max;

float dx, totaltime, dt, time;
float E0;
int bc; 
int timesteps;
int ngrainsmax;
int N, M;
int vel_count;

float mu_anisotropy, sigma_anisotropy, mu_e0anisotropy, sigma_e0anisotropy;
float E0_array[4], mu_e0anisotropy_array[4], mu_anisotropy_array[4], sd_delta_array[4], mob_ratio_array[4];

float dEpct;
int randomseed;

float c_length, c_energy, c_time, c_mobility, c_diffusivity;
float beta_scaled;
float beta_table[2];
float k_dt;
float threshold, p_threshold;
float *grid, *newgrid, *c, *newc, *curvature;
float *dphidtold;
float **orientation, **rotationmatrix;
float **masterorientation;
int  *flaggrid, *newflaggrid;
int *activegrain;

int grid2, grid3;
int flag_geometry, flag_curvature, flag_anisotropy, flag_e0anisotropy, flag_gbeanisotropy;
int flag_sd, flag_tt;
int t_output, t_print, t_createcheckpoint;
int graincontrol;
float initaverageradius;
float expected_velocity, pvelocity;

int readfromrestart, writetorestart;
int trackid;
int numbergrains, ngrainsinit, ngrainsor;


float **fullmisorientation;
int symm_length;
float **symmetryoperator;


float **sdinfo;
float sdmin, sdmax;
int sdcount;

float **ttinfo;
int ttcount;

float iter_E0;
FILE *fs;

