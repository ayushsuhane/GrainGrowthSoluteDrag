#ifndef GLOBAL_H
#define GLOBAL_H

#define PI 3.14159
#define RADIAN 0.01745
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix_double.h"

extern gsl_rng *URN;
extern gsl_rng *PDF;


extern int NUM_THREADS;

/************************************/
/************************************/
extern int dim, nx, ny, nz;
extern float m0, qm, d0, qd;
extern float eps, omega, E0, T; 
extern float c_init, vf, dalpha;
extern float R, Vm, Tref;
extern float gbenergy, gbpoints;
extern float beta;
extern float mob;
extern float averagecurvature;
extern float sd_delta;
extern float aver_radius;
extern float mob_ratio, gbe_ratio;
extern float maxratio;
extern float eps_min, eps_max, gamma_min, gamma_max;

extern float real_dx, real_omega, real_eps, real_totaltime, sim_time;
extern float dx, totaltime, dt, time;
extern float E0;
extern int bc; 
extern int timesteps;
extern int ngrainsmax;
extern int N, M;
extern int vel_count;

extern float mu_anisotropy, sigma_anisotropy, mu_e0anisotropy, sigma_e0anisotropy;
extern float E0_array[4], mu_e0anisotropy_array[4], mu_anisotropy_array[4], sd_delta_array[4], mob_ratio_array[4];

extern float dEpct;
extern int randomseed;

extern float c_length, c_energy, c_time, c_mobility, c_diffusivity;
extern float beta_scaled;
extern float beta_table[2];
extern float k_dt;
extern float threshold, p_threshold;
extern float *grid, *newgrid, *c, *newc, *curvature;
extern float *dphidtold;
extern float **orientation, **rotationmatrix;
extern float **masterorientation;
extern int  *flaggrid, *newflaggrid;
extern int *activegrain;

extern int grid2, grid3;
extern int flag_geometry, flag_curvature, flag_anisotropy, flag_e0anisotropy, flag_gbeanisotropy;
extern int flag_sd, flag_tt;
extern int t_output, t_print, t_createcheckpoint;
extern int graincontrol;
extern float initaverageradius;
extern float expected_velocity, pvelocity;

extern int readfromrestart, writetorestart;
extern int trackid;
extern int numbergrains, ngrainsinit, ngrainsor;


extern float **fullmisorientation;
extern int symm_length;
extern float **symmetryoperator;


extern float **sdinfo;
extern float sdmin, sdmax;
extern int sdcount;

extern float **ttinfo;
extern int ttcount;

extern float iter_E0;
extern FILE *fs;

#endif
