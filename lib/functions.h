#ifndef FUNCTIONS_H
#define FUNCTIONS_H


// Setup functions
void simulationsetup();
void readparamsfromfile();
void nondimensionalize();
int initializegeometry();
float mob_int(float, float, float);
float mob_nd(float);
void freememory();
void uniformrandomnumbergenerator();
void readSDinfo();
void free_sdinfo();

void readTTinfo();



// SOlver functions
void pfsolver(int t);
void evaluate_constants();
float calc_temperature();
int pos2ind(int , int , int );
int boundarycondition(int , int , int, int );
float calculatecurvature(float, float);
float calculatelaplacian(int , int , int , int , int , int , int , int , int );
float searchphi(int, int );
void updatefornextstep();
void assembleactivegrains(int , int , int );
float checkgradient(int ,int );
void gradientvector(int, int, int, int, int, int, int, int, int, float *);
float calculategradient(int, int, int);
float findvelocity(float *, float, int);
float dG_drag(float, float, float, float, float);
float adrag(float, float, float);
float bdrag(float, float, float);
float delta_fn(float, int);
float dG_drag2(float, float, float, float);
float dG_dragfit(float, float, float, float);
float dG_drag_lookup(float);

float iw();
void printinfo_gridpoint(int , int , int );
float anisotropic_mobility(float, float, int);
float anisotropic_segregation(float, float, int);
float anisotropic_eps(float, float, int);
float anisotropic_omega(float, float, int);
float anisotropic_gamma(float, float, int);

float evaluatemisorientation(int grain1, int grain2);
float effective_velocity(float E0, float T);



//Postprocess functions
void output(float);
void printstats(int);
void write2dgnuplot();
void writedisplay(int);
void writegrainstatistics(int );
int evaluatenumberofgrains();
float evaluatevolumefraction();
float evaluateaverageradius();
float checkcurvature();
void printstatistics(int);
void printmicrostructure();

//misorientation
void populatemisorientation();
void create_gsl_matrix(float **matrix, int num, gsl_matrix *RM);
void inverse_gsl(gsl_matrix *RM, gsl_matrix *Inv);
void multiply_matrix(gsl_matrix *RM, gsl_matrix *Inv, gsl_matrix *m);
float extractmisorientation_from_misorientationmatrix(gsl_matrix *m);
void freemisorientation();

#endif
