#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include"Global.h"
#include"functions.h"
#include "gsl/gsl_matrix_double.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"


float S[24][9] = {
											{1, 0, 0, 0, 1, 0, 0, 0, 1},
											{-1, 0, 0, 0, 1, 0, 0, 0, -1},
											{-1, 0, 0, 0, -1, 0, 0, 0, 1},
											{1, 0, 0, 0, -1, 0, 0, 0, -1},
											{0, 1, 0, 0, 0, 1, 1, 0, 0},
											{0, -1, 0, 0, 0, 1, -1, 0, 0},
											{0, -1, 0, 0, 0, -1, 1, 0, 0},
											{0, 1, 0, 0, 0, -1, -1, 0, 0},
											{0, 0, 1, 1, 0, 0, 0, 1, 0},
											{0, 0, -1, 1, 0, 0, 0, -1, 0},
											{0, 0, -1, -1, 0, 0, 0, 1, 0},
											{0 , 0, 1, -1, 0, 0, 0, -1, 0},
											{0, 0, -1, 0, -1, 0, -1, 0, 0},
											{0, 0, 1, 0, -1, 0, 1, 0, 0},
											{0, 0, 1, 0, 1, 0, -1, 0, 0},
											{0, 0, -1, 0, 1, 0, 1, 0, 0},
											{-1, 0, 0, 0, 0, -1, 0, -1, 0},
											{1, 0, 0, 0, 0, -1, 0, 1, 0},
											{1, 0, 0, 0, 0, 1, 0, -1, 0},
											{-1, 0, 0, 0, 0, 1, 0, 1, 0},
											{0, -1, 0, -1, 0, 0, 0, 0, -1},
											{0, 1, 0, -1, 0, 0, 0, 0, 1},
											{0, 1, 0, 1, 0, 0, 0, 0, -1},
											{0, -1, 0, 1, 0, 0, 0, 0, 1},
								};


void populatemisorientation()
{

	printf("!!!POPULATING MISORIENTATION!!!!\n");
	extern float **fullmisorientation;
	extern float **rotationmatrix;
	extern int ngrainsinit, ngrainsor;
	extern float **symmetryoperator;
	extern int symm_length;

	symm_length=24;


	symmetryoperator=(float **)malloc((symm_length) * sizeof(float *)); // 9 rotation matrix 00, 01, 02, 10, 11, 12, 20, 21, 22
    for(int i=0; i<symm_length; i++) symmetryoperator[i] = (float *)malloc(9 * sizeof(float));

    //for(int i=0; i<symm_length; i++) for(int j=0; j<9; j++) printf("%f \n", S[i][j]);

    for(int i=0; i<symm_length; i++) for(int j=0; j<9; j++) symmetryoperator[i][j]=S[i][j];


    printf("Read Symmetry operators : Cubic Symmetry\n");

	int grain1, grain2;
	int size=3;
	float min_angle, angle;
	
	//fullmisorientation=(float **)malloc((ngrainsinit) * sizeof(float *)); // 9 rotation matrix 00, 01, 02, 10, 11, 12, 20, 21, 22
  //  for(int i=0; i<ngrainsinit; i++) fullmisorientation[i] = (float *)malloc(ngrainsinit* sizeof(float));
	fullmisorientation=(float **)malloc((ngrainsor) * sizeof(float *)); // 9 rotation matrix 00, 01, 02, 10, 11, 12, 20, 21, 22
  for(int i=0; i<ngrainsor; i++) fullmisorientation[i] = (float *)malloc(ngrainsor* sizeof(float));

    //initialize
  //for(int i=0; i<ngrainsinit;i++) for(int j=0; j<ngrainsinit; j++) fullmisorientation[i][j]=0.0;
 	for(int i=0; i<ngrainsor;i++) for(int j=0; j<ngrainsor; j++) fullmisorientation[i][j]=0.0;

    gsl_matrix *temp_matrix=gsl_matrix_alloc(size, size);
	gsl_matrix *Symm_matrix = gsl_matrix_alloc(size, size);
	gsl_matrix *RM1 = gsl_matrix_alloc(size, size);
	gsl_matrix *RM2 = gsl_matrix_alloc(size, size);
	gsl_matrix *Inv_RM1 = gsl_matrix_alloc(size, size);
	gsl_matrix *misorientationmatrix = gsl_matrix_alloc(size, size);
		
		

    //for(int i=0; i<ngrainsinit; i++)
		for(int i=0; i<ngrainsor; i++)
    {
      for(int j=i+1; j<ngrainsor; j++)
      {

      	min_angle=360;
      	grain1=i; grain2=j;
				create_gsl_matrix(rotationmatrix, grain1, RM1);
				create_gsl_matrix(rotationmatrix, grain2, RM2);
				inverse_gsl(RM1, Inv_RM1);
				for(int k=0; k<symm_length; k++)
				{
					create_gsl_matrix(symmetryoperator, k, Symm_matrix);
					multiply_matrix(RM2, Symm_matrix, temp_matrix);
					multiply_matrix(temp_matrix, Inv_RM1, misorientationmatrix);
					angle = extractmisorientation_from_misorientationmatrix(misorientationmatrix);
					//printf("Angle : %f\n", angle);
					if(angle<min_angle) min_angle=angle;
				}
				fullmisorientation[grain1%ngrainsor][grain2%ngrainsor]=min_angle;
				fullmisorientation[grain2%ngrainsor][grain1%ngrainsor]=min_angle;
      }
    }
    gsl_matrix_free(RM1);
    gsl_matrix_free(RM2);
	gsl_matrix_free(Inv_RM1);
	gsl_matrix_free(misorientationmatrix);
	gsl_matrix_free(Symm_matrix);
	gsl_matrix_free(temp_matrix);
}


void create_gsl_matrix(float **matrix, int num, gsl_matrix *RM)
{
	gsl_matrix_set(RM, 0, 0, matrix[num][0]);
	gsl_matrix_set(RM, 0, 1, matrix[num][1]);
	gsl_matrix_set(RM, 0, 2, matrix[num][2]);
	gsl_matrix_set(RM, 1, 0, matrix[num][3]);
	gsl_matrix_set(RM, 1, 1, matrix[num][4]);
	gsl_matrix_set(RM, 1, 2, matrix[num][5]);
	gsl_matrix_set(RM, 2, 0, matrix[num][6]);
	gsl_matrix_set(RM, 2, 1, matrix[num][7]);
	gsl_matrix_set(RM, 2, 2, matrix[num][8]);
}

void inverse_gsl(gsl_matrix *RM, gsl_matrix *Inv)
{
	int size=3;
	gsl_permutation *p = gsl_permutation_alloc(size);
	int s;
	gsl_linalg_LU_decomp(RM, p, &s);
	gsl_linalg_LU_invert(RM, p, Inv);
	gsl_permutation_free(p);
}

void multiply_matrix(gsl_matrix *RM, gsl_matrix *Inv, gsl_matrix *m)
{
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, RM, Inv, 0.0, m);
}

float extractmisorientation_from_misorientationmatrix(gsl_matrix *m)
{
	float g11, g22, g33;
	g11 = gsl_matrix_get(m, 0, 0);
	g22 = gsl_matrix_get(m, 1, 1);
	g33 = gsl_matrix_get(m, 2, 2);
	return acos((g11+g22+g33-1)/2)/RADIAN;
}


void freemisorientation()
{
	extern int ngrainsinit;
	extern int ngrainsmax;
	extern float **orientation, **rotationmatrix, **fullmisorientation;
	for(int i=0; i<ngrainsinit; i++) free(orientation[i]);
	//for(int i=0; i<ngrainsinit; i++) free(rotationmatrix[i]);
	//for(int i=0; i<ngrainsinit; i++) free(fullmisorientation[i]);	
	for(int i=0; i<ngrainsor; i++) free(rotationmatrix[i]);
	for(int i=0; i<ngrainsor; i++) free(fullmisorientation[i]);	
}