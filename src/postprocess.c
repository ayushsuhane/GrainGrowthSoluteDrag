#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"Global.h"
#include"functions.h"

void output(float local_time)
{
	extern float *grid, *c;
	extern int *flaggrid;
	extern int nx, ny, nz, dim, grid2, grid3, ngrainsmax;
	extern int writetorestart;

	int index, count, gind;
	FILE *fo;
	char *outputfile;

	outputfile = "output/restartfile.txt";
	
	fo = fopen(outputfile, "w");

	printf("Inside outputfunction\n");
	
	fprintf(fo, "%e\n", local_time);
	for(int i=0; i<nx; i++)
	{
	  for(int j=0; j<ny; j++)
	  {
	    for(int k=0; k<nz; k++)
		  {
		  	index = k*grid2 + j*nx + i;
		  	count = 0; 
		  	for(int l=0; l<ngrainsmax; l++) {gind=index+l*grid3; if(flaggrid[gind]==-1) break; count+=1;}
		  	fprintf(fo, "%d %d %d %d", i, j, k, count);
		  	for(int l=0;l<ngrainsmax;l++)
		  	{
		  	  if(flaggrid[l*grid3 + index] >= 0) fprintf(fo, " %d %e", flaggrid[l*grid3 + index], grid[l*grid3 + index]); 
		  	}
		  	fprintf(fo, "\n");
		  }
	  }
	}
	fclose(fo);
	//free(outputfile);
	return;
}

void printstats(int t)
{
	extern float dt, c_time;
	extern int trackid, ngrainsmax;
	extern int *flaggrid;
	extern float sim_time;
	extern int nx, ny, nz, grid2, grid3;
	extern int dim;
	
	//printf("\n");
	printf("Timestep: %d, Time : %e ", t, sim_time);
	if(dim==1) printf("volumefraction : %lf \n", evaluatevolumefraction());
	else printf("volumefraction : %lf Radius : %lf, average radius : %lf\n", evaluatevolumefraction(), sqrt(evaluatevolumefraction()*grid3/PI),  evaluateaverageradius());
	/*
	extern float *curvature;
	extern float threshold;
	int gind;
	for(int i=0; i<nx; i++)
	{
	  for(int j=0; j<ny; j++)
	  {
	  	for(int k=0; k<nz; k++)
	  	{
	  	  for(int l=0; l<ngrainsmax; l++)
	  	  {
	  	  	gind = l*grid3 + k*grid2 + j*nx + i;
	  	  	if(flaggrid[gind]==-1) continue;
	  	  	if(curvature[gind]>threshold) printf("i: %d, j : %d, k: %d, l : %d, curvature : %lf\n", i, j, k, l, curvature[gind]);
	  	  }	
	  	}
	  }
	  //printgrainstatistics();
	}
	*/
	int i=trackid;
	//for(int l=0; l<ngrainsmax; l++) printf("flaggrid: %d ", flaggrid[l*grid3 + i]);
    //printf("\n");
	
}

void write2dgnuplot()
{
	extern float *grid, *c;
	extern int *flaggrid;
	extern int nx, ny, nz, dim, grid2, grid3, ngrainsmax;

	int index,  id=1, l=0;
	printf("Printing file for gnuplot\n");
	FILE *fo;
	fo = fopen("output/gnuplot.txt", "w");

	for(int i=0; i<nx; i++)
	{
	  for(int j=0; j<ny; j++)
	  {
	  	//printf("%d %d\n",i, j);
	    index =  j*nx + i;
	    id=1;
		fprintf(fo, "%d %d", i, j);
		for(l=0;l<ngrainsmax;l++)
		{
		  if(flaggrid[l*grid3 + index] == 1) {id=0;break;} 
		}
		if(id==0) fprintf(fo, " %e", grid[l*grid3 + index]);
		else fprintf(fo, " 0.0");
		fprintf(fo, "\n");
	  }
	  fprintf(fo, "\n");
	}
	fclose(fo);
}

void writedisplay(int t)
{
	extern int nx, ny, nz, grid2, grid3, ngrainsmax;
	extern int flag_curvature, bc;
	extern int flag_anisotropy, flag_e0anisotropy, flag_gbenanisotropy;
	extern float averagecurvature;
	extern int t_output;
	extern float threshold;
	extern float sim_time;
	float misvariable=0;
	averagecurvature=0;

	char FILENAME[100];
	 //sprintf(FILENAME, "output/visual.txt");
	int check = sprintf(FILENAME, "output/visual_%0.2e.txt", sim_time);
	//printf("Check : %d\n", check);
	FILE *fd; 
	fd = fopen(FILENAME, "w");
	int i, j, k, l, index, gind, id;
	float sum=0, sumv=0;
	float mag_kappa=0;
	int count = 0, countgrain=0;
	float misorientation;
	int grain1=0, grain2=0, lead_grain=0;
	float sumg = 0;
	float lead_grain_phi=0;
	for(i=0; i<nx; i++)
	{
	  for(j=0; j<ny; j++)
	  {
	  	for(k=0; k<nz; k++)
	  	{
	  		misorientation = 0;
	  		countgrain=0;
	  		mag_kappa=0;
	  		sum=0;
	  		sumv=0;
	  		sumg = 0;
	  		grain1=0;
	  		grain2=0;
	  		index = k*grid2 + j*nx + i;
	  		lead_grain_phi = 0;
	  		lead_grain = 0;
	  		for(l=0; l<ngrainsmax; l++)
	  		{
	  			gind = index + l*grid3;
	  			if(grid[gind]>lead_grain_phi) {lead_grain_phi = grid[gind]; lead_grain = flaggrid[gind];}
	  			if(flaggrid[gind]==-1) break;
	  			countgrain += 1;
	  			sum += grid[gind]*grid[gind];
	  			sumv += (flaggrid[gind]+1)*grid[gind];
	  			sumg += (flaggrid[gind]+1)*(flaggrid[gind]+1);
	  		}
	  		if((flag_anisotropy||flag_e0anisotropy||flag_gbeanisotropy) && countgrain>1)
	  		{
	  		  	/***Misorientation calculations *****/  
	  	  	  int maxindex=0, maxindex1=1;
	          float maxphi=0.0, maxphi1=0.0;
	          for(l=0; l<ngrainsmax;l++)
	          {
		        if(grid[index+l*grid3]>maxphi) {maxindex=index+l*grid3; maxphi=grid[maxindex]; grain1 = flaggrid[maxindex];}
	          }
	          for(l=0; l<ngrainsmax;l++)
	          {
		        if(grid[index+l*grid3]>maxphi1 && ((index+l*grid3)!=maxindex)) {maxindex1=index+l*grid3; maxphi1=grid[maxindex1]; grain2=flaggrid[maxindex1];}
	          }
	          //if(grid[maxindex1]>0.9) misorientation=0;//Can uncomment
	          if(sum>0.9) misorientation=0;//Here, classifying based on square of phi values. above phi=0.95 and phi2=0.05, it assumes a bulk
	          
	          else misorientation = evaluatemisorientation(grain1, grain2);

	          ////misorientation = fabs(orientation[grain1] - orientation[grain2]);
	          //printf("i : %d, j: %d, grain1 : %d, grain2 :%d, misorientation : %lf\n", i, j, grain1, grain2, misorientation);
	  	      /********/
	  		}
	  		if((flag_anisotropy||flag_e0anisotropy||flag_gbeanisotropy))
	  		{
	  			if(countgrain>1){misvariable = sqrt(orientation[grain1%ngrainsor][0]*orientation[grain1%ngrainsor][0] + orientation[grain1%ngrainsor][1]*orientation[grain1%ngrainsor][1] + orientation[grain1%ngrainsor][2]*orientation[grain1%ngrainsor][2]);}
	  			else{misvariable = sqrt(orientation[flaggrid[index]%ngrainsor][0]*orientation[flaggrid[index]%ngrainsor][0] + orientation[flaggrid[index]%ngrainsor][1]*orientation[flaggrid[index]%ngrainsor][1] + orientation[flaggrid[index]%ngrainsor][2]*orientation[flaggrid[index]%ngrainsor][2]);}
	  		}
	  		fprintf(fd, "%d %d %d %lf %lf %lf %lf %d ", i, j, k, sum, sumv, sumg, misvariable, lead_grain);
	  		if(flag_curvature) {mag_kappa = checkcurvature(index); fprintf(fd, "%le ", mag_kappa);  if(fabs(mag_kappa)>threshold) {averagecurvature+=mag_kappa; count+=1;/*printf("i :%d, j: %d, k:%d, Kappa : %le\n", i, j, k, mag_kappa);*/}}
	  		if(flag_anisotropy||flag_e0anisotropy||flag_gbeanisotropy) fprintf(fd, "%lf %d %d", misorientation, grain1, grain2);
	  		if(bc==0){id = boundarycondition(i, j, k, dim);if(!id) {fprintf(fd, "\n"); continue;}}
	  		fprintf(fd, "\n");
	  	}
	  }
	  if(dim!=1) fprintf(fd, "\n");
	}
	fclose(fd);
	averagecurvature = averagecurvature/count;			
	return;
}

float evaluatevolumefraction()
{
	extern int *flaggrid;
	extern float *grid;
	extern int nx, ny, nz, grid2, grid3;
	float vf=0;
	int index, gind;

	for(int i=0; i<nx; i++)
	{
	  for(int j=0; j<ny; j++)
	  {
	  	for(int k=0; k<nz; k++)
	  	{
	  		index = k*grid2 + j*nx + i;
	  		for(int l=0; l<ngrainsmax; l++)
	  		{
	  			gind = index + l*grid3;
	  			if(flaggrid[gind]==-1) break;
	  			if(flaggrid[gind]==0) vf += grid[gind]/(grid3);
	  		}

	  	}
	  }
	}
	return vf;
}

int evaluatenumberofgrains()
{
	extern float *grid;
	extern float threshold, p_threshold;
	extern int *flaggrid;
	extern int nx, ny, nz, grid2, grid3, ngrainsmax;
	extern int ngrainsinit;
	int i, j, k, l, a, index, gind, count=0, id;
	int grains[ngrainsinit];
	//printf("Calculating Number of grains!");
	for(i=0; i<nx; i++)
	{
	  for(j=0; j<ny; j++)
	  {
	  	for(k=0; k<nz;k++)
	  	{
	  	  index = k*grid2 + j*nx + i;
	  	  for(l=0; l<ngrainsmax; l++)
	  	  {
	  	  	id=1;
	  	  	gind = index + l*grid3;
	  	  	if(flaggrid[gind]==-1) break;
	  	  	//if(grid[gind]<p_threshold) continue;//Can uncomment 
	  	  	if(grid[gind]<1e-3) continue;
	  		for(a=0; a<count; a++) {if(flaggrid[gind]==grains[a]) {id=0; break;}}
	  		if(id) {grains[count]=flaggrid[gind]; count+=1;}
	  	  }
	  	}

	  }
	}
	return count;
}

float evaluateaverageradius()
{
	/* I think its correct now*/
	/*Wrong, hve to correct*/
	extern int dim;
	int count = evaluatenumberofgrains();
	if(dim==1) return 0.0;
	if(dim==2) return sqrt(grid3/(count*PI));
	if(dim==3) return pow(3*grid3/(count*4*PI), (1/3.0));
	return 0;
}


void writegrainstatistics(int t)
{
	extern float *grid;
	extern int *flaggrid;
	extern int numbergrains;
	extern int nx, ny, nz, grid2, grid3, ngrainsmax, dim;
	extern float sim_time;
	

	int i, j, k, l, a, index, gind, count, id;
	//int grainsdata[numbergrains]; //need to know the number of grains, so have to call evaluatenumberofgrains before this function
	int grainsdata[20000];
	//int *grainsdata;
	//grainsdata=malloc(numbergrains *sizeof(int));


	printf("Initial number of grains : %d\n", numbergrains);
	
	for(i=0; i<20000; i++) grainsdata[i]=0;

	for(i=0; i<nx; i++)
	{
	  for(j=0; j<ny; j++)
	  {
	  	for(k=0; k<nz; k++)
	  	{
	  		index=k*grid2 + j*nx + i;
	  		for(l=0; l<ngrainsmax; l++)
	  		{
	  			gind = index + l*grid3;
	  			if(flaggrid[gind]==-1) break;
	  			if(grid[gind]>0.5) grainsdata[flaggrid[gind]] += 1;  
	  		}
	  	}
	  }
	}

	printf("Ngrains remaining : %d\n", evaluatenumberofgrains());
	char FILENAME[200];
	int check = sprintf(FILENAME, "output/gstat_%0.2e.txt", sim_time);
	
	FILE *fg = fopen(FILENAME, "w");
	if (fg == NULL) {printf("Cannot open file\n"); perror("Failed: ");return;}
	//printf("fg address : %d\n", &fg);

	//printf("Check : %d\n", check);
	float maxvol_indgrain = 0;
	float avgrad_indgrain = 0;
	for(i=0; i<20000; i++) 
	{
	  if(grainsdata[i]==0) continue;
	  fprintf(fg, "%d ",i);
	  if(grainsdata[i]>maxvol_indgrain) maxvol_indgrain=grainsdata[i];
	  if(dim==2) avgrad_indgrain=sqrt(grainsdata[i]/(PI));
	  if(dim==3) avgrad_indgrain=pow((3*grainsdata[i]/(4*PI)), (1/3.0));
	  fprintf(fg, "%d %le\n", grainsdata[i], avgrad_indgrain);
	  //printf("After statistics writing\n");
	  //printf("After statistics writing2\n");
	  
	}
	fclose(fg);
	//free(grainsdata);
	//for(i=0; i<numbergrains; i++) {if(grainsdata[i]!=0) fprintf(fgsd, "%d : %lf\n", i, grainsdata[i]);}
	return;

}

void printstatistics(int t)
{
	extern FILE *fs;
	extern float *grid;
	extern int *flaggrid;
	extern int dim;
	extern int numbergrains;
	extern float averagecurvature;
	extern int flag_curvature;
	extern float expected_velocity;
	extern float sim_time;
	extern int flag_geometry;
	extern float c_length;
	extern int grid3;
	extern float aver_radius;
	extern int readfromrestart; // if reading from restart, append to the statistics file if sim time is not 0
	extern float maxratio; //describes the ratio of maximum grain with the average grain radius
	int i, j, k, l, gind, index;
	int grainsdata[20000];
	for(i=0; i<numbergrains; i++) grainsdata[i]=0;
	for(i=0; i<nx; i++)
	{
	  for(j=0; j<ny; j++)
	  {
	  	for(k=0; k<nz; k++)
	  	{
	  		index=k*grid2 + j*nx + i;
	  		for(l=0; l<ngrainsmax; l++)
	  		{
	  			gind = index + l*grid3;
	  			if(flaggrid[gind]==-1) break;
	  			//if(grid[gind]>0.5) grainsdata[flaggrid[gind]] += 1;  // Uncomment it for sharp interface calculations
	  			grainsdata[flaggrid[gind]] += grid[gind];
	  		}
	  	}
	  }
	}
	float maxvol_indgrain = 0;
	float avgrad_indgrain = 0, maxrad_indgrain=0;
	for(i=0; i<numbergrains; i++) 
	{
	  if(grainsdata[i]==0) continue;
	  if(dim==2) avgrad_indgrain=sqrt(grainsdata[i]/(PI));
	  if(dim==3) avgrad_indgrain=pow((3*grainsdata[i]/(4*PI)), (1/3.0));
	  if(grainsdata[i]>maxvol_indgrain){ maxvol_indgrain=grainsdata[i]; maxrad_indgrain = avgrad_indgrain;}
	}


	if(flag_geometry!=2) aver_radius =  evaluateaverageradius()*dx*c_length; //circular geometry
	else aver_radius = sqrt(evaluatevolumefraction()*grid3/PI);
	maxratio = maxrad_indgrain/(aver_radius/(dx*c_length)); //Only works for multiple grains

	
	if(flag_geometry!=2){
		fprintf(fs, "%e %e %e ",sim_time, evaluateaverageradius(), maxratio);
	}
	else{
		fprintf(fs, "%e %e %e ",sim_time, aver_radius, maxratio);
	}
	fflush(fs);

	fprintf(fs, "%e ", evaluatevolumefraction());
	fflush(fs);
	if(flag_curvature) fprintf(fs, "%e ", averagecurvature);
	fprintf(fs, "\n");
	printf("Ngrains remaining : %d\n", evaluatenumberofgrains());
	fflush(fs);
	//printf("Statistics written\n");
	return;
}

float checkcurvature(int index)
{
	extern int ngrainsmax;
	extern int *flaggrid;
	extern float *curvature;
	int l, count=0;
	float kappa1, kappa2;
	// Check if the magnitude of curvature is within 10%

	for(l=0; l<ngrainsmax;l++) {if(flaggrid[index + l*grid3]==-1) break; count+=1;}
	if(count!=2) return 0.0;

	int maxindex=0, maxindex1=1;
	float maxphi=0.0, maxphi1=0.0;
	for(l=0; l<ngrainsmax;l++)
	{
		if(grid[index+l*grid3]>maxphi) {maxindex=index+l*grid3; maxphi=grid[maxindex];}
	}
	for(l=0; l<ngrainsmax;l++)
	{
		if(grid[index+l*grid3]>maxphi1 && ((index+l*grid3)!=maxindex)) {maxindex1=index+l*grid3; maxphi1=grid[maxindex1];}
	}


	kappa1 = curvature[maxindex];
	kappa2 = curvature[maxindex1];

	//else return fabs(kappa1);
	if(((fabs(kappa1)-fabs(kappa2))/fabs(kappa1) < 0.01) && !(grid[maxindex]<0.15 || grid[maxindex]>0.85)) {/*printf("index : %d, k1 : %lf, k2 :%lf\n", index, kappa1, kappa2)*/; return fabs(kappa1);}
	else return 0.0;
}

void printmicrostructure()
{
	extern float *grid, *c;
	extern int *flaggrid;
	extern int nx, ny, nz, dim, grid2, grid3, ngrainsmax;
	extern int writetorestart;

	int index, count, gind;
	FILE *fo2;
	char *outputfile2;
	
//	int micr[ny][nx];
	int maxgrain, maxindex;
	float maxphi;
//	for(int i=0;i<nx; i++){for(int j=0;j<ny;j++) micr[j][i] = 0;}

	outputfile2 = "output/microstructure.txt";
	
	fo2 = fopen(outputfile2, "w");

	printf("Inside outputfunction\n");
	
	for(int j=0; j<ny; j++)
	{
	  for(int i=0; i<nx; i++)
	  {
	  	maxgrain = 0; maxphi=0.0; maxindex=0;
		index = j*nx + i;
		for(int l=0; l<ngrainsmax;l++){ if(grid[index+l*grid3]>maxphi) {maxindex=index+l*grid3; maxphi=grid[maxindex];}}
		fprintf(fo2, "%d ", flaggrid[maxindex]);
	  }
	  fprintf(fo2, "\n");
	}

	//fclose(fo2);
	//free(outputfile);
	return;

}
