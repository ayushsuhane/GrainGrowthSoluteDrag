#include<stdio.h>
#include<omp.h>
#include"Global.h"
#include"functions.h"
#include"stdio.h"
#include"string.h"


void main()
{
	extern float dt, c_time, real_totaltime, sim_time;
	extern int t_output, t_print, t_createcheckpoint; // will create checkpoint after every t_createcheckpoint steps

	float runtime = omp_get_wtime();

	sim_time=0;
	int i=0;
	printf("NOTE: The Simulation will overwrite the restart file everytime, make sure to re-copy the file before starting the simulations");
	
	simulationsetup();
	while(sim_time<real_totaltime)
	{
	  if(i%(t_print)==0){printstats(i);printstatistics(i);/*printmicrostructure();*/} 
	  if(i%(t_output)==0)
	  {
	  	writedisplay(i);
	  	writegrainstatistics(i);
	  	
	  } 
	  if(i%(t_createcheckpoint)==0) output(sim_time); // sim_time: simulation time, 0: create restartfile with time =0
	  // This will overwrite the restart file, make sure to copy the correct restart file every time you run the simulations

	  pfsolver(i);
//	  Postprocessdata();
	  i+=1;
	  sim_time += dt*c_time;
	}
//	Printstats();
	output(0);
	write2dgnuplot();
	freememory();
	runtime = omp_get_wtime() - runtime;
	extern int NUM_THREADS;
	printf("Total Runtime with %d threads : %lf\n", NUM_THREADS, runtime);
	return;
}
