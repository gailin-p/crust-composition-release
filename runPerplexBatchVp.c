/******************************************************************************
 * FILE: runPerplexBatchVp.c
 * COMPILATION: mpicc -std=c99 -o runPerplexBatch runPerplexBatch.c
 * USAGE: mpiexec -np N ./runPerplexBatchVp ignmajors.csv
 * DESCRIPTION:  
 *   Configures and runs PerpleX seismic velocity calculations on n processors
 *   for each bulk composition in ignmajors.csv, along a specified 
 *   geothermal gradient.
 *
 ******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <mpi.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include "arrays.h"


int main(int argc, char **argv){
	uint32_t datarows, datacolumns;
	uint32_t i, j, k;
	int world_size, world_rank, rc;


	//Check input arguments
	if (argc != 2) {
		fprintf(stderr,"USAGE: %s <input_filename>\n", argv[0]);
		exit(1);
	}

	// Start MPI
	rc = MPI_Init(&argc,&argv); 
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n"); MPI_Abort(MPI_COMM_WORLD, rc);
	}

	// Get world size (number of MPI processes) and world rank (# of this process)
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);


	if (world_rank==0){
		// Declare variables used only on the root node
		int buf[world_size-1], nextReady;
		MPI_Request reqs[world_size-1];
		MPI_Status stats[world_size-1];

		// Print format of output
//		printf("");

		// Import 2-d source data array as a flat double array. Format:
		// KV, SiO2, TiO2, Al2O3, FeO, MgO, CaO, Na2O, K2O, H2O, CO2, tc1Crust
		double** const data = csvparse(argv[1],',', &datarows, &datacolumns);

		// Listen for task requests from the worker nodes
		for (i=1; i<world_size; i++){
			//        *buf, count, datatype, dest, tag, comm, *request
			MPI_Irecv(&buf[i-1], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &reqs[i-1]);
		}

		// Once any worker asks for a new task, send next task to that worker and keep listening
		for (i=0; i<datarows; i++){
			MPI_Waitany(world_size-1, reqs, &nextReady, stats);
			//       *buf, count, datatype, dest, tag, comm
			MPI_Send(data[i], 12, MPI_DOUBLE, nextReady+1, 1, MPI_COMM_WORLD);
			//        *buf, count, datatype, source, tag, comm, *request
			MPI_Irecv(&buf[nextReady], 1, MPI_INT, nextReady+1, 0, MPI_COMM_WORLD, &reqs[nextReady]);
		}

		// Wait for all workers to complete, then send the stop signal
		MPI_Waitall(world_size-1, reqs, stats);	
		double stop[12] = {-1};
		for (i=1; i<world_size; i++){
			MPI_Send(&stop, 12, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);	
		}
	}

	else {
		// Declare variables used only on the worker nodes
		MPI_Request sReq;
		MPI_Status sStat;
		double ic[12];
		FILE *fp;
		char* prefix = malloc(500*sizeof(char));
		char* cmd_string = malloc(1000*sizeof(char));
		char* path_string = malloc(500*sizeof(char));

		// Simulation parameters
		/**********************************************************/


		// Variables that control size and location of the simulation
		/***********************************************************/	
		// Location of scratch directory (ideally local scratch for each node)
		// This location may vary on your system - contact your sysadmin if unsure
//		const char scratchdir[]="/scratch/";
		const char scratchdir[]="./"; // Just use local directory for now

		/***********************************************************/



		double dpdz = 2900. * 9.8 / 1E5 * 1E3; // Pressure gradient (bar/km)


		while (1) {
			// Ask root node for new task
			//       *buf, count, datatype, dest, tag, comm, *request
			MPI_Isend(&world_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &sReq);
			//       *buf, count, datatype, source, tag, comm, *status
			MPI_Recv(&ic, 12, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &sStat);

			// Exit loop if stop signal recieved
			if (ic[0]<0) break;


			usleep(1000);

			for (i=0; i<12; i++){
				printf("%g\t", ic[i]);
			}
			printf("\n");

			
			
//			//Set water
//			ic[9]=0.1;
//			//Set CO2
//			ic[10]=0.1;


			//Configure working directory
			sprintf(prefix,"%sout%.0f_%i/", scratchdir, world_rank, ic[0]);
			sprintf(cmd_string,"rm -rf %s; mkdir %s", prefix, prefix);
			system(cmd_string);

			// Place required data files
			sprintf(cmd_string,"cp ./*.dat %s", prefix);
			system(cmd_string);


			// Create build batch file
			sprintf(path_string, "%sbuild.txt", prefix);
			fp=fopen(path_string,"w");
			// Name, components, and basic options. Holland and Powell (2002) thermodynamic dataset and (1998) fluid equation state.
			fprintf(fp,"%.0f\nhp02ver.dat\nperplex_option.dat\nn\nn\nn\nn\nSIO2\nTIO2\nAL2O3\nFEO\nMGO\nCAO\nNA2O\nK2O\nH2O\nCO2\n\n5\n", ic[0]);
			// Pressure gradient details
			fprintf(fp,"3\nn\ny\n2\n1\n273.15\n%g\n1\n25000\ny\n", 550.0/ic[11]/dpdz);
			// Whole-rock composition
			for(i=1; i<11; i++){
				fprintf(fp,"%g ",ic[i]);
			}
			//Solution model
			fprintf(fp,"\nn\nn\ny\nsolution_model.dat\nmelt(HP)\nO(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nAnth\nChum\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar\nNeph(FB)\nScap\nDo(HP)\nF\n\nClosed System");	
			fclose(fp);

			// build PerpleX problem definition
			sprintf(cmd_string,"cd %s; build < build.txt >/dev/null", prefix);
			system(cmd_string);

			// Run PerpleX vertex calculations
			sprintf(cmd_string,"cd %s; echo %.0f | vertex > /dev/null", prefix, ic[0]);
			system(cmd_string);


			// Create werami batch file
			sprintf(path_string, "%swerami.txt", prefix);
			fp=fopen(path_string,"w");
			fprintf(fp,"%.0f\n3\n1\n1\n25000\n100\n2\nn\nn\n13\nn\nn\n15\nn\nn\n0\n0\n", ic[0]);
			fclose(fp);

			// Extract Perplex results with werami
			sprintf(cmd_string,"cd %s; werami < werami.txt", prefix);
			system(cmd_string);

			
			// If results can't be found, clean up scratch directory and move on to next simulation
			sprintf(cmd_string,"%s%.0f_1.tab", prefix, ic[0]);
			if ((fp = fopen(cmd_string, "r")) == NULL) {
				fprintf(stderr, "%s : Simulation output could not be found.\n", prefix);
				sprintf(cmd_string,"rm -r %s", prefix);
				system(cmd_string);
				continue;
			}

//			// Import results, if they exist. Format:
//			// P(bar) T(K) rho Vp(km/s) Vp/Vs
//
//			// Can delete temp files after we've read them
//			sprintf(cmd_string,"rm -r %s", prefix);
//			system(cmd_string);
//
//			// Print results. Format:
//			// Kv, P, T, Rho, Vp, Vp/Vs
//			printf("");

		}
	}
	MPI_Finalize();
	return 0;
}

