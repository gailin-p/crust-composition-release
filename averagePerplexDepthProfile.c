#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "arrays.h"

// Format of PerplexResults.tsv
#define INDEX 0
#define P 1
#define T 2
#define RHO 3
#define VP 4
#define VPVS 5

int main(int argc, char **argv){
	uint32_t datarows, datacolumns;
	uint32_t indexrows, indexcolumns;
	uint32_t i, j, k, k1, k2;

	//Check input arguments
	if (argc != 3) {
		fprintf(stderr,"USAGE: %s <PerplexResults.tsv> <index.tsv>\n", argv[0]);
		exit(1);
	}
	
	
	// Import data
	double* const data = csvparseflat(argv[1],'\t', &datarows, &datacolumns);
	double* const index = csvparseflat(argv[2],'\t', &indexrows, &indexcolumns);

	// Convert indices to unsigned intergers for faster comparision
	uint32_t* restrict target = malloc(indexrows * sizeof(uint32_t));
	uint32_t* restrict source = malloc(datarows * sizeof(uint32_t));
	for (i=0; i<indexrows; i++){
		target[i] = (uint32_t) round(index[i]);
	}
	for (j=0; j<datarows; j++){
		source[j] = (uint32_t) round(data[INDEX*datarows + j]);
	}

	
	uint32_t jmin, step;
	// Find first line of real output
	for (jmin=0; isnan(data[INDEX*datarows+jmin]); jmin++){
	}
	// Find output step size
	for (step=1; source[jmin+step]==source[jmin]; step++) {
	}
//	printf("jmin: %i\n",jmin);
//	printf("Step size: %i\n",step);



	// Make target variables and initialize with NANs
	double* restrict Rho = malloc(step * sizeof(double));
	for (i=0; i<step; i++) Rho[i]=0;
	double* restrict Vp = malloc(step * sizeof(double));
	for (i=0; i<step; i++) Vp[i]=0;
	double* restrict VpVs = malloc(step * sizeof(double));
	for (i=0; i<step; i++) VpVs[i]=0;
	uint32_t* restrict nRho = malloc(step * sizeof(double));
	for (i=0; i<step; i++) nRho[i]=0;
	uint32_t* restrict nVp = malloc(step * sizeof(double));
	for (i=0; i<step; i++) nVp[i]=0;
	uint32_t* restrict nVpVs = malloc(step * sizeof(double));
	for (i=0; i<step; i++) nVpVs[i]=0;


	
	for (i=0; i<indexrows; i++){
		for (j=jmin; j<datarows; j += step){
			if (source[j]==target[i]){
				// Add up rows
				for (k=0; k<step; k++){
					if (!isnan(data[RHO*datarows + j + k])){
						Rho[k] += data[RHO*datarows + j + k];
						nRho[k]++;
					}
					if (!isnan(data[VP*datarows + j + k])){
						Vp[k] += data[VP*datarows + j + k];
						nVp[k]++;
					}
					if (!isnan(data[VPVS*datarows + j + k])){
						VpVs[k] += data[VPVS*datarows + j + k];
						nVpVs[k]++;
					}
				}

				// And move to the next index
				break;
			} 
		}
	}




	FILE *fp;
	// Write results to file
	fp = fopen("Rho.csv","w");
	for (i=0; i<step; i++) fprintf(fp, "%g\n", Rho[i]/(double)nRho[i]);
	fclose(fp);
	fp = fopen("Vp.csv","w");
	for (i=0; i<step; i++) fprintf(fp, "%g\n", Vp[i]/(double)nVp[i]);
	fclose(fp);
	fp = fopen("VpVs.csv","w");
	for (i=0; i<step; i++) fprintf(fp, "%g\n", VpVs[i]/(double)nVpVs[i]);
	fclose(fp);

		
return 0;
}
