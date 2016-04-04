#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "arrays.h"
//#include "pcg_variants.h"
//#include "gauss.h"

int main(int argc, char **argv){
	uint32_t datarows, datacolumns;
	uint32_t indexrows, indexcolumns;
	uint32_t prows, pcolumns;
	uint32_t i, j, k, kmin;

	//Check input arguments
	if (argc != 4) {
		fprintf(stderr,"USAGE: %s <PerplexResults.tsv> <index.tsv> <pressure.tsv>\n", argv[0]);
		exit(1);
	}
	
	
	// Import data
	double* const data = csvparseflat(argv[1],'\t', &datarows, &datacolumns);
	double* const index = csvparseflat(argv[2],'\t', &indexrows, &indexcolumns);
	double* const pressure = csvparseflat(argv[3],'\t', &prows, &pcolumns);
//	printf("Data rows: %i, Data columns: %i\n",datarows,datacolumns);
//	printf("Index rows: %i, Index columns: %i\n",indexrows,indexcolumns);
//	printf("Pressure rows: %i, Pressure columns: %i\n",prows,pcolumns);


	// Convert indices to unsigned intergers for faster comparision
	uint32_t* restrict target = malloc(indexrows * sizeof(uint32_t));
	uint32_t* restrict source = malloc(datarows * sizeof(uint32_t));
	for (i=0; i<indexrows; i++){
		target[i] = (uint32_t) round(index[i]);
	}
	for (j=0; j<datarows; j++){
		source[j] = (uint32_t) round(data[0*datarows + j]);
	}

	
	uint32_t jmin, step;
	// Find first line of real output
	for (jmin=0; isnan(data[0*datarows+jmin]); jmin++){
	}
	// Find output step size
	for (step=1; source[jmin+step]==source[jmin]; step++) {
	}
//	printf("jmin: %i\n",jmin);
//	printf("Step size: %i\n",step);



	// Make target variables and initialize with NANs
	double* restrict Rho = malloc(indexrows * sizeof(double));
	for (i=0; i<indexrows; i++) Rho[i]=NAN;
	double* restrict Vp = malloc(indexrows * sizeof(double));
	for (i=0; i<indexrows; i++) Vp[i]=NAN;
	double* restrict VpVs = malloc(indexrows * sizeof(double));
	for (i=0; i<indexrows; i++) VpVs[i]=NAN;
	uint32_t* restrict kmins = malloc(indexrows * sizeof(uint32_t));


	double r, rmin, p;

//	pcg32_random_t rng;
//	pcg32_srandom_r(&rng,time(NULL), clock());
	
	for (i=0; i<indexrows; i++){
		for (j=jmin; j<datarows; j += step){
			if (source[j]==target[i]){
				// Find best pressure match
				p = pressure[i];
				rmin = p;
				kmin = 0;
				for (k=0; k<100; k++){
					r=fabs(p - data[1*datarows + j + k]);
					if (r < rmin) {
						rmin = r;
						kmin = k;
					}
				}

				// Save fitting results
				Rho[i] =  data[3*datarows + j + kmin];// * (1 + pcg_gaussian_ziggurat(&rng, 0.005)); //Rho
				Vp[i] = data[4*datarows + j + kmin];// * (1 + pcg_gaussian_ziggurat(&rng, 0.005)); //Vp
				VpVs[i] = data[5*datarows + j + kmin];// * (1 + pcg_gaussian_ziggurat(&rng, 0.005)); //VpVs

				// And move to the next index
				break;
			} 
		}
	}

	FILE *fp;
	// Write results to file
	fp = fopen("Rho.csv","w");
	for (i=0; i<indexrows; i++) fprintf(fp, "%g\n", Rho[i]);
	fclose(fp);
	fp = fopen("Vp.csv","w");
	for (i=0; i<indexrows; i++) fprintf(fp, "%g\n", Vp[i]);
	fclose(fp);
	fp = fopen("VpVs.csv","w");
	for (i=0; i<indexrows; i++) fprintf(fp, "%g\n", VpVs[i]);
	fclose(fp);

		
return 0;
}
