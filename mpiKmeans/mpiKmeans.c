#include "mpi.h"
#include <stdio.h>
#include <string.h>

#include "Util.h"

#define MASTER 0
#define DNA_MIN 1
#define TD_MIN 0.1
#define MAX_DIST 2147483647 


char bases[4] = {'A','C','G','T'};

int main(int argc,char** argv)
{
	/* number of task, task rank */
	int numtasks, rank, sendcnt, recvcnt;

	/* point dimension, number of points*/
	int dimension, lines, cnt;
	
	/* file to be read */
	FILE* fp;

	char* filename;
	
	/* content read from file */
	char* dnaSource;
	double* tdSource;

	/* data points for every process */
	char* dnaBuf;
	double* tdBuf;

	/* centroid */
	char* dnaC;
	double* tdC;


	/* labels for each point */
	int* labels;

	int dna;
	int cluster;

	if (argc < 5) {
		printf("Usage: mpiKmeans <dna or 2d> <input file name> <lines> <clusters> <dimension>\n");
		exit(1);
	}

	if (!strcmp(argv[1],"dna") && argc < 6) {
		printf("Usage: mpiKmeans <dna or 2d> <input file name> <lines> <clusters> <dimension>\n");
		exit(1);
	}
	
	filename = argv[2];
	lines = atoi(argv[3]);
	dna = !strcmp(argv[1], "dna");
	cluster = atoi(argv[4]);


	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	if (dna) {
		dimension = atoi(argv[5]);
		sendcnt = (lines + numtasks - 1)/numtasks * dimension;
		recvcnt = sendcnt;
		dnaSource = malloc(sizeof(char) * dimension * lines);
		dnaBuf = malloc(sizeof(char) * recvcnt);
		dnaC = malloc(sizeof(char) * cluster * dimension);
	} else {
		dimension = 2;
		sendcnt = (lines + numtasks - 1)/numtasks * dimension;
		recvcnt = sendcnt;
		tdSource = malloc(sizeof(double) * dimension * lines);
		tdBuf = malloc(sizeof(double) * recvcnt);
		tdC = malloc(sizeof(double) * cluster * dimension);
	}

	
	if (rank == MASTER) {
		fp = fopen(filename, "r");	
		if(fp == NULL) {
			printf("Cannot open file %s\n", filename);
			exit(1);
		}
		
		if (dna) {
			readDNA(fp, dnaSource, dimension);
			generateDNACenter(dnaC, dnaSource, lines, dimension, cluster);
		} else {
			read2D(fp, tdSource);
			generate2DCenter(tdC, tdSource, lines, dimension, cluster);
		}	
		fclose(fp);
		
	}

	/* compute the number and offset of data elements to send */
	int* sendcnts = malloc(sizeof(int)*numtasks);
	int* strides = malloc(sizeof(int)*numtasks);
	int i;
	for(i = 0;i<numtasks-1;i++) {
		sendcnts[i] = sendcnt;
	}
	sendcnts[numtasks-1] = lines * dimension - sendcnt * (numtasks - 1);
	int base = 0;
	for(i = 0;i<numtasks;i++) {
		strides[i] = base;
		base += sendcnts[i];
	}
	/* scatter the data and broadcast the center */
	if(dna) {
		MPI_Scatterv (dnaSource,sendcnts,strides,MPI_CHAR,dnaBuf,recvcnt,MPI_CHAR,MASTER,MPI_COMM_WORLD);
		MPI_Bcast (dnaC,cluster * dimension,MPI_CHAR,MASTER,MPI_COMM_WORLD); 
	} else {
		MPI_Scatterv (tdSource,sendcnts,strides,MPI_DOUBLE,tdBuf,recvcnt,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
		MPI_Bcast (tdC,cluster * dimension,MPI_DOUBLE,MASTER,MPI_COMM_WORLD); 
	}

	int num = sendcnts[rank]/dimension;
	labels = malloc(sizeof(int) * num);
	
	do {
		int j;
		if(dna) {
			int* newcnts = malloc(sizeof(int) * cluster * dimension * 4);
			int* recvcnts = malloc(sizeof(int) * cluster * dimension * 4);
			memset(newcnts, 0, sizeof(int) * cluster * dimension * 4);
			char* recvctrs = malloc(sizeof(char) * cluster * dimension);
			for(i = 0; i < num; i++) {
				int lb = -1;
				int dDNA = MAX_DIST;
				for(j = 0; j < cluster;j++) {
					int tmpDist = distDNA(dnaC+j*dimension, dnaBuf+i*dimension, dimension);		
					//printf("tmpDist%d:%d\n",j,tmpDist);
					if (tmpDist < dDNA) {
						dDNA = tmpDist;
						lb = j;
					}
				}	
				for(j = 0; j < dimension; j++) {
					switch(dnaBuf[i*dimension+j]) {
						case 'A':
							newcnts[4 * dimension * lb + 4 * j + 0]++;
							break;
						case 'C':
							newcnts[4 * dimension * lb + 4 * j + 1]++;
							break;
						case 'G':
							newcnts[4 * dimension * lb + 4 * j + 2]++;
							break;
						case 'T':
							newcnts[4 * dimension * lb + 4 * j + 3]++;
							break;
						default:
							break;
					}
				}
				labels[i] = lb;
			}
			
			MPI_Reduce(newcnts, recvcnts, cluster * 4 * dimension, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
			int flag = 0;
			if(rank == MASTER) {
				for (i = 0; i < cluster; i++) {
					for (j = 0; j < dimension; j++) {
						int left = 
						recvcnts[i*dimension*4+4*j+0] > recvcnts[i*dimension*4+4*j+1]?0:1;
						int right = 
						recvcnts[i*dimension*4+4*j+2] > recvcnts[i*dimension*4+4*j+3]?2:3;

						int final = recvcnts[i*dimension*4+4*j+left] > 
						recvcnts[i*dimension*4+4*j+right]?left:right;

						recvctrs[i * dimension + j] = bases[final];
					}
				}
				int dsum = 0;
				for(i =0;i<cluster;i++) {
					dsum += distDNA(recvctrs+dimension*i,dnaC+dimension*i,dimension);
				}
				if(dsum <= DNA_MIN)
					flag = 1;
				free(dnaC);
				dnaC = recvctrs;
				/*for(i=0;i<cluster;i++) {
					for(j=0; j<dimension;j++) {
						printf("%c,",dnaC[dimension*i+j]);	
					}
					printf("\n");
				}*/
			}
			MPI_Bcast (&flag,1,MPI_INT,MASTER,MPI_COMM_WORLD); 
			free(newcnts);
			free(recvcnts);
			if(flag)
				break;
			
			MPI_Bcast (dnaC,cluster * dimension,MPI_CHAR,MASTER,MPI_COMM_WORLD); 
			
		} else {
			double* newctrs = malloc(sizeof(double) * cluster * dimension);
			double* recvctrs = malloc(sizeof(double) * cluster * dimension);
			int* newcnts = malloc(sizeof(int) * cluster);
			int* recvcnts = malloc(sizeof(int) * cluster);
			memset(newctrs, 0, sizeof(double) * cluster * dimension);
			memset(newcnts, 0, sizeof(int) * cluster);
			for(i = 0; i < num; i++) {
				int lb = -1;
				double d2D = MAX_DIST;
				for(j = 0; j < cluster;j++) {
					double tmpDist = 
						dist2D(tdC+j*dimension, tdBuf+i*dimension);		
					if (tmpDist < d2D) {
						d2D = tmpDist;
						lb = j;
					}
							
				}	
				labels[i] = lb;
				newcnts[lb]++;
				for (j = 0; j < dimension; j++) {
					newctrs[lb * dimension + j] += tdBuf[i * dimension + j];
				}
			}
			MPI_Reduce(newctrs, recvctrs, cluster * dimension, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
			MPI_Reduce(newcnts, recvcnts, cluster, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
			int flag = 0;
			if (rank == MASTER){
				for (i = 0; i < cluster; i++) {
					for (j = 0; j < dimension; j++) {
						recvctrs[i * dimension + j] /= recvcnts[i];	
					}
				}
				double sum = 0;
				for (i = 0; i < cluster; i++) {
					sum += dist2D(recvctrs+i*dimension,tdC+i*dimension);
				}
				if (sum < TD_MIN)		
					flag = 1;
				free(tdC);
				tdC = recvctrs;
				/*for(i = 0; i<cluster ;i ++) {
					printf("cnts:%d\n",recvcnts[i]);
					printf("%lf,%lf\n",tdC[2*i],tdC[2*i+1]);
				}*/
			}
			MPI_Bcast (&flag,1,MPI_INT,MASTER,MPI_COMM_WORLD); 

			free(newcnts);
			free(recvcnts);
			free(newctrs);

			if(flag)
				break;


			MPI_Bcast (tdC,cluster * dimension,MPI_DOUBLE,MASTER,MPI_COMM_WORLD); 
			
		}
			
	}while(1);

	int* reclabels = malloc(sizeof(int)*lines);
	int offset = 0;	
	for(i = 0;i<numtasks;i++) {
		sendcnts[i] /= dimension;
		strides[i] = offset;
		offset += sendcnts[i];
	}
	MPI_Gatherv (labels,num,MPI_INT,reclabels,sendcnts,strides,MPI_INT,MASTER,MPI_COMM_WORLD);

	if(rank == MASTER && dna) {
		/* collect and output each cluster */

		for(i = 0;i<cluster;i++) {
			printf("Clster%d:\n",i);	
			printf("Center:\n");
			int j;
			for(j = 0;j < dimension; j++) {
				printf("%c,", dnaC[i*dimension+j]);
			}
			printf("\nPoints:\n");
			for(j = 0;j < lines;j++) {
				if(reclabels[j] == i) {
					int k;
					for(k = 0;k < dimension;k++) {
						printf("%c,",dnaSource[j*dimension+k]);
					}
					printf("\n");
				}
			}
		}
	} else if (rank == MASTER) {
		for(i = 0;i<cluster;i++) {
			printf("Clster%d:\n",i);	
			printf("Center:\n");
			int j;
			for(j = 0;j < dimension; j++) {
				printf("%lf,", tdC[i*dimension+j]);
			}
			printf("\nPoints:\n");
			for(j = 0;j < lines;j++) {
				if(reclabels[j] == i) {
					int k;
					for(k = 0;k < dimension;k++) {
						printf("%lf,",tdSource[j*dimension+k]);
					}
					printf("\n");
				}
			}
		}
	}
	
	free(labels);
	free(reclabels);
	free(sendcnts);
	free(strides);
	if(dna) {
		free(dnaSource);
		free(dnaBuf);
		free(dnaC);
	} else {
		free(tdSource);
		free(tdBuf);
		free(tdC);
	}

	MPI_Finalize();

}
