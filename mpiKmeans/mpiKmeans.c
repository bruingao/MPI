#include "mpi.h"
#include <stdio.h>
#include <string.h>

#include "Util.h"

/* master node */
#define MASTER 0
#define DNA_MIN 1
#define TD_MIN 0.000001
#define MAX_DIST 2147483647 


/* DNA Strands characters */
char bases[4] = {'A','C','G','T'};

int main(int argc,char** argv)
{
	/* number of task, task rank */
	int numtasks, rank, sendcnt, recvcnt;

	/* point dimension, number of points*/
	int dimension, lines;
	
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

	/* whether this task is a dna clustering */
	int dna;
	/* how many clusters */
	int cluster;

	/* check the arguments */
	if (argc < 5) {
		printf("Usage: mpiKmeans <dna or 2d> <input file name> <lines> <clusters> <dimension>\n");
		exit(1);
	}

	/* check the arguments */
	if (!strcmp(argv[1],"dna") && argc < 6) {
		printf("Usage: mpiKmeans <dna or 2d> <input file name> <lines> <clusters> <dimension>\n");
		exit(1);
	}
	
	/* file name */
	filename = argv[2];
	/* how many points */
	lines = atoi(argv[3]);
	/* dna clustring or not */
	dna = !strcmp(argv[1], "dna");
	/* how many clusters */
	cluster = atoi(argv[4]);

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	if (dna) {
		/* initialize dna clustring parameters */
		dimension = atoi(argv[5]);
		sendcnt = (lines + numtasks - 1)/numtasks * dimension;
		recvcnt = sendcnt;
		dnaSource = malloc(sizeof(char) * dimension * lines);
		dnaBuf = malloc(sizeof(char) * recvcnt);
		dnaC = malloc(sizeof(char) * cluster * dimension);
	} else {
		/* initialize 2d clustring parameters */
		dimension = 2;
		sendcnt = (lines + numtasks - 1)/numtasks * dimension;
		recvcnt = sendcnt;
		tdSource = malloc(sizeof(double) * dimension * lines);
		tdBuf = malloc(sizeof(double) * recvcnt);
		tdC = malloc(sizeof(double) * cluster * dimension);
	}

	
	/* if it's master, read the data source and genterate centroids */
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

	/* compute the number and offset of data elements to send for each process */
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

	/* compute the number of points for each process */
	int num = sendcnts[rank]/dimension;
	/* allocate an array for labeling each point */
	labels = malloc(sizeof(int) * num);
	
	/* do k-means until converging */
	do {
		int j;
		/* dna clustering */
		if(dna) {
			/* the count of each character in every position of each point in every process */
			int* newcnts = malloc(sizeof(int) * cluster * dimension * 4);
			/* the count of each character in every position of each point in sum */
			int* recvcnts = malloc(sizeof(int) * cluster * dimension * 4);
			/* initialize to zero */
			memset(newcnts, 0, sizeof(int) * cluster * dimension * 4);
			/* the new centroids for every cluster computed */
			char* recvctrs = malloc(sizeof(char) * cluster * dimension);
			/* for every point, compute its distance from each centroid, choose the minimum one and label it */
			for(i = 0; i < num; i++) {
				/* initialize label and distance */
				int lb = -1;
				int dDNA = MAX_DIST;
				/* iterate each cluster centroid and update the label and distance */
				for(j = 0; j < cluster;j++) {
					int tmpDist = distDNA(dnaC+j*dimension, dnaBuf+i*dimension, dimension);		
					if (tmpDist < dDNA) {
						dDNA = tmpDist;
						lb = j;
					}
				}	
				/* iterate each position of the point and update the counts of each character */
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
				/* update the point's label */
				labels[i] = lb;
			}
			
			/* reduce the character counts in every position of each point on the master node */
			MPI_Reduce(newcnts, recvcnts, cluster * 4 * dimension, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
			/* set the terminate flag */
			int flag = 0;
			/* on the master node, compute the new centroid for each cluster */
			if(rank == MASTER) {
				/* iterate each cluster */
				for (i = 0; i < cluster; i++) {
					/* iterate each position of a centroid */
					for (j = 0; j < dimension; j++) {
						/* choose the character with most counts */
						int left = 
						recvcnts[i*dimension*4+4*j+0] > recvcnts[i*dimension*4+4*j+1]?0:1;
						int right = 
						recvcnts[i*dimension*4+4*j+2] > recvcnts[i*dimension*4+4*j+3]?2:3;

						int final = recvcnts[i*dimension*4+4*j+left] > 
						recvcnts[i*dimension*4+4*j+right]?left:right;

						/* update the centroid */
						recvctrs[i * dimension + j] = bases[final];
					}
				}
				/* initialize the difference of centroids */
				int dsum = 0;
				/* iterate each cluster centroid and update the difference of centroids */
				for(i =0;i<cluster;i++) {
					dsum += distDNA(recvctrs+dimension*i,dnaC+dimension*i,dimension);
				}
				/* see if it is time to terminate */
				if(dsum <= DNA_MIN)
					flag = 1;
				free(dnaC);
				/* update the centroids */
				dnaC = recvctrs;
			}
			/* broadcast the terminate flag */
			MPI_Bcast (&flag,1,MPI_INT,MASTER,MPI_COMM_WORLD); 
			free(newcnts);
			free(recvcnts);
			if(flag)
				break;
			
			/* broadcast the new centroids */
			MPI_Bcast (dnaC,cluster * dimension,MPI_CHAR,MASTER,MPI_COMM_WORLD); 
			
		} else {
			/* the new centroids on each process */
			double* newctrs = malloc(sizeof(double) * cluster * dimension);
			/* the received centroids on the master */
			double* recvctrs = malloc(sizeof(double) * cluster * dimension);
			/* the number of points on each cluster */
			int* newcnts = malloc(sizeof(int) * cluster);
			/* the received number of points on each cluster */
			int* recvcnts = malloc(sizeof(int) * cluster);
			/* initialization */
			memset(newctrs, 0, sizeof(double) * cluster * dimension);
			memset(newcnts, 0, sizeof(int) * cluster);
			/* iterate each point on one process */
			for(i = 0; i < num; i++) {
				/* initialize the label and distance */
				int lb = -1;
				double d2D = MAX_DIST;
				/* iterate each centroid and update the label to the one with minimum distance */
				for(j = 0; j < cluster;j++) {
					double tmpDist = 
						dist2D(tdC+j*dimension, tdBuf+i*dimension);		
					if (tmpDist < d2D) {
						d2D = tmpDist;
						lb = j;
					}
							
				}	
				/* update the label */
				labels[i] = lb;
				/* update the count */
				newcnts[lb]++;
				/* sum each point */
				for (j = 0; j < dimension; j++) {
					newctrs[lb * dimension + j] += tdBuf[i * dimension + j];
				}
			}
			/* reduce the new centroid sum to the master process */
			MPI_Reduce(newctrs, recvctrs, cluster * dimension, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
			/* reduce the cluster counts to the master process */
			MPI_Reduce(newcnts, recvcnts, cluster, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
			/* initilize the termination flag */
			int flag = 0;
			/* on the master process, compute the new centroids */
			if (rank == MASTER){
				/* iterate each cluster */
				for (i = 0; i < cluster; i++) {
					for (j = 0; j < dimension; j++) {
						recvctrs[i * dimension + j] /= recvcnts[i];	
					}
				}
				/* intilize the difference sum */
				double sum = 0;
				/* iterate each cluster and update the sum */
				for (i = 0; i < cluster; i++) {
					sum += dist2D(recvctrs+i*dimension,tdC+i*dimension);
				}
				/* if the difference is less than a threshold, set the termination flag */
				if (sum < TD_MIN)		
					flag = 1;
				free(tdC);
				tdC = recvctrs;
			}
			/* broadcast the temination flag */
			MPI_Bcast (&flag,1,MPI_INT,MASTER,MPI_COMM_WORLD); 

			free(newcnts);
			free(recvcnts);
			free(newctrs);

			if(flag)
				break;

			/* broadcast the new centroids */
			MPI_Bcast (tdC,cluster * dimension,MPI_DOUBLE,MASTER,MPI_COMM_WORLD); 
			
		}
			
	}while(1);

	/* gather all the labels of all the points on the master process */
	int* reclabels = malloc(sizeof(int)*lines);
	int offset = 0;	
	for(i = 0;i<numtasks;i++) {
		sendcnts[i] /= dimension;
		strides[i] = offset;
		offset += sendcnts[i];
	}
	MPI_Gatherv (labels,num,MPI_INT,reclabels,sendcnts,strides,MPI_INT,MASTER,MPI_COMM_WORLD);

	/* output each cluster including the centroid and all the points */
//	 if(rank == MASTER && dna) {
//	 	for(i = 0;i<cluster;i++) {
//	 		printf("Clster%d:\n",i);	
//	 		printf("Center:\n");
//	 		int j;
//			for(j = 0;j < dimension; j++) {
//				printf("%c,", dnaC[i*dimension+j]);
//	 		}
//	 		printf("\nPoints:\n");
//	 		for(j = 0;j < lines;j++) {
//	 			if(reclabels[j] == i) {
//	 				int k;
//	 				for(k = 0;k < dimension;k++) {
//	 					printf("%c,",dnaSource[j*dimension+k]);
//	 				}
//	 				printf("\n");
//	 			}
//	 		}
//	 	}
//	 } else if (rank == MASTER) {
//	 	for(i = 0;i<cluster;i++) {
//	 		printf("Clster%d:\n",i);	
//	 		printf("Center:\n");
//	 		int j;
//	 		for(j = 0;j < dimension; j++) {
//	 			printf("%lf,", tdC[i*dimension+j]);
//	 		}
//	 		printf("\nPoints:\n");
//	 		for(j = 0;j < lines;j++) {
//	 			if(reclabels[j] == i) {
//	 				int k;
//	 				for(k = 0;k < dimension;k++) {
//	 					printf("%lf,",tdSource[j*dimension+k]);
//	 				}
//	 				printf("\n");
//	 			}
//	 		}
//	 	}
//	 }
	
	/* gabage collection */
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
