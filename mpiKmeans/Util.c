#include "Util.h"

/*
 * read 2D points file 
 * file: file pointer 
 * array: store the content read from the file 
 */
void read2D(FILE* file, double* array) {
	int index = 0;
	/* file not open */
	if (file == NULL) {
		printf("File not open\n");
		exit(1);
	}
	
	/* scan the whole file line by line */	
	while(fscanf(file,"%lf,%lf\n",array + 2 * index,array + 2 * index + 1) != EOF) {
		index++;
	}
	
}

/*
 * read DNA points file 
 * file: file pointer 
 * array: store the content read from the file 
 */
void readDNA(FILE* file, char* array, int size) {
	int index = 0;

	if (file == NULL) {
		printf ("File not open\n");
		exit(1);
	}
	
	int cnt = 1;
	while (1) {
		int res;
		/* scan the whole file number by number */	
		if(cnt%size == 0) {
			res = fscanf(file,"%c\n", array+cnt-1);
		} else {
			res = fscanf(file,"%c,", array+cnt-1);	
		}
		if (res == EOF)
			break;
		cnt++;
	}
}

/*
 * generate DNA cluster centroids 
 * centers: centers generated 
 * source: sources in which centers generated from 
 * size: how many points in source
 * dimension: the dimension of each point 
 * cluster: the number of clusters 
 */
void generateDNACenter(char* centers, char* source, int size, int dimension, int cluster) {
	int i;
	int index = ((int)rand())%size;
	memcpy(centers,source+index*dimension,dimension);
	index = ((int)rand())%size;
	for (i = 1;i < cluster;i++) {
		while(tooSimilar(centers, source+index*dimension, dimension, i)) {
			index = ((int)rand())%size;
		}
		memcpy(centers+i*dimension, source+index*dimension, dimension);
	}
}

/*
 * generate 2D cluster centroids 
 * centers: centers generated 
 * source: sources in which centers generated from 
 * size: how many points in source
 * dimension: the dimension of each point 
 * cluster: the number of clusters 
 */
void generate2DCenter(double* centers, double* source, int size, int dimension, int cluster) {
	int i;
	int index = ((int)rand())%size;
	memcpy(centers, source+index*dimension,dimension * sizeof(double));
	for (i = 0;i < cluster;i++) {
		while(tooClose(centers, source+index*dimension, i)) {
			index = ((int)rand())%size;
		}
		memcpy(centers + i*dimension, source+index*dimension,dimension * sizeof(double));
	}
}

/*
 * compute the distance between two DNA points 
 * centers: one point 
 * point: another point 
 * dimension: the dimension of each point 
 */
int distDNA(char* center, char* point, int dimension) {
        int i;
        int sum = 0;
        for(i = 0;i<dimension;i++)
                sum += (center[i] != point[i]);
        return sum;
}

/*
 * compute the distance between two 2D points 
 * centers: one point 
 * point: another point 
 */
double dist2D(double* center, double* point) {
        return sqrt(pow(center[0]-point[0],2) + pow(center[1] - point[1],2));
}

/*
 * determin whether the source point are too similar to one of the centroids
 * centers: the centroids 
 * source: the source point 
 * dimension: the dimension of each point 
 * num: the number of points in the center array
 */
int tooSimilar(char* centers, char* source, int dimension, int num) {
	int minDist = (int)(0.2 * dimension);
	int i;
	for(i = 0;i<num;i++){
		if(distDNA(centers+i*dimension, source, dimension) < minDist)
			return 1;
	}
	return 0;
	
}

/*
 * determin whether the source point are too close to one of the centroids
 * centers: the centroids 
 * source: the source point 
 * num: the number of points in the center array
 */
int tooClose(double* centers, double* source, int num) {
	double minDist = 0.1; 
	int i;
	for(i = 0;i<num;i++){
		if(dist2D(centers+i*2, source) < minDist)
			return 1;
	}
	return 0;
}

