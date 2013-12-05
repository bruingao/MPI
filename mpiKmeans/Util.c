#include "Util.h"

void read2D(FILE* file, double* array) {
	int index = 0;

	if (file == NULL) {
		printf("File not open\n");
		exit(1);
	}
	
	
	while(fscanf(file,"%lf,%lf\n",array + 2 * index,array + 2 * index + 1) != EOF) {
		index++;
	}
	
}


void readDNA(FILE* file, char* array, int size) {
	int index = 0;

	if (file == NULL) {
		printf ("File not open\n");
		exit(1);
	}
	
	int cnt = 1;
	while (1) {
		int res;
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

int distDNA(char* center, char* point, int dimension) {
        int i;
        int sum = 0;
        for(i = 0;i<dimension;i++)
                sum += (center[i] != point[i]);
        return sum;
}

double dist2D(double* center, double* point) {
        return sqrt(pow(center[0]-point[0],2) + pow(center[1] - point[1],2));
}

int tooSimilar(char* centers, char* source, int dimension, int num) {
	int minDist = (int)(0.2 * dimension);
	int i;
	for(i = 0;i<num;i++){
		if(distDNA(centers+i*dimension, source, dimension) < minDist)
			return 1;
	}
	return 0;
	
}

int tooClose(double* centers, double* source, int num) {
	double minDist = 0.1; 
	int i;
	for(i = 0;i<num;i++){
		if(dist2D(centers+i*2, source) < minDist)
			return 1;
	}
	return 0;
}

