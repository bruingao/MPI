#ifndef _UTIL_H_
#define _UTIL_H_

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define MAX_LINE 9999

/*
 * read 2D points file 
 * file: file pointer 
 * array: store the content read from the file 
 */
void read2D(FILE* file, double* array);
/*
 * read DNA points file 
 * file: file pointer 
 * array: store the content read from the file 
 */
void readDNA(FILE* file, char* array, int dimension);


/*
 * generate DNA cluster centroids 
 * centers: centers generated 
 * source: sources in which centers generated from 
 * size: how many points in source
 * dimension: the dimension of each point 
 * cluster: the number of clusters 
 */
void generateDNACenter(char* centers, char* source, int size, int dimension, int cluster);
/*
 * generate 2D cluster centroids 
 * centers: centers generated 
 * source: sources in which centers generated from 
 * size: how many points in source
 * dimension: the dimension of each point 
 * cluster: the number of clusters 
 */
void generate2DCenter(double* centers, double* source, int size, int dimension, int cluster);


/*
 * compute the distance between two DNA points 
 * centers: one point 
 * point: another point 
 * dimension: the dimension of each point 
 */
int distDNA(char* center, char* point, int dimension); 
/*
 * compute the distance between two 2D points 
 * centers: one point 
 * point: another point 
 */
double dist2D(double* center, double* point);


/*
 * determin whether the source point are too similar to one of the centroids
 * centers: the centroids 
 * source: the source point 
 * dimension: the dimension of each point 
 * num: the number of points in the center array
 */
int tooSimilar(char* centers, char* source, int dimension, int num);
/*
 * determin whether the source point are too close to one of the centroids
 * centers: the centroids 
 * source: the source point 
 * num: the number of points in the center array
 */
int tooClose(double* centers, double* source, int num);

#endif 

