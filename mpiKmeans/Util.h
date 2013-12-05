#ifndef _UTIL_H_
#define _UTIL_H_

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define MAX_LINE 9999

void read2D(FILE* file, double* array);
void readDNA(FILE* file, char* array, int size);

void generateDNACenter(char* centers, char* source, int size, int dimension, int cluster);
void generate2DCenter(double* centers, double* source, int size, int dimension, int cluster);

int distDNA(char* center, char* point, int dimension); 
double dist2D(double* center, double* point);

int tooSimilar(char* centers, char* source, int dimension, int num);

#endif 

