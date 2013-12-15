/*
 * MPIMain.cpp
 *
 *  Created on: Oct 31, 2012
 *      Author: caominhvu
 */
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "DijkstraAlgo/DijkstraAlgo.h"

#define MAX_INTEGER 32767

using namespace std;

/**

 * Ansi C "itoa" based on Kernighan & Ritchie's "Ansi C"

 * with slight modification to optimize for specific architecture:

 */


void strreverse(char* begin, char* end)
{

	char aux;

	while (end > begin)

		aux = *end, *end-- = *begin, *begin++ = aux;

}

void itoa(int value, char* str, int base)
{

	static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";

	char* wstr = str;

	int sign;

	div_t res;

	// Validate base

	if (base < 2 || base > 35)
	{
		*wstr = '\0';
		return;
	}

	// Take care of sign

	if ((sign = value) < 0)
		value = -value;

	// Conversion. Number is reversed.

	do
	{

		res = div(value, base);

		*wstr++ = num[res.rem];

	} while (value = res.quot);

	if (sign < 0)
		*wstr++ = '-';

	*wstr = '\0';

	// Reverse string

	strreverse(str, wstr - 1);

}
void createFile(char* file_path, int size)
{

	int** adj_matrix;

	adj_matrix = (int**) malloc(sizeof(*adj_matrix) * size);
	for (unsigned int i = 0; i < size; i++)
	{
		*(adj_matrix + i) = (int*) malloc(
				sizeof(*(adj_matrix + i)) * size);
	}

	for(int i = 0; i<size; i++)
	{
		for(int j = 0; j< size; j++)
		{
			adj_matrix[i][j] = (i == j) ? 0 :  (std::rand() % MAX_INTEGER) /size;
		}
	}

	ofstream myfile;
	myfile.open (file_path, ios::out | ios::binary);
	myfile.write((char*) &size, sizeof(size));
	for(int i = 0; i< size; i++)
		myfile.write((char*) *(adj_matrix + i), sizeof(int) * size);
	myfile.close();
	free(adj_matrix);
}

void readFile(char* file_path) {
	int size;
	fstream myfile;
	myfile.open(file_path, ios::in | ios:: binary);
	myfile.read((char*) &size, sizeof(size));
	printf("%s Matrix size: %d\n",file_path,  size);

	int** adj_matrix;
	adj_matrix = (int**) malloc(sizeof(*adj_matrix) * size);
	for (unsigned int i = 0; i < size; i++)
	{
		*(adj_matrix + i) = (int*) malloc(sizeof(*(adj_matrix + i)) * size);
	}

	for(int i = 0; i<size; i++)
		myfile.read((char*)  *(adj_matrix+i), sizeof(int) * size);

	for(int i = 0; i<size; i++)
	{
		for(int j = 0; j<size; j++) {
			printf("%d ", adj_matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	myfile.close();
	free(adj_matrix);
}

void generateTestCase()
{
	for (int i = 1; i < 14; i++)
	{
		char str[8];
		int val = pow(2, i);
		itoa(val, str, 10);
		printf("%d %s \n", val, str);
		char link[8];
		strcpy(link, "./");
		strcpy(link, str);
		createFile(link, val);
		//readFile(link);
	}
}


int main(int argc, char** argv) {
//	generateTestCase();

	int rank, numprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

	if(argv[1] == NULL) {
		printf("No data\n");
		return EXIT_FAILURE;
	}
	Algo *algo = new DijkstraAlgo(argv[1], argv[2]);
	algo->perform(rank, numprocs);
	free(algo);

	MPI_Finalize();
	return EXIT_SUCCESS;
}
