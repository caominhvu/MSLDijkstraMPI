/*
 * IO.h
 *
 *  Created on: Oct 31, 2012
 *      Author: caominhvu
 */

#ifndef __DIJKSTRA_ALGO_H__
#define __DIJKSTRA_ALGO_H__


/*------------------------------------  Includes   -------------------------------------*/
/*******************************************|********************************************/
#include <queue>
#include "../Algo.h"
/*------------------------------------- Data Types -------------------------------------*/
/*******************************************|********************************************/

/*----------------------------------  Generic Macros   ---------------------------------*/
/*******************************************|********************************************/

/*-----------------------------------  Definitions  ------------------------------------*/
/*******************************************|********************************************/

/*---------------------------  Data Structure Declarations   ---------------------------*/
/*******************************************|********************************************/

/*---------------------------  Public Function Declarations  ---------------------------*/
/*******************************************|********************************************/

#ifdef __cplusplus
extern "C" {
#endif

class DijkstraAlgo:public Algo {
public:
	DijkstraAlgo(const char* in, const char* out);
	~DijkstraAlgo();
	int perform(int rank, int numprocs);
private:
	int** adj_matrix_;
	int adj_matrix_size_;
	char* input_path_;
	char* output_path_;

protected:
//	void input_text();
//	void input_binary();
	void input_binary(int rank, int numprocs, int*, int*, int*);
//	void output(int rank, int numprocs, int u, int* l);
	int dijkstra_single_sequence(unsigned int adj_matrix_size, int** adj_matrix, unsigned int s);
	int dijkstra_single_parallel(unsigned int adj_matrix_size, int** adj_matrix, unsigned int s, int rank, int numprocs, MPI_Comm comm);
};

#ifdef __cplusplus
}
#endif


#endif /* Dijkstra_ALGO_H_ */
