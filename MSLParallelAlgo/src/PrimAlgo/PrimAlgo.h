/*
 * IO.h
 *
 *  Created on: Oct 31, 2012
 *      Author: caominhvu
 */

#ifndef __PRIM_ALGO_H__
#define __PRIM_ALGO_H__


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

class PrimAlgo:public Algo {
public:
	PrimAlgo(const char* in, const char* out);
	~PrimAlgo();
	int perform(int rank, int numprocs);
private:

	int** adj_matrix_;
	unsigned int adj_matrix_size_;
	std::queue<unsigned int> mst_;
	unsigned int mst_weight_;
	char* input_path_;
	char* output_path_;

protected:
	void input(int rank, int numprocs);
	void output(int rank, int numprocs);
	int prim_mst(unsigned int adj_matrix_size, int** adj_matrix, unsigned int* weight, std::queue<unsigned int>* mst, int rank, int numprocs);
	int prim_mst(unsigned int adj_matrix_size, int** adj_matrix, unsigned int* weight, std::queue<unsigned int>* mst);
};

#ifdef __cplusplus
}
#endif


#endif /* PRIM_ALGO_H_ */
