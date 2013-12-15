/*
 * PrimAlgo.cpp
 *
 *  Created on: Oct 31, 2012
 *      Author: caominhvu
 */


/*------------------------------------  Includes   -------------------------------------*/
/*******************************************|********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <string.h>
#include <mpi.h>

#include "PrimAlgo.h"
/*------------------------------------- Data Types -------------------------------------*/
/*******************************************|********************************************/

/*----------------------------------  Generic Macros   ---------------------------------*/
/*******************************************|********************************************/

/*-----------------------------------  Definitions  ------------------------------------*/
/*******************************************|********************************************/
#define MAX_INTEGER 32767
#define ROOT_RANK 0

/*---------------------------  Data Structure Declarations   ---------------------------*/
/*******************************************|********************************************/

/*---------------------------  Public Function Declarations  ---------------------------*/
/*******************************************|********************************************/

void PrimAlgo::input(int rank, int numprocs) {
	FILE* file;
	if((file = fopen(input_path_, "r"))==NULL)
	{
		printf("ERR: can't read file");
		return;
	}

	fscanf(file, "%d", &adj_matrix_size_); //Load matrix size

	adj_matrix_ = (int**) malloc(sizeof(*adj_matrix_) * adj_matrix_size_);
	for(unsigned int i=0; i< adj_matrix_size_; i++)
	{
		*(adj_matrix_ + i) = (int*) malloc(sizeof(*(adj_matrix_+i)) * adj_matrix_size_);
		for(unsigned int j=0; j<adj_matrix_size_; j++)
		{
			fscanf(file, "%d", *(adj_matrix_ + i) + j); //Load value for each elements in the matrix
		}
	}
	fclose(file);

//	if(rank == ROOT_RANK) {
//		printf("Matrix size: %d\n", adj_matrix_size_);
//		for(int i=0; i< adj_matrix_size_; i++)
//		{
//			for(int j=0; j< adj_matrix_size_; j++)
//			{
//				if (adj_matrix_[i][j] != MAX_INTEGER)
//				{
//					printf("%d   ", adj_matrix_[i][j]);
//				} else {
//					printf("∞   ");
//				}
//			}
//			printf("\n");
//		}
//	}
}

void PrimAlgo::output(int rank, int numprocs)
{
	if(output_path_ == NULL) {
		printf("MST: ");

		while (!mst_.empty()) {
			printf("%d ", mst_.front());
			mst_.pop();
		}
		printf("\nWeight: %d\n", mst_weight_);
	} else {
		std::ofstream myfile;
		myfile.open (output_path_);
		myfile << "MST: ";

		while (!mst_.empty()) {
			myfile << mst_.front();
			myfile << " ";
			mst_.pop();
		}
		myfile <<"\n";
		myfile << "Weight: ";
		myfile << mst_weight_;
		myfile <<"\n";
		myfile.close();
	}
}

PrimAlgo::PrimAlgo(const char* in, const char* out)
{
	if(in != NULL)
	{
		input_path_ = (char*) malloc(sizeof(input_path_) * (strlen(in) + 1));
		strcpy(input_path_, in);
	}
	else
	{
		input_path_ = NULL;
	}

	if (out != NULL)
	{
		output_path_ = (char*) malloc(sizeof(output_path_) * (strlen(out) + 1));
		strcpy(output_path_, out);
	}
	else
	{
		output_path_ = NULL;
	}

	adj_matrix_ = NULL;
	adj_matrix_size_ = 0;
	mst_weight_ = 0;
}

PrimAlgo::~PrimAlgo()
{

	if(adj_matrix_ != NULL)
	{
		free(adj_matrix_);
	}
}

int PrimAlgo::perform(int rank, int numprocs)
{
	double start_time, end_time;

	if(numprocs == 1) {
		input(rank, numprocs);
		start_time = MPI_Wtime();
		int ret = prim_mst(adj_matrix_size_, adj_matrix_, &mst_weight_, &mst_);
		end_time = MPI_Wtime();

		printf("Time executing: %f\n", end_time - start_time);

		if (ret == EXIT_SUCCESS) {
			output(rank, numprocs);
		}
	} else {
		if(rank == 0) {
			input(rank, numprocs);
		}
		MPI_Bcast(&adj_matrix_size_, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

		if(rank == 0) {
			for(int r=1; r< numprocs; r++)
			{
				int base = (adj_matrix_size_/numprocs) *r;
				int num_vertex_pi = adj_matrix_size_/numprocs;
				int remainder = r!= numprocs -1 ? 0 : adj_matrix_size_ % numprocs;
				for(int n = base; n<base + num_vertex_pi + remainder; n++)
					MPI_Send(adj_matrix_[n], adj_matrix_size_, MPI_INT, r, r, MPI_COMM_WORLD );
			}
		} else {
			MPI_Status status;
			int base = (adj_matrix_size_/numprocs) *rank;
			int num_vertex_pi = adj_matrix_size_/numprocs;
			int remainder = rank!= numprocs -1 ? 0 : adj_matrix_size_ % numprocs;
			adj_matrix_ = (int**) malloc(sizeof(*adj_matrix_) * adj_matrix_size_);
			for(unsigned int i = 0; i< adj_matrix_size_; i++)
			{
				*(adj_matrix_ + i) = (int*) malloc(sizeof(*(adj_matrix_+i)) * adj_matrix_size_);
			}
			for(int n = base; n<base+num_vertex_pi+ remainder; n++)
				MPI_Recv(adj_matrix_[n], adj_matrix_size_, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		if(rank == 0)
		{
			start_time = MPI_Wtime();
		}
		prim_mst(adj_matrix_size_, adj_matrix_, &mst_weight_, &mst_, rank, numprocs);
		if(rank == 0)
		{
			end_time = MPI_Wtime();
			printf("Time executing: %f\n", end_time - start_time);
			output(rank, numprocs);
		}


	}
	return EXIT_SUCCESS;
}

int PrimAlgo::prim_mst(unsigned int adj_matrix_size, int** adj_matrix, unsigned int* mst_weight, std::queue<unsigned int>* mst)
{
	bool V[adj_matrix_size]; //List of vertex V[v] = true if v ∈ MST, false if vice versa
	int d[adj_matrix_size]; //List of the minimum weight from one vertex to MST, d[u] = min{w(u,v)|v ∈ V - MST}

	std::fill_n(V, adj_matrix_size, false);
	int r = random() % adj_matrix_size;
	V[r] = true;
	d[r] = 0;
	mst->push(r); //Just store for final result

	for (unsigned int v = 0; v < adj_matrix_size; v++)
	{
		if (!V[v])
		{
			d[v] = adj_matrix[r][v];

		}
	}


	for (unsigned int MST_size = 1; MST_size < adj_matrix_size; MST_size++)
	{

		int min_weight = MAX_INTEGER;
		int u = -1;
		for (unsigned int v = 0; v < adj_matrix_size; v++) //Finding vertex that have weight to MST is smallest
		{
			if (!V[v])
			{
				if (min_weight > d[v])
				{
					min_weight = d[v];
					u = v;
				}
			}
		}

		if (u == -1) return EXIT_FAILURE;

		V[u] = true; //put u into MST
		mst->push(u);
		(*mst_weight) += d[u];

		for (unsigned int v = 0; v < adj_matrix_size; v++) //Update d[] values
		{
			if (!V[v])
			{
				d[v] = std::min(d[v], adj_matrix[u][v]);
			}
		}
	}
	return EXIT_SUCCESS;
}

int PrimAlgo::prim_mst(unsigned int adj_matrix_size, int** adj_matrix, unsigned int* mst_weight, std::queue<unsigned int>* mst, int rank, int numprocs)
{
	const unsigned int base = (adj_matrix_size/numprocs) *rank;
	const unsigned int num_vertex_pi = adj_matrix_size/numprocs;
	const unsigned int remainder = rank!= numprocs -1 ? 0 : adj_matrix_size % numprocs;
	bool V[adj_matrix_size]; //List of vertex V[v] = true if v ∈ MST, false if vice versa
	int d[adj_matrix_size]; //List of the minimum weight from one vertex to MST, d[u] = min{w(u,v)|v ∈ V - MST}

	std::fill_n(V, adj_matrix_size, false);
	unsigned int r;
	if(rank ==0)
	{
		r = 0;
		mst->push(r); //Just store for final result
	}
	MPI_Bcast(&r, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
	if(r >= base && r<base + num_vertex_pi + remainder)
	{
		V[r] = true;
		d[r] = 0;
	}

	for (unsigned int v = base; v < base + adj_matrix_size/numprocs + remainder; v++)
	{
		if (!V[v])
		{
			d[v] = adj_matrix[v][r];
		}
	}
	for (unsigned int MST_size = 1; MST_size < adj_matrix_size; MST_size++)
	{
		int min_weight = MAX_INTEGER;
		unsigned int u = MAX_INTEGER;
		for (unsigned int v = base; v < base + num_vertex_pi + remainder; v++) //Finding vertex that have weight to MST is smallest
		{
			if (!V[v])
			{
				if (min_weight > d[v])
				{
					min_weight = d[v];
					u = v;
				}
			}
		}

		int buf[2];
		buf[0] = u< MAX_INTEGER ? d[u] : MAX_INTEGER;
		buf[1] = u;
		MPI_Allreduce(buf, buf, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

		u = buf[1];
		if(u >= base && u<base + num_vertex_pi + remainder)
		{
			V[u] = true;
		}
		if(rank == 0)
		{
			mst->push(u);
			(*mst_weight) += buf[0];
		}

		for (unsigned int v = base; v < base + num_vertex_pi + remainder; v++) //Update d[] values
		{
			if (!V[v])
			{
				d[v] = std::min(d[v], adj_matrix[v][u]);
			}
		}
	}
	return EXIT_SUCCESS;
}



