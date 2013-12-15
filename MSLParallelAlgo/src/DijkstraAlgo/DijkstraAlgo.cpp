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
#include <malloc/malloc.h>

#include "DijkstraAlgo.h"
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

DijkstraAlgo::DijkstraAlgo(const char* in, const char* out)
{
	if (in != NULL)
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
}

DijkstraAlgo::~DijkstraAlgo()
{
	if (adj_matrix_ != NULL)
	{
		free(adj_matrix_);
	}

}

void DijkstraAlgo::input_binary(int rank, int numprocs, int* numprocs_for_vi, int* vi, int* ranki)
{

	std::ifstream myfile;
	myfile.open(input_path_, std::ios::in | std::ios::binary);
	myfile.read((char*) &adj_matrix_size_, sizeof(adj_matrix_size_));

	if(numprocs <= adj_matrix_size_) {
		adj_matrix_ = (int**) malloc(sizeof(*adj_matrix_) * adj_matrix_size_);
		for (unsigned int i = 0; i < adj_matrix_size_; i++)
		{
			*(adj_matrix_ + i) = (int*) malloc(
					sizeof(*(adj_matrix_ + i)) * adj_matrix_size_);
		}
		for (int i = 0; i < adj_matrix_size_; i++)
			myfile.read((char*) *(adj_matrix_ + i),
					sizeof(int) * adj_matrix_size_);

		myfile.close();
	} else {
		//Recomputing for vertex vi
		*numprocs_for_vi = numprocs / adj_matrix_size_;
		*vi = rank / (*numprocs_for_vi);
		*ranki = rank % (*numprocs_for_vi);


		const unsigned int num_vertex_pi = adj_matrix_size_ / (*numprocs_for_vi); //Number of vertices that one pi holding
		const unsigned int base = num_vertex_pi * (*ranki); // base index of process #rank
		const unsigned int extras_vertex =
				(*ranki) != (*numprocs_for_vi) - 1 ? 0 : adj_matrix_size_ % (*numprocs_for_vi);

		adj_matrix_ = (int**) malloc(sizeof(*adj_matrix_) * adj_matrix_size_);
		for (unsigned int i = 0; i < adj_matrix_size_; i++)
		{
			*(adj_matrix_ + i) = (int*) malloc(
					sizeof(*(adj_matrix_ + i)) * adj_matrix_size_);
		}


		for (int i = 0; i < adj_matrix_size_; i++) {
			myfile.seekg(sizeof(int) + sizeof(int) * (i*adj_matrix_size_ + base));
			myfile.read((char*) (*(adj_matrix_ + i) + base), sizeof(int) * (num_vertex_pi + extras_vertex));
		}
		//printf("Rank #%d handle column:%d -> %d\n", rank, base, base+ num_vertex_pi + extras_vertex);

	//	for(int i = 0; i<adj_matrix_size_; i++) {
	//		for(int j= base; j<base + num_vertex_pi+extras_vertex; j++) {
	//			printf("%d ", *(*(adj_matrix_ +i) +j));
	//		}
	//		printf("\n");
	//	}
		myfile.close();
	}
}

//void DijkstraAlgo::output(int rank, int numprocs, int u, int* l)
//{
//	if (output_path_ == NULL)
//	{
//		for (unsigned int v = 0; v < adj_matrix_size_; v++)
//		{
//			printf("Path: %d -> %d: %d \n", u, v, l[v]);
//		}
//
//	}
//	else
//	{
////		std::ofstream myfile;
////		myfile.open (output_path_);
////		myfile << "MST: ";
////
////		while (!mst_.empty()) {
////			myfile << mst_.front();
////			myfile << " ";
////			mst_.pop();
////		}
////		myfile <<"\n";
////		myfile << "Weight: ";
////		myfile << mst_weight_;
////		myfile <<"\n";
////		myfile.close();
//	}
//}

int DijkstraAlgo::perform(int rank, int numprocs)
{
	double start_time, end_time;

	//Sequence
	if (numprocs == 1)
	{

		start_time = MPI_Wtime();
		input_binary(rank, numprocs, NULL, NULL, NULL);

		for (unsigned int v = 0; v < adj_matrix_size_; v++)
		{
			dijkstra_single_sequence(adj_matrix_size_, adj_matrix_, v);
		}

		end_time = MPI_Wtime();
		printf("Time executing: %f\n", end_time - start_time);

	}
	else //Parallel
	{
		start_time = MPI_Wtime();
		int numprocs_for_vi; //Total task assign for vertex i
		int vi; //Vertex i
		int ranki; //Task will handle vertex i

		input_binary(rank, numprocs, &numprocs_for_vi, &vi, &ranki);

		//0-> n process - computers
		if (numprocs <= adj_matrix_size_)
		{

			int num_vertex_pi = adj_matrix_size_ / numprocs;
			int base = num_vertex_pi * rank;
			int extra_vertex =
					rank != numprocs - 1 ? 0 : adj_matrix_size_ % numprocs;

			for (int v = base; v < base + num_vertex_pi + extra_vertex; v++)
				dijkstra_single_sequence(adj_matrix_size_, adj_matrix_, rank);

			end_time = MPI_Wtime();
			printf("Time executing for rank %d: %f\n", rank,
					end_time - start_time);
		}
		else //n->n^2/ logn processes - computers
		{

			MPI_Comm comm;
			MPI_Comm_split( MPI_COMM_WORLD, vi, ranki, &comm );

			dijkstra_single_parallel(adj_matrix_size_, adj_matrix_, vi, ranki,
					numprocs_for_vi, comm);

			end_time = MPI_Wtime();
			printf("Time executing for rank %d: %f\n", rank,
					end_time - start_time);
		}

	}
	return EXIT_SUCCESS;
}

int DijkstraAlgo::dijkstra_single_sequence(unsigned int adj_matrix_size,
		int** adj_matrix, unsigned int s)
{
	int l[adj_matrix_size];
	bool V[adj_matrix_size];
	int total_vertex = 0;

	std::fill_n(V, adj_matrix_size, false);

	V[s] = true;
	l[s] = 0;
	total_vertex++;
	for (unsigned int v = 0; v < adj_matrix_size; v++)
	{
		if (!V[v])
		{
			l[v] = adj_matrix[s][v];
		}
	}

	while (total_vertex < adj_matrix_size)
	{
		int min_weight = MAX_INTEGER;
		int u = -1;
		for (unsigned int v = 0; v < adj_matrix_size; v++)
		{
			if (!V[v])
			{
				if (min_weight > l[v])
				{
					min_weight = l[v];
					u = v;
				}
			}
		}

		if (u == -1)
		{
			printf("Fail...");
			return EXIT_FAILURE;
		}
		V[u] = true;
		total_vertex++;
		for (unsigned int v = 0; v < adj_matrix_size; v++)
		{
			if (!V[v])
			{
				l[v] = std::min(l[v], l[u] + adj_matrix[u][v]);
			}
		}
	}

//	for (unsigned int v = 0; v < adj_matrix_size_; v++)
//	{
//		printf("%d -> %d : %d\n", s, v, l[v]);
//	}
	return EXIT_SUCCESS;
}

int DijkstraAlgo::dijkstra_single_parallel(unsigned int adj_matrix_size,
		int** adj_matrix, unsigned int s, int rank, int numprocs, MPI_Comm comm)
{

	const unsigned int num_vertex_pi = adj_matrix_size / numprocs; //Number of vertices that one pi holding
	const unsigned int base = num_vertex_pi * rank; // base index of process #rank
	const unsigned int extras_vertex =
			rank != numprocs - 1 ? 0 : adj_matrix_size % numprocs;

	int l[adj_matrix_size];
	bool V[adj_matrix_size];

	std::fill_n(V, adj_matrix_size, false);

	if (s >= base && s < base + num_vertex_pi + extras_vertex)
	{
		V[s] = true;
	}
	l[s] = 0;

	for (unsigned int v = base; v < base + num_vertex_pi + extras_vertex; v++)
	{
		if (!V[v])
		{
			l[v] = adj_matrix[s][v];
		}
	}

	for (unsigned int total_vertex = 1; total_vertex < adj_matrix_size;
			total_vertex++)
	{
		int min_weight = MAX_INTEGER;
		int u = MAX_INTEGER;
		for (unsigned int v = base; v < base + num_vertex_pi + extras_vertex;
				v++)
		{
			if (!V[v])
			{
				if (min_weight > l[v])
				{
					min_weight = l[v];
					u = v;
				}
			}
		}

		int buf[2]; // value : vertex
		buf[0] = u < MAX_INTEGER ? l[u] : MAX_INTEGER;
		buf[1] = u;

		MPI_Allreduce(buf, buf, 1, MPI_2INT, MPI_MINLOC, comm );

		u = buf[1];
		l[u] = buf[0];
		if (u >= base && u < base + num_vertex_pi + extras_vertex)
		{
			V[u] = true;
		}

		for (unsigned int v = base; v < base + num_vertex_pi + extras_vertex;
				v++)
		{
			if (!V[v])
			{
				l[v] = std::min(l[v], l[u] + adj_matrix[u][v]);
			}
		}
	}
//	if (rank == 0)
//		for (unsigned int v = 0; v < adj_matrix_size_; v++)
//		{
//			printf("%d -> %d : %d\n", s, v, l[v]);
//		}
	return EXIT_SUCCESS;
}

