#pragma once

#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

extern "C" {
	/* PARDISO prototype. */
	void pardisoinit(void   *, int    *, int *, int *, double *, int *);
	void pardiso(void   *, int    *, int *, int *, int *, int *,
		double *, int    *, int *, int *, int *, int *,
		int *, double *, double *, int *, double *);
	void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
	void pardiso_chkvec(int *, int *, double *, int *);
	void pardiso_printstats(int *, int *, double *, int *, int *, int *,
		double *, int *);
}

class PardisoSolver
{
public:
	PardisoSolver();
	~PardisoSolver();
	
	void pardiso_init();
	bool factorize();
	void pardiso_solver();
	void free_numerical_factorization_memory();

	 vector<double> result;

	 vector<int> ia;
	 vector<int> ja;
	 vector<double> a;
	 vector<double> rhs;

	 int nnz;
	 int num;

protected:
	
	//pardiso stuff
	/*
	1: real and structurally symmetric, supernode pivoting
	2: real and symmetric positive definite
	-2: real and symmetric indefinite, diagonal or Bunch-Kaufman pivoting
	11: real and nonsymmetric, complete supernode pivoting
	*/
	int mtype;       /* Matrix Type */

	// Remember if matrix is symmetric or not, to
	// decide whether to eliminate the non-upper-
	// diagonal entries from the input II,JJ,SS
	bool is_symmetric;
	bool is_upper_half;

	int nrhs;     /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	int iparm[64];
	double   dparm[64];
	int maxfct, mnum, phase, error, msglvl, solver;
	/* Number of processors. */
	int      num_procs;
	/* Auxiliary variables. */
	char    *var;
	int i, k;
	double ddum;          /* Double dummy */
	int idum;         /* Integer dummy. */
};

