#include "PardisoSolver.h"

//#define PLOTS_PARDISO 

PardisoSolver::PardisoSolver()
{

}


PardisoSolver::~PardisoSolver()
{
	if (mtype == -1)
		return;
	/* -------------------------------------------------------------------- */
	/* ..  Termination and release of memory.                               */
	/* -------------------------------------------------------------------- */
	/*for (i = 0; i < num + 1; i++) {
		ia[i] -= 1;
		}
		for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}*/
	//ia.clear(); ja.clear();

	phase = -1;                 /* Release internal memory. */

	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&num, &ddum, ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);

}

void PardisoSolver::pardiso_init()
{
	mtype = -2;
	if (mtype == -1)
		throw std::runtime_error("Pardiso mtype not set.");

	/* -------------------------------------------------------------------- */
	/* ..  Setup Pardiso control parameters.                                */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}


	//printf("Variable number : %d\n", num);
	result.clear();
	result.resize(num, 0);

	nrhs = 1;
	error = 0;
	solver = 0;/* use sparse direct solver */
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

	if (error != 0)
	{
		if (error == -10)
			throw std::runtime_error("No license file found \n");
		if (error == -11)
			throw std::runtime_error("License is expired \n");
		if (error == -12)
			throw std::runtime_error("Wrong username or hostname \n");
	}
	// else
	//   printf("[PARDISO]: License check was successful ... \n");


	/* Numbers of processors, value of OMP_NUM_THREADS */
	 var = getenv("OMP_NUM_THREADS");
	 if(var != NULL)
	   sscanf( var, "%d", &num_procs );
	 else 
	  throw std::runtime_error("Set environment OMP_NUM_THREADS to 1");

	//num_procs = 12;
	iparm[2] = num_procs;

	maxfct = 1;		/* Maximum number of numerical factorizations.  */
	mnum = 1;         /* Which factorization to use. */

	msglvl = 0;         /* Print statistical information  */
	error = 0;         /* Initialize error flag */

	for (i = 0; i < num + 1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	//  /* -------------------------------------------------------------------- */
	//  /* Initialize the internal solver memory pointer. This is only */
	//  /* necessary for the FIRST call of the PARDISO solver. */
	//  /* -------------------------------------------------------------------- */
	/*for ( i = 0; i < 64; i++ )
	{
	pt[i] = 0;
	}*/
	phase = 11;
	//cout << "err" << endl;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&num, a.data(), ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
}
bool PardisoSolver::factorize()
{
	if (mtype == -1)
		throw std::runtime_error("Pardiso mtype not set.");
	/* -------------------------------------------------------------------- */
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */
	phase = 22;
	//  iparm[32] = 1; /* compute determinant */

	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&num, a.data(), ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);

	
#ifdef PLOTS_PARDISO
	printf("\nFactorization completed ... ");
#endif
	return (error == 0);
}


void PardisoSolver::pardiso_solver()
{
	

#ifdef PLOTS_PARDISO
	/* -------------------------------------------------------------------- */
	/* ..  pardiso_chkvec(...)                                              */
	/*     Checks the given vectors for infinite and NaN values             */
	/*     Input parameters (see PARDISO user manual for a description):    */
	/*     Use this functionality only for debugging purposes               */
	/* -------------------------------------------------------------------- */

	pardiso_chkvec(&numRows, &nrhs, rhs.data(), &error);
	if (error != 0) {
		printf("\nERROR  in right hand side: %d", error);
		exit(1);
	}

	/* -------------------------------------------------------------------- */
	/* .. pardiso_printstats(...)                                           */
	/*    prints information on the matrix to STDOUT.                       */
	/*    Use this functionality only for debugging purposes                */
	/* -------------------------------------------------------------------- */

	pardiso_printstats(&mtype, &numRows, a.data(), ia.data(), ja.data(), &nrhs, rhs.data(), &error);
	if (error != 0) {
		printf("\nERROR right hand side: %d", error);
		exit(1);
	}

#endif
	/* -------------------------------------------------------------------- */
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */
	phase = 33;

	iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&num, a.data(), ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, rhs.data(), result.data(), &error, dparm);

#ifdef PLOTS_PARDISO
	printf("\nSolve completed ... ");
	printf("\nThe solution of the system is: ");
	for (i = 0; i < numRows; i++) {
		printf("\n x [%d] = % f", i, result.data()[i]);
	}
	printf("\n\n");
#endif
}

void PardisoSolver::free_numerical_factorization_memory()
{
	phase = 0;                 /* Release internal memory. */

	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&num, &ddum, ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);
}
