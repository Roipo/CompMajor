//
//  PardisoSolver.cpp
//
//  Created by Olga Diamanti on 07/01/15.
//  Copyright (c) 2015 Olga Diamanti. All rights reserved.
//

#include "PardisoSolver.h"
#include <igl/sortrows.h>
#include <igl/unique.h>
#include <igl/unique_rows.h>
#include <igl/matlab_format.h>


using namespace std;
//#define PLOTS_PARDISO


template <typename vectorTypeI, typename vectorTypeS>
PardisoSolver<vectorTypeI,vectorTypeS>::PardisoSolver():
mtype(-1)
{}

template <typename vectorTypeI, typename vectorTypeS>
void PardisoSolver<vectorTypeI,vectorTypeS>::set_type(int _mtype , bool is_upper_half)
{
  if ((_mtype !=-2) && (_mtype !=2) && (_mtype !=1) && (_mtype !=11))
	  throw std::runtime_error(std::string("Pardiso mtype not supported. mtype = ")+std::to_string( _mtype));
 
  mtype = _mtype;
  // As per https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/283738
  // structurally symmetric need the full matrix
  // structure to be passed, not only the upper
  // diagonal part. So is_symmetric should be set to
  // false in that case.
  is_symmetric = (mtype ==2) ||(mtype ==-2);
  this->is_upper_half = is_upper_half;
  if(!is_symmetric && is_upper_half)
	  throw std::runtime_error("Using upper half is only possible if the matrix is symmetric.");
  init();
}

template <typename vectorTypeI, typename vectorTypeS>
void PardisoSolver<vectorTypeI,vectorTypeS>::init()
{
   if (mtype ==-1)
     throw std::runtime_error("Pardiso mtype not set.");

   /* -------------------------------------------------------------------- */
   /* ..  Setup Pardiso control parameters.                                */
   /* -------------------------------------------------------------------- */

   error = 0;
   solver=0;/* use sparse direct solver */
   pardisoinit (pt,  &mtype, iparm);

   if (error != 0)
   {
     if (error == -10 )
		 throw std::runtime_error("No license file found \n");
     if (error == -11 )
		 throw std::runtime_error("License is expired \n");
     if (error == -12 )
		 throw std::runtime_error("Wrong username or hostname \n");
		 }
 // else
 //   printf("[PARDISO]: License check was successful ... \n");
  
  
  /* Numbers of processors, value of OMP_NUM_THREADS */
  // var = getenv("OMP_NUM_THREADS");
  // if(var != NULL)
  //   sscanf( var, "%d", &num_procs );
  // else 
 	//  throw std::runtime_error("Set environment OMP_NUM_THREADS to 1");
	 
  iparm[2]  = num_procs;
  
  maxfct = 1;		/* Maximum number of numerical factorizations.  */
  mnum   = 1;         /* Which factorization to use. */
  
  msglvl = 0;         /* Print statistical information  */
  error  = 0;         /* Initialize error flag */
  
  
  //  /* -------------------------------------------------------------------- */
  //  /* .. Initialize the internal solver memory pointer. This is only */
  //  /* necessary for the FIRST call of the PARDISO solver. */
  //  /* -------------------------------------------------------------------- */
  //  for ( i = 0; i < 64; i++ )
  //  {
  //    pt[i] = 0;
  //  }
  
  
}

template <typename vectorTypeI, typename vectorTypeS>
void PardisoSolver<vectorTypeI,vectorTypeS>::update_a(const vectorTypeS &SS_)
 {
   if (mtype ==-1)
	   throw std::runtime_error("Pardiso mtype not set.");
 vectorTypeS SS;
  MKL_INT numel = SS_.size();
  MKL_INT ntotal = numel +numRows;
  SS.resize(ntotal);
  for (MKL_INT k = 0; k < numel; ++k)
    SS[k] = SS_[k];
  for (MKL_INT k = 0; k < numRows; ++k)
    SS[numel+k] = 0;
   vectorTypeS SS_true = SS;

   //if the matrix is symmetric, only store upper triangular part
   if (is_symmetric && !is_upper_half)
   {
    SS_true.resize(lower_triangular_ind.size(),1);
    for (MKL_INT i = 0; i<lower_triangular_ind.size();++i)
      SS_true[i] = SS[lower_triangular_ind[i]];
  }
  
	      
  for (MKL_INT i=0; i<a.rows(); ++i)
  {
    a(i) = 0;
    for (MKL_INT j=0; j<iis[i].size(); ++j)
      a(i) += SS_true[iis[i](j)];
  }
}

template <typename vectorTypeI, typename vectorTypeS>
void PardisoSolver<vectorTypeI,vectorTypeS>::set_pattern(const vectorTypeI &II_,
                                                         const vectorTypeI &JJ_,
                                                         const vectorTypeS &SS_)


{
  if (mtype ==-1)
  	throw std::runtime_error("Pardiso mtype not set.");
  numRows = 0;
  for (MKL_INT i=0; i<II_.size(); ++i)
    if (II_[i] > numRows )
      numRows = II_[i];
  numRows ++;
  
  
  //make sure diagonal terms are included, even as zeros (pardiso claims this is necessary for best performance, and it also prevents it from occasionally hanging with zero diagonal elements)
  vectorTypeI II;
  vectorTypeI JJ;
  vectorTypeS SS;
  MKL_INT numel = II_.size();
  MKL_INT ntotal = numel +numRows;
  II.resize(ntotal);
  JJ.resize(ntotal);
  SS.resize(ntotal);
  for (MKL_INT k = 0; k < numel; ++k)
  {
    II[k] = II_[k];
    JJ[k] = JJ_[k];
    SS[k] = SS_[k];
  }
  for (MKL_INT k = 0; k < numRows; ++k)
  {
    II[numel+k] = k;
    JJ[numel+k] = k;
    SS[numel+k] = 0;
  }
  
  
  vectorTypeS SS_true = SS;
  Eigen::Matrix<MKL_INT, Eigen::Dynamic, Eigen::Dynamic> M0;
  //if the matrix is symmetric, only store upper triangular part
  if (is_symmetric && !is_upper_half)
  {
    lower_triangular_ind.reserve(II.size()/2 + numRows / 2 + 1);
    for (MKL_INT i = 0; i<II.size();++i)
      if (II[i]<=JJ[i])
        lower_triangular_ind.push_back(i);
    M0.resize(lower_triangular_ind.size(),3);
    SS_true.resize(lower_triangular_ind.size(),1);
    for (MKL_INT i = 0; i<lower_triangular_ind.size();++i)
    {
      M0.row(i)<< II[lower_triangular_ind[i]], JJ[lower_triangular_ind[i]], i;
      SS_true[i] = SS[lower_triangular_ind[i]];
    }
  }
  else
  {
    M0.resize(II.size(),3);
    for (MKL_INT i = 0; i<II.size();++i)
      M0.row(i)<< II[i], JJ[i], i;
  }
  
  //temps
  Eigen::Matrix<MKL_INT, Eigen::Dynamic, Eigen::Dynamic> t;
  Eigen::Matrix<MKL_INT, Eigen::Dynamic, 1> tI;
  
  Eigen::Matrix<MKL_INT, Eigen::Dynamic, Eigen::Dynamic> M_;
  igl::sortrows(M0, true, M_, tI);
  
  MKL_INT si,ei,currI;
  si = 0;
  while (si<M_.rows())
  {
    currI = M_(si,0);
    ei = si;
    while (ei<M_.rows() && M_(ei,0) == currI)
      ++ei;
    igl::sortrows(M_.block(si, 1, ei-si, 2).eval(), true, t, tI);
    M_.block(si, 1, ei-si, 2) = t;
    si = ei;
  }
  
  Eigen::Matrix<MKL_INT, Eigen::Dynamic, Eigen::Dynamic> M;
  Eigen::Matrix<MKL_INT, Eigen::Dynamic, 1> IM_;
  igl::unique_rows(M_.leftCols(2).eval(), M, IM_, tI);
  MKL_INT numUniqueElements = M.rows();
  iis.resize(numUniqueElements);
  for (MKL_INT i=0; i<numUniqueElements; ++i)
  {
    si = IM_(i);
    if (i<numUniqueElements-1)
      ei = IM_(i+1);
    else
      ei = M_.rows();
    iis[i] = M_.block(si, 2, ei-si, 1);
  }
  
  a.resize(numUniqueElements, 1);
  for (MKL_INT i=0; i<numUniqueElements; ++i)
  {
    a(i) = 0;
    for (MKL_INT j=0; j<iis[i].size(); ++j)
      a(i) += SS_true[iis[i](j)];
  }
  
  // now M_ and elements in sum have the row, column and indices in sum of the
  // unique non-zero elements in B1
  ia.setZero(numRows+1,1);ia(numRows) = numUniqueElements+1;
  ja = M.col(1).array()+1;
  currI = -1;
  for (MKL_INT i=0; i<numUniqueElements; ++i)
  {
    if(currI != M(i,0))
    {
      ia(M(i,0)) = i+1;//do not subtract 1
      currI = M(i,0);
    }
  }
  
#ifdef PLOTS_PARDISO
  printf("ia: ");
  for (int i=0; i<ia.size(); ++i)
    printf("%d ",ia[i]);
  printf("\n\n");
  
  printf("ja: ");
  for (int i=0; i<ja.size(); ++i)
    printf("%d ",ja[i]);
  printf("\n\n");
  
  printf("a: ");
  for (int i=0; i<a.size(); ++i)
    printf("%.2f ",a[i]);
  printf("\n\n");
#endif
  
  // matrix in CRS can be expressed with ia, ja and iis
  
}

template <typename vectorTypeI, typename vectorTypeS>
void PardisoSolver<vectorTypeI,vectorTypeS>::analyze_pattern()
{
  if (mtype ==-1)
  	throw std::runtime_error("Pardiso mtype not set.");
  
#ifdef PLOTS_PARDISO
  /* -------------------------------------------------------------------- */
  /*  .. pardiso_chk_matrix(...)                                          */
  /*     Checks the consistency of the given matrix.                      */
  /*     Use this functionality only for debugging purposes               */
  /* -------------------------------------------------------------------- */
  
  pardiso_chkmatrix  (&mtype, &numRows, a.data(), ia.data(), ja.data(), &error);
  if (error != 0) 
  	throw std::runtime_error(std::string("\nERROR in consistency of matrix: ") + std::to_string(error));
#endif
  /* -------------------------------------------------------------------- */
  /* ..  Reordering and Symbolic Factorization.  This step also allocates */
  /*     all memory that is necessary for the factorization.              */
  /* -------------------------------------------------------------------- */
  phase = 11;
  
  //pardiso (pt, &maxfct, &mnum, &mtype, &phase,
  //         &numRows, a.data(), ia.data(), ja.data(), &idum, &nrhs,
  //         iparm, &msglvl, &ddum, &ddum, &error);
  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &numRows, a.data(), ia.data(),
      ja.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);


  if (error != 0) 
  	throw std::runtime_error(std::string("\nERROR during symbolic factorization: ") + std::to_string(error));
#ifdef PLOTS_PARDISO
  printf("\nReordering completed ... ");
  printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
#endif
  
}

template <typename vectorTypeI, typename vectorTypeS>
bool PardisoSolver<vectorTypeI,vectorTypeS>::factorize()
{
  if (mtype ==-1)
  	throw std::runtime_error("Pardiso mtype not set.");
  /* -------------------------------------------------------------------- */
  /* ..  Numerical factorization.                                         */
  /* -------------------------------------------------------------------- */
  phase = 22;
  //  iparm[32] = 1; /* compute determinant */
  
  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
           &numRows, a.data(), ia.data(), ja.data(), &idum, &nrhs,
           iparm, &msglvl, &ddum, &ddum, &error);
  
  if (error != 0) 
  	throw std::runtime_error(std::string("\nERROR during numerical factorization: ")+std::to_string(error));
#ifdef PLOTS_PARDISO
  printf ("\nFactorization completed ... ");
#endif
  return (error ==0);
}

template <typename vectorTypeI, typename vectorTypeS>
void PardisoSolver<vectorTypeI,vectorTypeS>::solve(Eigen::VectorXd &rhs,
                                                   Eigen::VectorXd &result)
{
  if (mtype ==-1)
  	throw std::runtime_error("Pardiso mtype not set.");
  
#ifdef PLOTS_PARDISO
  /* -------------------------------------------------------------------- */
  /* ..  pardiso_chkvec(...)                                              */
  /*     Checks the given vectors for infinite and NaN values             */
  /*     Input parameters (see PARDISO user manual for a description):    */
  /*     Use this functionality only for debugging purposes               */
  /* -------------------------------------------------------------------- */
  
  pardiso_chkvec (&numRows, &nrhs, rhs.data(), &error);
  if (error != 0) {
    printf("\nERROR  in right hand side: %d", error);
    exit(1);
  }
  
  /* -------------------------------------------------------------------- */
  /* .. pardiso_printstats(...)                                           */
  /*    prints information on the matrix to STDOUT.                       */
  /*    Use this functionality only for debugging purposes                */
  /* -------------------------------------------------------------------- */
  
  pardiso_printstats (&mtype, &numRows, a.data(), ia.data(), ja.data(), &nrhs, rhs.data(), &error);
  if (error != 0) {
    printf("\nERROR right hand side: %d", error);
    exit(1);
  }
  
#endif
  result.resize(numRows, 1);
  /* -------------------------------------------------------------------- */
  /* ..  Back substitution and iterative refinement.                      */
  /* -------------------------------------------------------------------- */
  phase = 33;
  
  iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
  
  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
           &numRows, a.data(), ia.data(), ja.data(), &idum, &nrhs,
           iparm, &msglvl, rhs.data(), result.data(), &error);
  
  if (error != 0) 
  	throw std::runtime_error(std::string("\nERROR during solution: ") + std::to_string(error));
	
#ifdef PLOTS_PARDISO
  printf("\nSolve completed ... ");
  printf("\nThe solution of the system is: ");
  for (i = 0; i < numRows; i++) {
    printf("\n x [%d] = % f", i, result.data()[i] );
  }
  printf ("\n\n");
#endif
}

template <typename vectorTypeI, typename vectorTypeS>
PardisoSolver<vectorTypeI,vectorTypeS>::~PardisoSolver()
{
  if (mtype == -1)
    return;
  /* -------------------------------------------------------------------- */
  /* ..  Termination and release of memory.                               */
  /* -------------------------------------------------------------------- */
  phase = -1;                 /* Release internal memory. */
  
  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
           &numRows, &ddum, ia.data(), ja.data(), &idum, &nrhs,
           iparm, &msglvl, &ddum, &ddum, &error);
}

template class PardisoSolver<std::vector<MKL_INT, std::allocator<MKL_INT> >, std::vector<double, std::allocator<double> > >;

template class PardisoSolver<Eigen::Matrix<MKL_INT, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >;

//extract II,JJ,SS (row,column and value vectors) from sparse matrix, Eigen version 
//Olga Diamanti's method for PARDISO
void extract_ij_from_matrix(const Eigen::SparseMatrix<double>& A, Eigen::Matrix<MKL_INT, Eigen::Dynamic, 1> & II, Eigen::Matrix<MKL_INT, Eigen::Dynamic, 1> & JJ, Eigen::VectorXd & SS)
{
	II.resize(A.nonZeros());
	JJ.resize(A.nonZeros());
	SS.resize(A.nonZeros());
	MKL_INT ind = 0;
	for (MKL_INT k = 0; k < A.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
		{
			double ss = it.value();
			MKL_INT ii = it.row();   // row index
			MKL_INT jj = it.col();   // row index
			{
				II[ind] = ii;
				JJ[ind] = jj;
				SS[ind] = ss;
				ind++;
			}
		}
}

//extract II,JJ,SS (row,column and value vectors) from sparse matrix, std::vector version
void extract_ij_from_matrix(const Eigen::SparseMatrix<double>& A, std::vector<MKL_INT>& II, std::vector<MKL_INT>& JJ, std::vector<double>& SS)
{
	II.clear();
	JJ.clear();
	SS.clear();
	for (MKL_INT k = 0; k < A.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
		{
			double ss = it.value();
			MKL_INT ii = it.row();   // row index
			MKL_INT jj = it.col();   // row index
			{
				II.push_back(ii);
				JJ.push_back(jj);
				SS.push_back(ss);
			}
		}
}
