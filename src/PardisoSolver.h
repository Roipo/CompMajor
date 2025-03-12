//
//  PardisoSolver.h
//
//  Created by Olga Diamanti on 07/01/15.
//  Copyright (c) 2015 Olga Diamanti. All rights reserved.
//

#ifndef _PardisoSolver__
#define _PardisoSolver__

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <mkl_pardiso.h>

//extract II,JJ,SS (row,column and value vectors) from sparse matrix, Eigen version 
//Olga Diamanti's method for PARDISO
void extract_ij_from_matrix(const Eigen::SparseMatrix<double> &A,
	Eigen::VectorX<MKL_INT> &II,
	Eigen::VectorX<MKL_INT> &JJ,
	Eigen::VectorXd &SS);

//extract II,JJ,SS (row,column and value vectors) from sparse matrix, std::vector version
void extract_ij_from_matrix(const Eigen::SparseMatrix<double> &A,
	std::vector<MKL_INT> &II,
	std::vector<MKL_INT> &JJ,
	std::vector<double> &SS);


template <typename vectorTypeI, typename vectorTypeS>
 class PardisoSolver
 {
 public:
   
   PardisoSolver() ;
   ~PardisoSolver();
   
   void set_type(int _mtype, bool is_upper_half = false);
   
   void init();

   void set_pattern(const vectorTypeI &II,
                    const vectorTypeI &JJ,
                    const vectorTypeS &SS);
   void analyze_pattern();
   
   bool factorize();
   
   void solve(Eigen::VectorXd &rhs,
              Eigen::VectorXd &result);

   void update_a(const vectorTypeS &SS);

 protected:
   //vector that indicates which of the elements II,JJ input will be
   //kept and read into the matrix (for symmetric matrices, only those
   //elements II[i],JJ[i] for which II[i]<<JJ[i] will be kept)
   std::vector<MKL_INT> lower_triangular_ind;

   Eigen::Matrix<MKL_INT, Eigen::Dynamic, 1> ia, ja;
   std::vector<Eigen::Matrix<MKL_INT, Eigen::Dynamic, 1>> iis;

   Eigen::VectorXd a;
   MKL_INT numRows;

   //pardiso stuff
   /*
    1: real and structurally symmetric, supernode pivoting
    2: real and symmetric positive definite
    -2: real and symmetric indefinite, diagonal or Bunch-Kaufman pivoting
    11: real and nonsymmetric, complete supernode pivoting
    */

   // Remember if matrix is symmetric or not, to
   // decide whether to eliminate the non-upper-
   // diagonal entries from the input II,JJ,SS
   bool is_symmetric;
   bool is_upper_half;
   
   void* pt[64];
   MKL_INT mtype;
   MKL_INT nrhs = 1;
   MKL_INT iparm[64];
   double dparm[64];
   MKL_INT maxfct, mnum, phase, error, msglvl, solver = 0;
   MKL_INT num_procs;
   MKL_INT idum = 0;
   double ddum = 0;
 };


#endif /* defined(_PardisoSolver__) */
