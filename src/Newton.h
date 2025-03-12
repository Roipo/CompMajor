#pragma once

#include "Solver.h"
#include "EigenTypes.h"
#include "PardisoSolver.h"

#include <iostream>
#include <Eigen/SparseCholesky>

using namespace std;

class Newton : public Solver
{
public:
	Newton();

	int step();
	void linesearch();
	bool test_progress();
	void internal_init();
	void internal_update_external_mesh();

private:
	// Wrapper function for flip_avoiding_line_search
	double eval_ls(Mat& x);

	// norm of the progress on the mesh
	double diff_norm=0.0;

	PardisoSolver<vector<MKL_INT>, vector<double>> pardiso;
};