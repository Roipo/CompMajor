#pragma once

#include "Solver.h"
#include "EigenTypes.h"
#include "PardisoSolver.h"

#include <iostream>
#include <Eigen/SparseCholesky>

using namespace std;

struct StepTiming {
	double eval_value_time;     // Time for computing function value
	double eval_gradient_time;  // Time for computing gradient
	double eval_hessian_time;   // Time for computing hessian
	double matrix_prep_time;
	double factorization_time;
	double solve_time;
	double analyze_pattern_time;  // One-time timing for pattern analysis
};

class Newton : public Solver
{
public:
	Newton();

	int step();
	void linesearch();
	void initialize();
	void internal_update_external_mesh();
	const StepTiming& get_last_step_timing() const { return last_step_timing; }

private:
	// Wrapper function for flip_avoiding_line_search
	double eval_ls(Mat& x);

	// norm of the progress on the mesh
	double diff_norm=0.0;

	PardisoSolver<vector<MKL_INT>, vector<double>> pardiso;
	
	// Store timing data for the last step
	StepTiming last_step_timing;
};