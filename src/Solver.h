#pragma once

#include "Energy.h"
#include <mkl_types.h>

#include <atomic>
#include <functional>
#include <shared_mutex>

class Solver
{
public:
	Solver();

	void init(const MatX3& V, const MatX3i& F);

	void computeTutte(const MatX3& V, const MatX3i& F, MatX2 &uv_init);
	// Pointer to the energy class
	shared_ptr<Energy> energy;

	Eigen::MatrixX2d uv;
	Eigen::MatrixX3i F;

	// External (interface) and internal working mesh
	Eigen::VectorXd m_x;

	// Descent direction evaluated in step
	Vec p;

	// Function pointers to the full and value-only energy evaluation
	function<void(const Eigen::VectorXd&, double&)> eval_f;
	function<void(const Eigen::VectorXd&, double&, Eigen::VectorXd&, SpMat&)> eval_fgh;
	
	// Current energy, gradient and hessian
	double f;
	Eigen::VectorXd g;
	Eigen::SparseMatrix<double> h;
	
	// pardiso variables
	vector<MKL_INT> IId, JJd, II, JJ;
	vector<double> SSd, SS;

private:
	virtual int step() = 0;
	virtual void linesearch() = 0;
	virtual void initialize() = 0;
};