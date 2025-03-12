#pragma once

#include "Energy.h"
#include "EigenTypes.h"
#include <mkl_types.h>

#include <atomic>
#include <functional>
#include <shared_mutex>
using namespace std;

class Solver
{
public:
	Solver();

	void init(const MatX3& V, const MatX3i& F);

	void computeTutte(const MatX3& V, const MatX3i& F, MatX2 &uv_init);
	// Pointer to the energy class
	shared_ptr<Energy> energy;

	MatX2 uv;
	MatX3i F;

	// External (interface) and internal working mesh
	Vec ext_x, m_x;

protected:
	// Descent direction evaluated in step
	Vec p;

	// Function pointers to the full and value-only energy evaluation
	function<void(const Vec&, double&)> eval_f;
	function<void(const Vec&, double&, Vec&, SpMat&)> eval_fgh;
	
	// Current energy, gradient and hessian
	double f;
	Vec g;
	SpMat h;
	
	// pardiso variables
	vector<MKL_INT> IId, JJd, IIp, JJp, II, JJ;
	vector<double> SSd, SSp, SS;

private:
	virtual int step() = 0;
	virtual void linesearch() = 0;
	virtual void internal_init() = 0;
};