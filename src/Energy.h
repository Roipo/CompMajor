#pragma once

#include "EigenTypes.h"
#include "EnergySymDir.h"


#include <memory>

using namespace std;

class Energy
{
public:
	Energy();

	void init(unsigned int nf, const MatX3& V, const MatX3i& F);
	void evaluate_f(const Vec& x, double& f);
	void evaluate_fgh(const Vec& x, double& f, Vec& g, SpMat& h);

	// helper functions
	inline void map_to_X(const Vec& x);
    double w=1.0;
	DistortionSymDir symDirichlet;
	// Internally work with two-column matrix instead
	// of a vector, which is used in the solver
	MatX2 X;
};