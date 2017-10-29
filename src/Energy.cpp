#include "Energy.h"
Energy::Energy()
	:
 	//position(make_unique<Position>()),
	symDirichlet(make_unique<DistortionSymDir>())
{
}

void Energy::init(unsigned int nf, const MatX3& V, const MatX3i& F)
{
 	//position->init(F, V.rows());
	symDirichlet->init(V, F);
}

void Energy::evaluate_f(const Vec& x, double& f)
{
	map_to_X(x);

	double fd;
	Vec gd;
	SpMat hd;
	symDirichlet->value(X, fd);


// 	double fp;
// 	Vec gp;
// 	SpMat hp;
// 	position->evaluate_fgh(X, fp, gp, hp, Position::eval_mode::F);

	f = fd;
}

void Energy::evaluate_fgh(const Vec& x, double& f, Vec& g, SpMat& h)
{
	map_to_X(x);

	double fs, fd;
	Vec gs, gd;
	SpMat hd;

	// Distortion
	symDirichlet->value(X, fd);
	symDirichlet->gradient(X, gd);
	symDirichlet->hessian(X);


	// Position
// 	double fp;
// 	Vec gp;
// 	SpMat hp;
// 	position->evaluate_fgh(X, fp, gp, hp, Position::eval_mode::FGH);

	f = fd;
	g = gd;
}

inline void Energy::map_to_X(const Vec& x)
{
	X = Eigen::Map<const MatX2>(x.data(), x.rows() / 2, 2);
}