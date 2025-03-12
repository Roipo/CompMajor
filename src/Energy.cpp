#include "Energy.h"
Energy::Energy()
{

}

void Energy::init(unsigned int nf, const MatX3& V, const MatX3i& F)
{
	symDirichlet.init(V, F);
}

void Energy::evaluate_f(const Vec& x, double& f)
{
	map_to_X(x);

	double fd;
	Vec gd;
	SpMat hd;
	symDirichlet.value(X, fd);

	f = fd;
}

void Energy::evaluate_fgh(const Vec& x, double& f, Vec& g, SpMat& h)
{
	int n = x.rows() / 2;
	X = MatX2::Zero(n, 2);
	X = Eigen::Map<const MatX2>(x.data(), x.rows() / 2, 2);

	double fs, fd;
	Vec gs, gd;
	SpMat hd;

	symDirichlet.value(X, fd);
	symDirichlet.gradient(X, gd);
	symDirichlet.hessian(X);

	f = fd;
	g = gd;
}

inline void Energy::map_to_X(const Vec& x)
{
	X = Eigen::Map<const MatX2>(x.data(), x.rows() / 2, 2);
}