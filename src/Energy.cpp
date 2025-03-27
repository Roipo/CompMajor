#include "Energy.h"
#include <chrono>

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

	// Time value computation
	auto value_start = std::chrono::high_resolution_clock::now();
	symDirichlet.value(X, fd);
	auto value_end = std::chrono::high_resolution_clock::now();
	last_eval_timing.value_time = std::chrono::duration<double>(value_end - value_start).count();

	// Time gradient computation
	auto gradient_start = std::chrono::high_resolution_clock::now();
	symDirichlet.gradient(X, gd);
	auto gradient_end = std::chrono::high_resolution_clock::now();
	last_eval_timing.gradient_time = std::chrono::duration<double>(gradient_end - gradient_start).count();

	// Time hessian computation
	auto hessian_start = std::chrono::high_resolution_clock::now();
	symDirichlet.hessian(X);
	auto hessian_end = std::chrono::high_resolution_clock::now();
	last_eval_timing.hessian_time = std::chrono::duration<double>(hessian_end - hessian_start).count();

	f = fd;
	g = gd;
}

inline void Energy::map_to_X(const Vec& x)
{
	X = Eigen::Map<const MatX2>(x.data(), x.rows() / 2, 2);
}