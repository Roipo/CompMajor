#include "Solver.h"

#include "Utils.h"
#include "Energy.h"

#include <iostream>
#include <igl/readOBJ.h>
#include <igl/file_dialog_open.h>
#include <igl/flipped_triangles.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>

Solver::Solver()
	:
	energy(make_shared<Energy>())
{
}

void Solver::init(const MatX3& V, const MatX3i& F)
{
	using namespace placeholders;
	this->F = F;
	computeTutte(V, F, uv);
	energy->init(F.rows(), V, F);
	m_x = Eigen::Map<Vec>(uv.data(), uv.rows() * uv.cols());
	eval_f = [this](const Vec& x, double& f) { energy->evaluate_f(x, f); };
	eval_fgh = [this](const Vec& x, double& f, Vec& g, SpMat& h) { energy->evaluate_fgh(x, f, g, h); };
	initialize();

	cout << "Solver Initialization Done: Energy: "<< f << endl;
}

void Solver::computeTutte(const MatX3& V, const MatX3i& F, MatX2 &uv_init)
{
	
	Eigen::VectorXi bnd; Eigen::MatrixXd bnd_uv;
	igl::boundary_loop(F,bnd);
	// igl::map_vertices_to_circle(V,bnd,bnd_uv);
	map_vertices_to_circle_area_normalize(V, F, bnd, bnd_uv);

	igl::harmonic(V,F,bnd,bnd_uv,1,uv_init);
	if (igl::flipped_triangles(uv_init,F).size() != 0) 
		igl::harmonic(F,bnd,bnd_uv,1,uv_init); // use uniform laplacian

}

void Solver::map_vertices_to_circle_area_normalize(
	const MatX3& V,
	const MatX3i& F,
	const Eigen::VectorXi& bnd,
	Eigen::MatrixXd& UV ) {

	double radius = sqrt(1.0 / (M_PI));	

	// Get sorted list of boundary vertices
	std::vector<int> interior,map_ij;
	map_ij.resize(V.rows());
	interior.reserve(V.rows()-bnd.size());
  
	std::vector<bool> isOnBnd(V.rows(),false);
	for (int i = 0; i < bnd.size(); i++)
	{
	  isOnBnd[bnd[i]] = true;
	  map_ij[bnd[i]] = i;
	}
  
	for (int i = 0; i < (int)isOnBnd.size(); i++)
	{
	  if (!isOnBnd[i])
	  {
		map_ij[i] = interior.size();
		interior.push_back(i);
	  }
	}
  
	// Map boundary to unit circle
	std::vector<double> len(bnd.size());
	len[0] = 0.;
  
	for (int i = 1; i < bnd.size(); i++)
	{
	  len[i] = len[i-1] + (V.row(bnd[i-1]) - V.row(bnd[i])).norm();
	}
	double total_len = len[len.size()-1] + (V.row(bnd[0]) - V.row(bnd[bnd.size()-1])).norm();
  
	UV.resize(bnd.size(),2);
	for (int i = 0; i < bnd.size(); i++)
	{
	  double frac = len[i] * (2. * M_PI) / total_len;
	  UV.row(map_ij[bnd[i]]) << radius*cos(frac), radius*sin(frac);
	}

}