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
	ext_x = m_x;
	eval_f = bind(&Energy::evaluate_f, energy, _1, _2);
	eval_fgh = bind(&Energy::evaluate_fgh, energy, _1, _2, _3, _4);
	internal_init();
}

void Solver::computeTutte(const MatX3& V, const MatX3i& F, MatX2 &uv_init)
{
	
	Eigen::VectorXi bnd; Eigen::MatrixXd bnd_uv;
	igl::boundary_loop(F,bnd);
	igl::map_vertices_to_circle(V,bnd,bnd_uv);

	igl::harmonic(V,F,bnd,bnd_uv,1,uv_init);
	if (igl::flipped_triangles(uv_init,F).size() != 0) 
		igl::harmonic(F,bnd,bnd_uv,1,uv_init); // use uniform laplacian

}