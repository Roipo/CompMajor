#include "Newton.h"

#include <chrono>
#include <igl/flip_avoiding_line_search.h>
#include <mkl_types.h>

Newton::Newton() {
	pardiso.set_type(2, true);
}

int Newton::step()
{
	eval_fgh(m_x, f, g, h);

 	SSd = energy->symDirichlet.SS;

	SS.clear();
	SS.insert(SS.end(), SSd.begin(), SSd.end());
	pardiso.update_a(SS);
	try
	{
		pardiso.factorize();
	}
	catch (runtime_error& err)
	{
		cout << err.what();
		return -1;
	}
	Vec rhs = -g;
	pardiso.solve(rhs, p);
	return 0;
}

bool Newton::test_progress()
{
	return true;
}

void Newton::internal_init()
{
	eval_fgh(m_x, f, g, h);
	
	
	IId = energy->symDirichlet.II;
	JJd = energy->symDirichlet.JJ;
	SSd = energy->symDirichlet.SS;

    // convert to eigen triplets
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(IId.size());
    for (size_t i = 0; i < IId.size(); ++i) {
        triplets.emplace_back(IId[i], JJd[i], SSd[i]);
    }
    
    // Create sparse matrix (assuming you want a square matrix)
    MKL_INT n = *std::max_element(IId.begin(), IId.end()) + 1;
    Eigen::SparseMatrix<double> sparse_matrix(n, n);
    sparse_matrix.setFromTriplets(triplets.begin(), triplets.end());

	// find pattern for initialization
	II.insert(II.end(), IId.begin(), IId.end());
	JJ.insert(JJ.end(), JJd.begin(), JJd.end());
	SS.insert(SS.end(), SSd.begin(), SSd.end());

	pardiso.set_pattern(II, JJ, SS);
	pardiso.analyze_pattern();
}

void Newton::internal_update_external_mesh()
{
	diff_norm = (ext_x - m_x).norm();
	ext_x = m_x;
}

void Newton::linesearch()
{
	Mat m_x2 = Eigen::Map<MatX2>(m_x.data(), m_x.rows() / 2, 2);
	Mat p2 = Eigen::Map<const MatX2>(p.data(), p.rows() / 2, 2);
	Mat m_plus_p = m_x2 + p2;
	double alpha = igl::flip_avoiding_line_search(F, m_x2, m_plus_p, bind(&Newton::eval_ls, this, placeholders::_1));
	m_x = Eigen::Map<Vec>(m_x2.data(), m_x2.rows() * m_x2.cols());
}

double Newton::eval_ls(Mat& x)
{
	double f;
	Vec g;
	SpMat h;
	Vec vec_x = Eigen::Map<Vec>(x.data(), x.rows()  * x.cols(), 1);
	eval_f(vec_x, f);
	return f;
}