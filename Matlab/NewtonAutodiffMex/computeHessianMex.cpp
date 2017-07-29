#include <mex.h>
#include <map>
#include <vector>
#include <string>
#include <Eigen\Core>
#include <Eigen\Sparse>

#include "autodiff.h"
#include <igl\cotmatrix_entries.h> 

using std::map;
using std::vector;
using std::string;




DECLARE_DIFFSCALAR_BASE();


typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
typedef Eigen::VectorXd                     Gradient;
typedef Eigen::MatrixXd                     Hessian;
typedef DScalar2<double, Gradient, Hessian> DScalar;


Eigen::MatrixXd	 V, F;
std::vector<double> m_cached_edges_1;
std::vector<double> m_cached_edges_2;
std::vector<double> m_cached_dot_prod;
Eigen::MatrixXd m_cot_entries;
Eigen::VectorXd m_dblArea_orig;
Eigen::VectorXd m_dbl_sqrt_area;

DScalar compute_face_energy_left_part(const const MapMat& uv, int f);
DScalar compute_face_energy_right_part(const const MapMat& uv, int f_idx, double orig_t_dbl_area);

 // Treat everything as a double
mxArray* EigenSparse2Matlab(const Eigen::SparseMatrix<double>& M)
{
	using namespace std;
	const int m = M.rows();
	const int n = M.cols();
	// THIS WILL NOT WORK FOR ROW-MAJOR
	assert(n == M.outerSize());
	const int nzmax = M.nonZeros();
	mxArray * mx_data = mxCreateSparse(m, n, nzmax, mxREAL);
	// Copy data immediately
	double * pr = mxGetPr(mx_data);
	mwIndex * ir = mxGetIr(mx_data);
	mwIndex * jc = mxGetJc(mx_data);

	// Iterate over outside
	int k = 0;
	for (int j = 0; j < M.outerSize(); j++)
	{
		jc[j] = k;
		// Iterate over inside
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, j); it; ++it)
		{
			pr[k] = it.value();
			ir[k] = it.row();
			k++;
		}
	}
	jc[M.outerSize()] = k;

	return mx_data;
}
double compute_energy_gradient_hessian(const const MapMat& uv, Eigen::VectorXd& grad, Eigen::SparseMatrix<double>& hessian, std::vector<Eigen::MatrixXd> &allHessians, bool bVanilla = false) {

	// can save some computation time
	hessian.resize(2 * V.rows(), 2 * V.rows());
	hessian.reserve(10 * 2 * V.rows());
	std::vector<Eigen::Triplet<double> > IJV;//(10*2*6*V.rows());
	IJV.reserve(36 * F.rows());
	grad.resize(2 * V.rows()); grad.setZero();

	double energy = 0;
	for (int i = 0; i < F.rows(); i++) {
		DiffScalarBase::setVariableCount(6); // 3 vertices with 2 rows for each
		DScalar l_part = compute_face_energy_left_part(uv, i);
		DScalar r_part = compute_face_energy_right_part(uv, i, m_dblArea_orig(i));

		DScalar temp = 0.5*l_part * r_part;
		energy += temp.getValue();

		Eigen::VectorXd local_grad = temp.getGradient();
		for (int v_i = 0; v_i < 3; v_i++) {
			int v_global = F(i, v_i);

			grad(v_global) = grad(v_global) + local_grad(v_i * 2); // x
			grad(v_global + V.rows()) = grad(v_global + V.rows()) + local_grad(v_i * 2 + 1); // y
		}

		Eigen::MatrixXd local_hessian = temp.getHessian();
		allHessians.push_back(local_hessian);
		
		if (!bVanilla) {
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> es(local_hessian);
			Eigen::VectorXd D = es.eigenvalues();
			Eigen::MatrixXd U = es.eigenvectors();
			D = D.unaryExpr([](double x) {return (x < 0) ? 0 : x; });
			local_hessian = U * D.asDiagonal()* U.transpose();
		}

		for (int v1 = 0; v1 < 6; v1++) {
			for (int v2 = 0; v2 < 6; v2++) {
				int v1_global = F(i, v1 / 2) + (v1 % 2)*V.rows();
				int v2_global = F(i, v2 / 2) + (v2 % 2)*V.rows();

				IJV.push_back(Eigen::Triplet<double>(v1_global, v2_global, local_hessian(v1, v2)));
			}
		}
	}
	hessian.setFromTriplets(IJV.begin(), IJV.end());
	return energy;
}
DScalar compute_face_energy_left_part(const const MapMat& uv, int f) {

	int v_1 = F(f, 0); int v_2 = F(f, 1); int v_3 = F(f, 2);

	// compute current triangle squared area
	DScalar x1(0 * 2, uv(v_1, 0)); DScalar y1(0 * 2 + 1, uv(v_1, 1)); //DScalar x0(F(f,0),0); DScalar y0(F(f,0),1);
	DScalar x2(1 * 2, uv(v_2, 0)); DScalar y2(1 * 2 + 1, uv(v_2, 1)); //DScalar x0(F(f,0),0); DScalar y0(F(f,0),1);
	DScalar x3(2 * 2, uv(v_3, 0)); DScalar y3(2 * 2 + 1, uv(v_3, 1)); //DScalar x0(F(f,0),0); DScalar y0(F(f,0),1);


	DScalar rx = x1 - x3;//uv(F(f,0),0)-uv(F(f,2),0);
	DScalar sx = x2 - x3;//uv(F(f,1),0)-uv(F(f,2),0);
	DScalar ry = y1 - y3;//uv(F(f,0),1)-uv(F(f,2),1);
	DScalar sy = y2 - y3;//uv(F(f,1),1)-uv(F(f,2),1);
	DScalar dblAd = rx*sy - ry*sx;
	DScalar uv_sqrt_dbl_area = dblAd*dblAd;


	return (1 + (m_dbl_sqrt_area(f) / uv_sqrt_dbl_area));
}
DScalar compute_face_energy_right_part(const const MapMat& uv, int f_idx, double orig_t_dbl_area) {
	int v_1 = F(f_idx, 0); int v_2 = F(f_idx, 1); int v_3 = F(f_idx, 2);

	DScalar x1(0 * 2, uv(v_1, 0)); DScalar y1(0 * 2 + 1, uv(v_1, 1)); //DScalar x0(F(f,0),0); DScalar y0(F(f,0),1);
	DScalar x2(1 * 2, uv(v_2, 0)); DScalar y2(1 * 2 + 1, uv(v_2, 1)); //DScalar x0(F(f,0),0); DScalar y0(F(f,0),1);
	DScalar x3(2 * 2, uv(v_3, 0)); DScalar y3(2 * 2 + 1, uv(v_3, 1)); //DScalar x0(F(f,0),0); DScalar y0(F(f,0),1);

																	  //DScalar part_1 = (uv.row(v_3)-uv.row(v_1)).squaredNorm() * m_cached_edges_1[f_idx];
	DScalar part_1 = (pow(x3 - x1, 2) + pow(y3 - y1, 2)) * m_cached_edges_1[f_idx];
	//part_1 += (uv.row(v_2)-uv.row(v_1)).squaredNorm()* m_cached_edges_2[f_idx];
	part_1 += (pow(x2 - x1, 2) + pow(y2 - y1, 2)) * m_cached_edges_2[f_idx];
	part_1 /= (2 * orig_t_dbl_area);

	//DScalar part_2_1 = (uv.row(v_3)-uv.row(v_1)).dot(uv.row(v_2)-uv.row(v_1));
	DScalar part_2_1 = (x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1);
	double part_2_2 = m_cached_dot_prod[f_idx];
	DScalar part_2 = -(part_2_1 * part_2_2) / (orig_t_dbl_area);

	return part_1 + part_2;
}

void precompute()
{
	// igl returns half a cot, so we need to double it by 2
	igl::cotmatrix_entries(V, F, m_cot_entries);
	m_cot_entries = m_cot_entries * 2;

	igl::doublearea(V, F, m_dblArea_orig);
	m_dbl_sqrt_area.resize(m_dblArea_orig.rows(), 1);
	for (int i = 0; i < m_dblArea_orig.rows(); i++) {
		m_dbl_sqrt_area(i) = pow(m_dblArea_orig(i), 2);
	}

	//m_cached_l_energy_per_face.resize(F.rows());
	//m_cached_r_energy_per_face.resize(F.rows());

	m_cached_edges_1.resize(F.rows());
	m_cached_edges_2.resize(F.rows());
	m_cached_dot_prod.resize(F.rows());
	for (int f = 0; f < F.rows(); f++) {
		int v_1 = F(f, 0); int v_2 = F(f, 1); int v_3 = F(f, 2);

		m_cached_edges_1[f] = (V.row(v_2) - V.row(v_1)).squaredNorm();
		m_cached_edges_2[f] = (V.row(v_3) - V.row(v_1)).squaredNorm();
		m_cached_dot_prod[f] = (V.row(v_3) - V.row(v_1)).dot(V.row(v_2) - V.row(v_1));
	}
}
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
		within an if statement, because it will never get to the else
		statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
		the MEX-file) */

	if (!mxIsChar(prhs[0]))
	{
		mexErrMsgIdAndTxt("MATLAB:inputNotChar", "First input must be a string.");
	}

	/* first input must be a double array */
	std::string input_buf = mxArrayToString(prhs[0]);
	if (input_buf == "Init")
	{
		if (nrhs != 3)
		{
			mexErrMsgIdAndTxt("MATLAB:invalidNumInputs", "Three inputs required.");
		}
		if (nlhs != 0)
		{
			mexErrMsgIdAndTxt("MATLAB:maxlhs", "No output exptected.");
		}
		if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
		{
			mexErrMsgIdAndTxt("MATLAB:inputNotDouble", " Second and third Inputs must be double matrices.");
		}

		V =  MapMat(mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]));
		F =  MapMat(mxGetPr(prhs[2]), mxGetM(prhs[2]), mxGetN(prhs[2])).array()-1;

		precompute();

	}
	else if(input_buf == "Compute")
	{
		bool bVanilla = false;
		if (nrhs < 2)
		{
			mexErrMsgIdAndTxt("MATLAB:invalidNumInputs", "Two inputs required.");
		}
		else {
			if (nrhs == 3)
				bVanilla = *mxGetPr(prhs[2]);
		}
		MapMat uv(mxGetPr(prhs[1]), mxGetM(prhs[1])/2, 2);

		Eigen::VectorXd grad;
		Eigen::SparseMatrix<double> hessian;
		std::vector<Eigen::MatrixXd> allHessians;
		compute_energy_gradient_hessian(uv, grad, hessian, allHessians,bVanilla);
		plhs[0] = EigenSparse2Matlab(hessian);
		if (nlhs == 2)
		{
			size_t dims[] = { 6,6,F.rows() };
			plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
			double *ptr = (double*)mxGetData(plhs[1]);
			for (int i = 0; i < F.rows(); i++)
			{
				memcpy(ptr+36*i, allHessians[i].data(), 36 * mxGetElementSize(plhs[1]));
			}
		}
	}
	else
		mexErrMsgIdAndTxt("MATLAB:inputNotDouble", "First parameter must be \"Init\" or \"Compute\".");
	/* Create N by 1 cell array such that N is the number of boundary's on the mesh */     
	return;
}