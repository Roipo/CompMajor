#include "Newton.h"
#include <igl\read_triangle_mesh.h>
#include <igl\writeOBJ.h>

using namespace std;

int main(int argc, char** argv)
{
	unsigned int num_threads = max(atoi(getenv("OMP_NUM_THREADS")), 1);
	omp_set_num_threads(num_threads);

	if (argc < 3)
	{
		cout << "Syntax: Parameterization_cmd.exe <Input file name> <Output file name> <Number of iterations>";
		return false;
	}
	MatX3 V;
	MatX3i F;
	if (!igl::read_triangle_mesh(argv[1], V, F))
	{
		cerr << "Failed to load mesh: " << argv[1] << endl;
		return false;
	}
	
	Newton solver;
	solver.init(V, F);
	for (int i = 0; i < atoi(argv[3]); i++)
	{
		solver.step();
		solver.linesearch();
	}

	igl::writeOBJ(argv[2], V, F, MatX3(0,0), MatX3i(0,0), solver.uv, F);
	return 0;
}