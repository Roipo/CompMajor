#include <igl\viewer\Viewer.h>
#include "SolverPlugin.h"

using namespace std;
using namespace igl::viewer;

int main(int argc, char** argv)
{
	unsigned int num_threads = max(atoi(getenv("OMP_NUM_THREADS")), 1);
	omp_set_num_threads(num_threads);

	Viewer viewer;

	SolverPlugin plugin;
	viewer.core.rotation_type = igl::viewer::ViewerCore::ROTATION_TYPE_TRACKBALL;
	viewer.plugins.push_back(&plugin);
	viewer.core.background_color << 1., 1., 1., 1.; // 0.25, 0.25, 0.25, 1.0;

	// set to ortographic view to synchronize pointer with models
	viewer.core.orthographic = true;

	// start viewer
	viewer.launch(true, false, 2400, 1350);

	return 0;
}