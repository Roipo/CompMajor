#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include "SolverPlugin.h"
#include "Energy.h"

using namespace std;
using namespace igl::opengl::glfw;

int main(int argc, char** argv)
{
	//unsigned int num_threads = max(atoi(getenv("OMP_NUM_THREADS")), 1);
	//omp_set_num_threads(num_threads);

	Viewer viewer;

	SolverPlugin solverPlugin;
//	viewer.core.rotation_type = viewer::ViewerCore::ROTATION_TYPE_TRACKBALL;
	viewer.plugins.push_back(&solverPlugin);

	// start viewer
	viewer.launch(true, false,1920,1080);
	return 0;
}