#include "Newton.h"
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
	// unsigned int num_threads = max(atoi(getenv("OMP_NUM_THREADS")), 1);
	// omp_set_num_threads(num_threads);

	if (argc < 3)
	{
		cout << "Syntax: Parameterization_cmd.exe <Input file name> <Output file name> <Number of iterations> [save_intermediate_meshes]\n"
			 << "  save_intermediate_meshes: 1 to save meshes for each iteration, 0 or omitted to save only final mesh";
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
	
	// Initialize timing variables
	auto total_start = std::chrono::high_resolution_clock::now();
	double total_step_time = 0.0;
	double total_linesearch_time = 0.0;
	
	// Create CSV file
	std::string output_filename = std::string(argv[2]) + "_timing.csv";
	std::ofstream file(output_filename);
	if (!file.is_open()) {
		cerr << "Error: Could not create timing output file: " << output_filename << endl;
		return false;
	}
	cout << "Writing timing data to: " << output_filename << endl;
	
	// Write CSV header
	file << "iteration,step_time,eval_value_time,eval_gradient_time,eval_hessian_time,matrix_prep_time,factorization_time,solve_time,linesearch_time,objective_value,gradient_norm,analyze_pattern_time\n";
	if (file.fail()) {
		cerr << "Error: Failed to write CSV header" << endl;
		file.close();
		return false;
	}
	
	// Check if we should save intermediate meshes
	bool save_intermediate = (argc > 4 && atoi(argv[4]) == 1);
	if (save_intermediate) {
		cout << "Will save intermediate meshes for each iteration" << endl;
	}
	
	for (int i = 0; i < atoi(argv[3]); i++)
	{
		cout << "Iteration " << i << endl;
		
		// Time the step() function
		auto step_start = std::chrono::high_resolution_clock::now();
		int res = solver.step();
		if (res != 0) return res;
		auto step_end = std::chrono::high_resolution_clock::now();
		double step_duration = std::chrono::duration<double>(step_end - step_start).count();
		total_step_time += step_duration;
		
		// Time the linesearch() function
		auto linesearch_start = std::chrono::high_resolution_clock::now();
		solver.linesearch();
		auto linesearch_end = std::chrono::high_resolution_clock::now();
		double linesearch_duration = std::chrono::duration<double>(linesearch_end - linesearch_start).count();
		total_linesearch_time += linesearch_duration;
		
		// Get internal step timing
		const auto& step_timing = solver.get_last_step_timing();
		
		// Write iteration data to CSV
		file << i << "," 
			 << std::fixed << std::setprecision(6) 
			 << step_duration << ","
			 << step_timing.eval_value_time << ","
			 << step_timing.eval_gradient_time << ","
			 << step_timing.eval_hessian_time << ","
			 << step_timing.matrix_prep_time << ","
			 << step_timing.factorization_time << ","
			 << step_timing.solve_time << ","
			 << linesearch_duration << ","
			 << solver.f << ","  // objective value
			 << solver.g.norm() << ","  // gradient norm
			 << step_timing.analyze_pattern_time << "\n";  // analyze_pattern time (will be 0 except for first iteration)
			 
		if (file.fail()) {
			cerr << "Error: Failed to write iteration " << i << " data" << endl;
			file.close();
			return false;
		}

		// Save mesh for this iteration if requested
		if (save_intermediate) {
			Mat uv = Eigen::Map<MatX2>(solver.m_x.data(), solver.m_x.rows() / 2, 2);
			std::stringstream ss;
			ss << argv[2] << "_iter" << std::setw(4) << std::setfill('0') << i << ".obj";
			if (!igl::writeOBJ(ss.str(), V, F, MatX3(), MatX3i(), uv, F)) {
				cerr << "Error: Failed to write mesh for iteration " << i << endl;
				file.close();
				return false;
			}
			cout << "Saved mesh to: " << ss.str() << endl;
		}
	}
	
	auto total_end = std::chrono::high_resolution_clock::now();
	double total_time = std::chrono::duration<double>(total_end - total_start).count();
	
	// Write summary statistics
	file << "\nSummary Statistics\n";
	file << "total_time," << total_time << "\n";
	file << "average_step_time," << total_step_time / atoi(argv[3]) << "\n";
	file << "average_linesearch_time," << total_linesearch_time / atoi(argv[3]) << "\n";
	file << "percentage_in_step," << (total_step_time / total_time) * 100 << "\n";
	file << "percentage_in_linesearch," << (total_linesearch_time / total_time) * 100 << "\n";
	file << "final_objective_value," << solver.f << "\n";
	file << "final_gradient_norm," << solver.g.norm() << "\n";
	file << "analyze_pattern_time," << solver.get_last_step_timing().analyze_pattern_time << "\n";
	
	if (file.fail()) {
		cerr << "Error: Failed to write summary statistics" << endl;
		file.close();
		return false;
	}
	
	file.close();
	if (file.fail()) {
		cerr << "Error: Failed to close timing file" << endl;
		return false;
	}
	cout << "Timing data successfully written to: " << output_filename << endl;
	
	// Save final mesh
	Mat uv = Eigen::Map<MatX2>(solver.m_x.data(), solver.m_x.rows() / 2, 2);
	if (!igl::writeOBJ(argv[2], V, F, MatX3(), MatX3i(), uv, F)) {
		cerr << "Error: Failed to write final mesh" << endl;
		return false;
	}
	cout << "Saved final mesh to: " << argv[2] << endl;
	
	return 0;
}
