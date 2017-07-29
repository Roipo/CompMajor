#pragma once

#ifndef SOLVERPLUGIN_H
#define SOLVERPLUGIN_H

#include "EigenTypes.h"
#include "SolverWrapper.h"

#include <igl\viewer\Viewer.h>
#include <igl\viewer\ViewerPlugin.h>
#include <thread>
#include <nanogui\slider.h>

#ifndef INT_INF
#define INT_INF numeric_limits<int>::max()
#endif

using namespace std;
using namespace nanogui;

class SolverPlugin : public igl::viewer::ViewerPlugin
{
public:
	SolverPlugin();
	void init(igl::viewer::Viewer *viewer);

	bool load(string filename);
	void add_color_clamp_slider(const string& name, const shared_ptr<double>& max_value, const shared_ptr<double>& value);
	void export_uv_to_obj();
	void add_texture_slider(nanogui::Window* window, double& var, const string& name);
	void initialize();
	
	bool mouse_move(int mouse_x, int mouse_y);
	bool process_mouse_move();

	bool mouse_down(int button, int modifier);
	bool mouse_up(int button, int modifier);
	bool mouse_scroll(float delta_y);

	void rotate(double phi_x, double phi_y);
	inline string removeTrailingZeros(string& s);
	void translate(double offset_x, double offset_y);
	void translate_uv_mesh(double offset_x, double offset_y);
	void translate_triangle(double offset_x, double offset_y);
	
	bool pre_draw();
	bool key_down(int key, int modifiers);
	bool key_up(int key, int modifiers);

	void update_mesh();
	void start_solver_thread();
	void stop_solver_thread();

	bool rotation_enabled = false;
	bool translation_enabled = false;
	bool uv_translation_enabled = false;

	int last_mouse_x, last_mouse_y;

	thread solver_thread;
	unique_ptr<SolverWrapper> solver_wrapper;

	// lambda slider
	Slider *slider;

	unsigned int uv_id = 0, mesh_id = 0;

	double texture_size = 0.5;
	double max_texture_val = 10.;

	int solver_iter_counter = 0;
	bool use_num_steps = false;

	int resolution_factor = 4;

	bool update_colors = false;

	bool store_3d_mesh = true;

	shared_ptr<double> dist_color_clamp = make_shared<double>(0.5);

	shared_ptr<double> max_dist_color_value = make_shared<double>(1.0);

	string mesh_filename;

	bool load_uv_from_file = false;

	Vec3 xh_center_3d;

private:
	// Pointer to the nano gui
	nanogui::FormHelper* bar;
	nanogui::Window* orig_window;
	// The 3d mesh
	MatX3 V;
	MatX3i F;

	// Rotation matrices
	Mat3 Rx, Ry;

	MatX3 uv_triangle_colors, mesh_triangle_colors;
	int hit_triangle = -1, last_hit_triangle = -1, hovered_triangle = -1, hovered_vertex = -1;
	//RVec3 last_hit_triangle_color;
	
	MatX3 mesh_pos_down, uv_mesh_pos_down;
	RVec3 last_mouse_pos, projected_mouse_down, point_pos_down;
	RVec3 point_pos;
	bool move_point = false;
	MatX3 mesh_3d_normals_down;

	bool mesh_loaded = false;
	bool mouse_on_uv_side = false;
	bool mouse_over_3d_mesh = false;
	bool mouse_over_2d_mesh = false;
	bool mouse_over_uv_mesh = false;
	bool release_when_translation_done = false;

	// minimalistic menu
	nanogui::Window* param_menu;
	nanogui::Window* bar_window;

	// coords used in process_mouse_move
	int curr_mouse_x, curr_mouse_y;

	// some mouse states
	bool MOUSE_LEFT = false;
	bool MOUSE_MID = false;
	bool MOUSE_RIGHT = false;

	bool show_distortion_error = false;
	bool colorByRGB=false;
	MatX3 RGBColors;

	//colors
	RVec3 C, C_hover, white, red, C_merge, zero3, ones3, black;

};

#endif