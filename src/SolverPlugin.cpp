#include "SolverPlugin.h"

#include "svg_exporter.h"

#include <queue>
#include <deque>
#include <array>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <igl\writeOBJ.h>
#include <igl\unproject.h>
#include <igl\png\writePNG.h>
#include <igl\unproject_ray.h>
#include <igl\file_dialog_save.h>
#include <igl\viewer\ViewerData.h>
#include <igl\unproject_onto_mesh.h>

SolverPlugin::SolverPlugin()
	: solver_wrapper(make_unique<SolverWrapper>()),
	  bar(nullptr)
{
	C << 0.0, 1.0, 0.0;
	C_hover << 0.854902, 0.647059, 0.12549;
	white << 0.7, 0.7, .4;
	red << 1.0, 0.0, 0.0;
	C_merge << 51. / 255., 204. / 255., 1.0;
	zero3 << 0., 0., 0.;
	ones3 << 1., 1., 1.;
	black << 0., 0., 0.;
}

void SolverPlugin::init(igl::viewer::Viewer *viewer)
{
	if (bar == nullptr)
	{
		ViewerPlugin::init(viewer);
		bar = viewer->ngui;
	}

	bar_window = viewer->ngui->window();

	bar->addGroup("Texture Colors");
	bar->addVariable<nanogui::Color>("Color 1",
		[&](nanogui::Color value) { this->viewer->get_mesh(0).tex_col1 = value; this->viewer->get_mesh(0).grid_texture(); },
		[&]() { return (nanogui::Color) this->viewer->get_mesh(0).tex_col1; },
		true);
	bar->addVariable<nanogui::Color>("Color 2",
		[&](nanogui::Color value) { this->viewer->get_mesh(0).tex_col2 = value; this->viewer->get_mesh(0).grid_texture(); },
		[&]() { return (nanogui::Color) this->viewer->get_mesh(0).tex_col2; },
		true);
	add_texture_slider(bar->window(), texture_size, "Texture size");

	param_menu = bar->addWindow(Vec2i(10, 10), "Parameterization");

	bar->addGroup("Solver");

	bar->addButton("Run", [&]()
	{
		start_solver_thread();
	});
		
	add_color_clamp_slider("Dist", max_dist_color_value, dist_color_clamp);
	bar->addVariable("Show RGB", colorByRGB);
	bar->addVariable<bool>("Show distortion",
		[&](bool value) { show_distortion_error = value; update_colors = true; },
		[&]() { return show_distortion_error; },
		true)->setFixedWidth(140);

	bar->addButton("Export UV to OBJ", [&]() { export_uv_to_obj(); });
		
	viewer->screen->performLayout();

	viewer->core.is_animating = true;
	viewer->core.animation_max_fps = 30.0;
	viewer->data.object_scale = 10.0;
}

bool SolverPlugin::load(string filename)
{
	if (solver_wrapper->solver->is_running)
		stop_solver_thread();

	bool read_obj = false;
	bool read_off = false;

	string file_format = filename.substr(filename.length() - 3, 3);
	if (file_format.compare("obj") == 0)
	{
		if (!igl::readOBJ(filename, V, F))
		{
			cerr << "Failed to load mesh: " << filename << endl;
			return false;
		}
	}
	else if (file_format.compare("off") == 0)
	{
		if (!igl::readOFF(filename, V, F))
		{
			cerr << "Failed to load mesh: " << filename << endl;
			return false;
		}
	}
	else
	{
		cerr << "Unknown file format " << filename << endl;
		return false;
	}

	initialize();

	mesh_loaded = true;	
	mesh_filename = filename;

	return true;
}

void SolverPlugin::initialize()
{
	if (solver_wrapper->solver->is_running)
		solver_wrapper->solver->stop();

	while (solver_wrapper->solver->is_running);

	if (V.rows() == 0 || F.rows() == 0)
		return;

	solver_wrapper->init(V, F);

	if (uv_id != 0) {
		viewer->get_mesh(uv_id).clear();
		viewer->get_mesh(mesh_id).clear();
	}
	else
	{
		uv_id = viewer->add_mesh("UV");
	}
	viewer->get_mesh(mesh_id).set_mesh(V, F);
	viewer->get_mesh(uv_id).set_mesh(solver_wrapper->solver->uv, solver_wrapper->solver->F);

	// setup mesh soup
	Mat3 face;
	RVec3i verts;

	RGBColors = V;
	RGBColors.rowwise() -= RGBColors.colwise().minCoeff();
	RGBColors *= RGBColors.colwise().maxCoeff().cwiseInverse().asDiagonal();
	//viewer->get_mesh(mesh_id).set_mesh(V, F);

	//viewer->get_mesh(mesh_id).set_colors(MatX3::Ones(F.rows(), 3));
	viewer->get_mesh(uv_id).F;
	viewer->get_mesh(uv_id).set_colors(MatX3::Ones(solver_wrapper->solver->F.rows(), 3));

	// init colors of uv mesh
	uv_triangle_colors = MatX3::Ones(solver_wrapper->solver->F.rows(), 3);

	viewer->get_mesh(mesh_id).set_colors(RVec3(1., 1., 1.));
}

void SolverPlugin::export_uv_to_obj()
{

}

void SolverPlugin::add_texture_slider(nanogui::Window* window, double& var, const string& name)
{
	Widget *panel = new Widget(window);
	panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));

	Slider* s = new Slider(panel);
	s->setValue(var);

	TextBox *textBox = new TextBox(panel);
	textBox->setValue(removeTrailingZeros(to_string(round(100 * max_texture_val * var) / 100)));
	textBox->setFixedWidth(50);
	textBox->setEditable(false);
	textBox->setFixedHeight(18);

	s->setCallback([&, textBox](float value)
	{
		textBox->setValue(removeTrailingZeros(to_string(round(100 * max_texture_val * value) / 100)));
		var = max_texture_val * value;
		update_colors = true;
	});

	var *= max_texture_val;

	viewer->ngui->addWidget(name, panel);
}

void SolverPlugin::add_color_clamp_slider(const string& name, const shared_ptr<double>& max_value, const shared_ptr<double>& value)
{
	Widget *panel = new Widget(bar->window());
	panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));

	Slider* s = new Slider(panel);
	s->setValue(*value);

	TextBox *textBox = new TextBox(panel);
	textBox->setValue(removeTrailingZeros(to_string(*value)));
	textBox->setFixedWidth(50);
	textBox->setEditable(false);
	textBox->setFixedHeight(18);

	s->setCallback([&, textBox](float new_value) {
		textBox->setValue(removeTrailingZeros(to_string(round(100 * new_value* *max_value) / 100)));
		*value = new_value * *max_value;
		update_colors = true;
	});

	viewer->ngui->addWidget(name, panel);
}

inline string SolverPlugin::removeTrailingZeros(string& s) {
	return s.erase(s.find_last_not_of('0') + 1, string::npos);
}

void SolverPlugin::start_solver_thread()
{
	cout << "start new solver" << endl;
	solver_thread = thread(&Solver::run, solver_wrapper->solver.get());
	solver_thread.detach();
}

void SolverPlugin::stop_solver_thread()
{
	
}

bool SolverPlugin::key_down(int key, int modifiers)
{
	switch (key)
	{
	case GLFW_KEY_T:
		// toggle texturing
		viewer->get_mesh(mesh_id).show_texture = !viewer->get_mesh(mesh_id).show_texture;
		return true; // dont trigger fill variable
	}
	return false;
}

bool SolverPlugin::key_up(int key, int modifiers)
{
	return false;
}

bool SolverPlugin::pre_draw()
{
	if (solver_wrapper->progressed() || update_colors)
		update_mesh();
	return false;
}

void SolverPlugin::update_mesh()
{
	MatX2 newX;
	solver_wrapper->solver->get_mesh(newX);
	viewer->get_mesh(uv_id).set_mesh(newX, solver_wrapper->solver->F);

	// set UV of 3d mesh with newX vertices
	// prepare first for 3d mesh soup
	viewer->get_mesh(mesh_id).set_uv(texture_size * newX);


	if (colorByRGB)
		uv_triangle_colors = RGBColors;
	else
		uv_triangle_colors = MatX3::Ones(3 * F.rows(), 3);

// 	viewer->get_mesh(uv_id).set_normals(viewer->get_mesh(mesh_id).F_normals);
// 	viewer->get_mesh(uv_id).dirty |= viewer->get_mesh(uv_id).DIRTY_NORMAL;

	MatX3 uv_sep_colors = MatX3::Ones(uv_triangle_colors.rows(), 3);
	Vec vals = Vec::Zero(3 * F.rows());


	MatX3 uv_dist_colors = MatX3::Ones(uv_triangle_colors.rows(), 3);
	if (show_distortion_error)
	{
		RVec dist_vals = solver_wrapper->solver->energy->symDirichlet->Efi;

		// new dist color impl
		RVec dist_err = dist_vals.transpose().array() - 4.;

		// scale to [0, dist_cutoff]
		dist_err = dist_err / *dist_color_clamp;
		// map > 1 -> 1
		dist_err= 1 - dist_err.unaryExpr([&](double val) { return (val > 1) ? 1 : val; }).array();

		// map from face to vertex coloring
		Mat3X mat_vertex_dist_vals;
		igl::repmat(dist_err, 3, 1, mat_vertex_dist_vals);
		Vec vertex_dist_vals = Eigen::Map<Vec>(mat_vertex_dist_vals.data(), 3 * mat_vertex_dist_vals.cols(), 1);
		
		uv_dist_colors.col(0) = uv_dist_colors.col(0).cwiseProduct(vertex_dist_vals);
		uv_dist_colors.col(1) = uv_dist_colors.col(1).cwiseProduct(vertex_dist_vals);
	}

#pragma omp parallel for num_threads(3)
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < uv_triangle_colors.rows(); ++j)
		{
			double val = (uv_dist_colors(j, i)+uv_sep_colors(j, i))/2;
			uv_triangle_colors(j, i) *= val;
		}
	}
	

	Vec3i face;
	// draw this first, so a later activation of the hover triangle
	// also overwrites the coloring
	if (hovered_triangle != -1)
	{
		// uv
		face = solver_wrapper->solver->F.row(hovered_triangle);
		uv_triangle_colors.row(face(0)) << C_hover;
		uv_triangle_colors.row(face(1)) << C_hover;
		uv_triangle_colors.row(face(2)) << C_hover;
		// 3d mesh
		mesh_triangle_colors.row(hovered_triangle) << C_hover;
	}


	viewer->get_mesh(uv_id).points.resize(0, Eigen::NoChange);

	viewer->get_mesh(uv_id).dirty |= viewer->get_mesh(uv_id).DIRTY_AMBIENT;
	viewer->get_mesh(mesh_id).dirty |= viewer->get_mesh(mesh_id).DIRTY_AMBIENT;

	update_colors = false;
}

bool SolverPlugin::mouse_down(int button, int modifier)
{
	viewer->core.viewport << 0, 0, 2400, 1350;
	igl::unproject(RVec3(viewer->down_mouse_x, viewer->screen->size()[1] - viewer->down_mouse_y, 0.), (viewer->core.view * viewer->get_mesh(uv_id).model).eval(), viewer->core.proj, viewer->core.viewport, projected_mouse_down);
	mesh_pos_down = V;
	igl::viewer::ViewerData dt = viewer->get_mesh(uv_id);
	uv_mesh_pos_down = viewer->get_mesh(uv_id).V;
	mesh_3d_normals_down = viewer->get_mesh(0).V_normals;

	bool LEFT_ALT = modifier == 4;
	bool LEFT_CTRL = modifier == 2;
	MOUSE_LEFT = button == 0;
	MOUSE_MID = button == 1;
	MOUSE_RIGHT = button == 2;

	if (mouse_on_uv_side)
	{
		if (MOUSE_MID)
		{
			// if the solver is running, stop it during translation, and release it when releaseing the mouse
			if (solver_wrapper->solver->is_running)
			{ // stop and remember to release when done
				solver_wrapper->get_slot();
				release_when_translation_done = true;
			}
			uv_translation_enabled = true;
		}
	}
	else
	{ // mouse on 3d mesh side
		if (MOUSE_LEFT)
		{
			rotation_enabled = true;
		}
		if (MOUSE_MID)
		{
			translation_enabled = true;
		}
	}
	return true;
}

bool SolverPlugin::mouse_up(int button, int modifier)
{
	rotation_enabled = false;
	translation_enabled = false;

	if (uv_translation_enabled)
	{
		MatX2 Vs_new = viewer->get_mesh(uv_id).V.block(0, 0, solver_wrapper->solver->uv.rows(), 2);
		solver_wrapper->solver->uv = Vs_new;
		solver_wrapper->solver->m_x = Eigen::Map<const Vec>(Vs_new.data(), Vs_new.rows() * Vs_new.cols());
		solver_wrapper->solver->ext_x = Eigen::Map<const Vec>(Vs_new.data(), Vs_new.rows() * Vs_new.cols());
		if (release_when_translation_done)
			solver_wrapper->release_slot();
	}
	uv_translation_enabled = false;

	return true;
}

bool SolverPlugin::mouse_scroll(float delta_y)
{
	if (mouse_on_uv_side)
	{ // on uv side
		viewer->core.uv_camera_zoom += delta_y > 0 ? .1 : -.1;
	}
	else
	{ // on mesh side
		viewer->core.mesh_camera_zoom += delta_y > 0 ? .1 : -.1;
	}
	return true;
}

bool SolverPlugin::mouse_move(int mouse_x, int mouse_y)
{
	if (mouse_x <= 1200)
		mouse_on_uv_side = true;
	else
		mouse_on_uv_side = false;
	curr_mouse_x = mouse_x;
	curr_mouse_y = mouse_y;
	return process_mouse_move();
}

bool SolverPlugin::process_mouse_move()
{
	int mouse_x = curr_mouse_x;
	int mouse_y = curr_mouse_y;

	RVec3 curr_mouse_pos, curr_screen_space_mouse_pos;
	curr_screen_space_mouse_pos << mouse_x, viewer->screen->size()[1] - mouse_y, 0.;
	viewer->core.viewport << 0, 0, 2400, 1350;
	igl::unproject(curr_screen_space_mouse_pos, (viewer->core.view * viewer->get_mesh(uv_id).model).eval(), viewer->core.proj, viewer->core.viewport, curr_mouse_pos);
	double x_diff = 2 * (curr_mouse_pos(0) - projected_mouse_down(0));
	double y_diff = curr_mouse_pos(1) - projected_mouse_down(1);
	if (rotation_enabled)
	{
		rotate(-y_diff, -x_diff);
		return true;
	}
	if (translation_enabled)
	{
		translate(x_diff, y_diff);
		return true;
	}
	if (uv_translation_enabled)
	{
		translate_uv_mesh(x_diff, y_diff);
		return true;
	}
	return false;
}

void SolverPlugin::translate(double offset_x, double offset_y)
{
 	MatX3 new_mesh_pos = mesh_pos_down;
	new_mesh_pos.col(0).array() += offset_x;
	new_mesh_pos.col(1).array() += offset_y;
	viewer->get_mesh(mesh_id).set_mesh(new_mesh_pos, F);
	V = new_mesh_pos;
}

void SolverPlugin::translate_uv_mesh(double offset_x, double offset_y)
{
	MatX3 new_mesh_pos = uv_mesh_pos_down;
	new_mesh_pos.col(0).array() += offset_x;
	new_mesh_pos.col(1).array() += offset_y;
	viewer->get_mesh(uv_id).set_vertices(new_mesh_pos);
}

void SolverPlugin::rotate(double phi_x, double phi_y)
{
	Rx.block<2, 2>(1, 1) << cos(phi_x), sin(phi_x), -sin(phi_x), cos(phi_x);
	Rx(0, 0) = 1;
	Ry.row(0) << cos(phi_y), 0, sin(phi_y);
	Ry.row(2) << -sin(phi_y), 0, cos(phi_y);
	Ry(1, 1) = 1;
	Mat3 R = Ry*Rx;
	MatX3 Vtemp = mesh_pos_down;
	RVec3 col_mean = Vtemp.colwise().mean();
	Vtemp.rowwise() -= col_mean;
	MatX3 Vrot = Vtemp * R;
	Vrot.rowwise() += col_mean;
	V = Vrot;
	viewer->get_mesh(mesh_id).set_mesh(Vrot, F);
	viewer->get_mesh(mesh_id).set_normals(mesh_3d_normals_down * R);
}