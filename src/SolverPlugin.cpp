#include "SolverPlugin.h"

#include "svg_exporter.h"

#include <queue>
#include <deque>
#include <array>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <igl/writeOBJ.h>
#include <igl/unproject.h>
#include <igl/unproject_ray.h>
#include <igl/file_dialog_save.h>
#include <igl/opengl/ViewerData.h>
#include <igl/unproject_onto_mesh.h>

SolverPlugin::SolverPlugin()
	: solver_wrapper(make_unique<SolverWrapper>())
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


void SolverPlugin::init(Viewer *viewer)
{
	leftView = 0;
	viewer->core->background_color << 1., 1., 1., 1.; // 0.25, 0.25, 0.25, 1.0;
	viewer->core->orthographic = true;
	viewer->core->viewport = Eigen::Vector4f(0, 0, 960, 1080);
	rightView = viewer->append_core(Eigen::Vector4f(960, 0, 960, 1080));

	this->viewer = viewer;
	viewer->plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&,viewer]()
	{
		menu.draw_viewer_menu();
		if (ImGui::CollapsingHeader("Texture Colors", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::ColorEdit4("Color 1", viewer->data_list[0].tex_col1.data(),
				ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel) ||
				ImGui::ColorEdit4("Color 2", viewer->data_list[0].tex_col2.data(),
					ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel))
				viewer->data_list[0].grid_texture();
			ImGui::DragFloat("Texture Size", &texture_size, 0.25f);
		}
	};
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"Solver", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);
		if (ImGui::Button("Start")) 
			start_solver_thread();
		ImGui::SameLine(); 
		if (ImGui::Button("Stop"))
			stop_solver_thread();


		ImGui::End();
	};

	viewer->core->is_animating = true;
	viewer->core->animation_max_fps = 30.0;
// 	viewer->data.object_scale = 10.0;
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
	viewer->coreDataPairs.clear();
	viewer->coreDataPairs.insert({ 0,1 });
	viewer->coreDataPairs.insert({ 1,0 });
	viewer->selected_data_index = 0;
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
		viewer->data_list[uv_id].clear();
		viewer->data_list[mesh_id].clear();
	}
	else
	{
		uv_id = viewer->append_mesh();
	}
	viewer->data_list[mesh_id].set_mesh(V, F);
	viewer->data_list[mesh_id].set_uv(solver_wrapper->solver->uv);
	viewer->data_list[uv_id].set_mesh(solver_wrapper->solver->uv, solver_wrapper->solver->F);

	RGBColors = V;
	RGBColors.rowwise() -= RGBColors.colwise().minCoeff();
	RGBColors *= RGBColors.colwise().maxCoeff().cwiseInverse().asDiagonal();

	viewer->data_list[uv_id].F;
	viewer->data_list[uv_id].set_colors(MatX3::Ones(solver_wrapper->solver->F.rows(), 3));

	// init colors of uv mesh
	uv_triangle_colors = MatX3::Ones(solver_wrapper->solver->F.rows(), 3);

	viewer->data_list[mesh_id].set_colors(RVec3(1., 1., 1.));
}

void SolverPlugin::export_uv_to_obj()
{

}

// void SolverPlugin::add_texture_slider(nanogui::Window* window, double& var, const string& name)
// {
// 	Widget *panel = new Widget(window);
// 	panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));
// 
// 	Slider* s = new Slider(panel);
// 	s->setValue(var);
// 
// 	TextBox *textBox = new TextBox(panel);
// 	textBox->setValue(removeTrailingZeros(to_string(round(100 * max_texture_val * var) / 100)));
// 	textBox->setFixedWidth(50);
// 	textBox->setEditable(false);
// 	textBox->setFixedHeight(18);
// 
// 	s->setCallback([&, textBox](float value)
// 	{
// 		textBox->setValue(removeTrailingZeros(to_string(round(100 * max_texture_val * value) / 100)));
// 		var = max_texture_val * value;
// 		update_colors = true;
// 	});
// 
// 	var *= max_texture_val;
// 
// 	viewer->ngui->addWidget(name, panel);
// }
// 
// void SolverPlugin::add_color_clamp_slider(const string& name, const shared_ptr<double>& max_value, const shared_ptr<double>& value)
// {
// 	Widget *panel = new Widget(bar->window());
// 	panel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));
// 
// 	Slider* s = new Slider(panel);
// 	s->setValue(*value);
// 
// 	TextBox *textBox = new TextBox(panel);
// 	textBox->setValue(removeTrailingZeros(to_string(*value)));
// 	textBox->setFixedWidth(50);
// 	textBox->setEditable(false);
// 	textBox->setFixedHeight(18);
// 
// 	s->setCallback([&, textBox](float new_value) {
// 		textBox->setValue(removeTrailingZeros(to_string(round(100 * new_value* *max_value) / 100)));
// 		*value = new_value * *max_value;
// 		update_colors = true;
// 	});
// 
// 	viewer->ngui->addWidget(name, panel);
// }
// 
// inline string SolverPlugin::removeTrailingZeros(string s) {
// 	return s.erase(s.find_last_not_of('0') + 1, string::npos);
// }

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
		viewer->data_list[mesh_id].show_texture = !viewer->data_list[mesh_id].show_texture;
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
	viewer->data_list[uv_id].set_mesh(newX, solver_wrapper->solver->F);

	// set UV of 3d mesh with newX vertices
	// prepare first for 3d mesh soup
	viewer->data_list[mesh_id].set_uv(texture_size * newX);


	if (colorByRGB)
		uv_triangle_colors = RGBColors;
	else
		uv_triangle_colors = MatX3::Ones(F.rows(), 3);

// 	viewer->get_mesh(uv_id).set_normals(viewer->get_mesh(mesh_id).F_normals);
// 	viewer->get_mesh(uv_id).dirty |= viewer->get_mesh(uv_id).DIRTY_NORMAL;

	MatX3 uv_sep_colors = MatX3::Ones(uv_triangle_colors.rows(), 3);
	Vec vals = Vec::Zero(3 * F.rows());


	MatX3 uv_dist_colors = MatX3::Ones(uv_triangle_colors.rows(), 3);
	if (show_distortion_error)
	{
		Vec dist_vals = solver_wrapper->solver->energy->symDirichlet->Efi;

		// new dist color impl
		Vec dist_err = dist_vals.transpose().array() - 4.;

		// scale to [0, dist_cutoff]
		dist_err = dist_err / *dist_color_clamp;
		// map > 1 -> 1
		dist_err= 1 - dist_err.unaryExpr([&](double val) { return (val > 1) ? 1 : val; }).array();
		
		uv_dist_colors.col(1) = uv_dist_colors.col(1).cwiseProduct(dist_err);
		uv_dist_colors.col(2) = uv_dist_colors.col(2).cwiseProduct(dist_err);
	}


	uv_triangle_colors.array() *= uv_dist_colors.array();
	

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


	viewer->data_list[uv_id].points.resize(0, Eigen::NoChange);

	viewer->data_list[uv_id].set_colors(uv_triangle_colors);
	viewer->data_list[mesh_id].set_colors(uv_triangle_colors);

	viewer->data_list[uv_id].dirty |= igl::opengl::MeshGL::DIRTY_AMBIENT;
	viewer->data_list[mesh_id].dirty |= igl::opengl::MeshGL::DIRTY_AMBIENT;

	update_colors = false;
}

bool SolverPlugin::mouse_down(int button, int modifier)
{
	return false;
}

bool SolverPlugin::mouse_up(int button, int modifier)
{
	rotation_enabled = false;
	translation_enabled = false;

	if (uv_translation_enabled)
	{
		MatX2 Vs_new = viewer->data_list[uv_id].V.block(0, 0, solver_wrapper->solver->uv.rows(), 2);
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
	return false;
}

bool SolverPlugin::mouse_move(int mouse_x, int mouse_y)
{

	return false;
}

bool SolverPlugin::process_mouse_move()
{
	return false;
}

void SolverPlugin::translate(double offset_x, double offset_y)
{
 	MatX3 new_mesh_pos = mesh_pos_down;
	new_mesh_pos.col(0).array() += offset_x;
	new_mesh_pos.col(1).array() += offset_y;
	viewer->data_list[mesh_id].set_mesh(new_mesh_pos, F);
	V = new_mesh_pos;
}

void SolverPlugin::translate_uv_mesh(double offset_x, double offset_y)
{
	MatX3 new_mesh_pos = uv_mesh_pos_down;
	new_mesh_pos.col(0).array() += offset_x;
	new_mesh_pos.col(1).array() += offset_y;
	viewer->data_list[uv_id].set_vertices(new_mesh_pos);
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
	viewer->data_list[mesh_id].set_mesh(Vrot, F);
	viewer->data_list[mesh_id].set_normals(mesh_3d_normals_down * R);
}