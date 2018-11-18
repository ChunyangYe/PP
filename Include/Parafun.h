#pragma once

#include "Eigen/Core"
#include<math.h>
#include<vector>
#include<iostream>
#include<fstream>
#include"Eigen/Dense"
#include"Eigen/Sparse"
#include"Eigen/SparseLU"
#include "MeshDefinition.h"
#include "OpenMesh/Core/IO/MeshIO.hh"

#include "PardisoSolver.h"
using namespace Eigen;
using namespace std;
double get_smallest_pos_quad_zero(double a, double b, double c);
class Parafun
{
public:
	Parafun();
	Parafun(string filename);
	~Parafun();
	bool writeobj(Mesh& out_mesh,const string out_filename);

	void BPE();
	void calc_gradient_norm(const VectorXd &x);

	void init();
	void recover_to_src();
	void init_area();

	void run_bpe();

	void Update_source_same_t();

	void Tutte();
	void Pre_calculate();

	void CM();
	void SLIM();

	void Energy(const VectorXd &x, double &energy);
	void Energysource();

	double newton_equation(const double & a, const double & b, const double & K);

	void backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha);

	void local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11);

	void max_step(const VectorXd &xx, const VectorXd &dd, double &step);


	Mesh mesh;
	double Intp_T_Min;
	double changetocm_flag;
	double convgence_con_rate;
	double time_consumption;
	int MAX_ITER_NUM;
	string file_str;
	string modelname;
	double originmesh_area_sqrt;

	int F_N;
	int V_N;

	VectorXd negative_grad_norm;
	double g_norm;

	vector<double> area;
	vector<double> area_uniform;
	vector<double> area_src;

	vector<double> source_p00;
	vector<double> source_p01;
	vector<double> source_p10;
	vector<double> source_p11;

	vector<double> update_p00;
	vector<double> update_p01;
	vector<double> update_p10;
	vector<double> update_p11;

	vector<int> F0;
	vector<int> F1;
	vector<int> F2;

	PardisoSolver* pardiso;
	vector<int> pardiso_ia;
	vector<int> pardiso_i;
	vector<int> pardiso_ja;
	vector<double> pardiso_a;
	vector<double> pardiso_b;

	double energy_uniform;
	double energy_area;

	double bound_distortion_K;
	VectorXd position_of_mesh;

	vector<int> id_h00; vector<int> id_h01; vector<int> id_h02; vector<int> id_h03; vector<int> id_h04; vector<int> id_h05;
	vector<int> id_h11; vector<int> id_h12; vector<int> id_h13; vector<int> id_h14; vector<int> id_h15;
	vector<int> id_h22; vector<int> id_h23; vector<int> id_h24; vector<int> id_h25;
	vector<int> id_h33; vector<int> id_h34; vector<int> id_h35;
	vector<int> id_h44; vector<int> id_h45;
	vector<int> id_h55;
};

