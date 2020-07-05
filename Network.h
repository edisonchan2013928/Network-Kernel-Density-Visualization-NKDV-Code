#define _CRT_SECURE_NO_WARNINGS
#pragma once
#ifndef NETWORK_H
#define NETWORK_H

#include "Library.h"

//#define PROFILING
#define STATISTICS

struct Point
{
	double dist_n1;
	double dist_n2;
	int edge_index;

	//Used in ADA method
	double aug_sum_dist_c_p_deg_1;
	double aug_sum_dist_c_p_deg_2;
	double aug_sum_dist_d_p_deg_1;
	double aug_sum_dist_d_p_deg_2;
};

struct Edge
{
	int n1;
	int n2;
	double length;
	vector<Point> PS;

	//Used in ADA method
	vector<double> dist_n1_vec; //which_vector = 0
	vector<double> dist_n2_vec; //which_vector = 1
	vector<double> aug_dist_diff_vec; //which_vector = 2

	//Used in IA method
	bool is_IA;
	int num_intervals;
	double min_dist_diff;
	double last_interval_size;
	vector<double> aug_sum_dist_c_I_deg_1_vec;
	vector<double> aug_sum_dist_c_I_deg_2_vec;
	vector<double> aug_sum_dist_d_I_deg_1_vec;
	vector<double> aug_sum_dist_d_I_deg_2_vec;
	vector<double> interval_point_vec;
	vector<double> weight_vec;
	vector<double> agg_weight_c_vec;
	vector<double> agg_weight_d_vec;
};

struct Lixel : Point
{
	double KDE_value;
};

struct sp_node
{
	int node_index;
	double cur_sp_value;
	bool is_opt;
};

//This is the minimum heap for dijkstra's algorithm
struct comparePriority
{
	bool operator()(sp_node& n1, sp_node& n2)
	{
		return n1.cur_sp_value > n2.cur_sp_value;
	}
};

struct model
{
	int method;
	int n, m; //number of nodes and edges
	double lixel_reg_length;
	char* network_fileName;
	char* out_NKDV_fileName;
	Edge* edge_set;
	vector<vector<int> > Network;
	vector<Lixel> lixel_set;

	//kernel functions
	int k_type; //k_type = 0: Gaussian kernel, k_type = 1: Triangular kernel 2: Epanechnikov kernel 3: Quartic kernel
	double bandwidth;
	double gamma; //bandwidth can determine the parameter gamma.

	//Used in Dijkstra's algorithm
	Lixel cur_l;
	vector<sp_node> sp_node_vec;

	//Used in method = 2
	int sel_node_index;
	int node_index_pi_a;
	int node_index_pi_b;
	vector<sp_node> sp_node_vec_node_a;
	vector<sp_node> sp_node_vec_node_b;

	//Used in IA
	double smallest_interval_size;

	#ifdef STATISTICS
		double num_of_augmentation;
		double num_of_edges_w_points;
	#endif
};

#endif