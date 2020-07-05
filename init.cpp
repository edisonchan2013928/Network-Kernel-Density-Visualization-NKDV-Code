#include "init.h"

void load_network(model& our_model)
{
	fstream network_file;
	int num_points;
	int degree;
	int edge_index;
	Point p;
	vector<int> tempVector;
	sp_node temp_sp_Node;

	network_file.open(our_model.network_fileName);
	if (network_file.is_open() == false)
	{
		cout << "Cannot open network file!" << endl;
		exit(0);
	}

	network_file >> our_model.m;
	our_model.edge_set = new Edge[our_model.m];

	for (int e = 0; e < our_model.m; e++)
	{
		network_file >> our_model.edge_set[e].n1;
		network_file >> our_model.edge_set[e].n2;
		network_file >> our_model.edge_set[e].length;
		network_file >> num_points;

		for (int i = 0; i < num_points; i++)
		{
			p.edge_index = e;
			network_file >> p.dist_n1;
			network_file >> p.dist_n2;
			our_model.edge_set[e].PS.push_back(p);
			our_model.edge_set[e].aug_dist_diff_vec.push_back(-inf);
			our_model.edge_set[e].dist_n1_vec.push_back(-1);
			our_model.edge_set[e].dist_n2_vec.push_back(-1);
		}
	}

	network_file >> our_model.n;
	for (int v = 0; v < our_model.n; v++)
	{
		our_model.Network.push_back(tempVector);
		network_file >> degree;
		for (int e = 0; e < degree; e++)
		{
			network_file >> edge_index;
			our_model.Network[v].push_back(edge_index);
		}
	}

	//init the Dijkstra's algorithm
	for (int v = 0; v < our_model.n; v++)
	{
		our_model.sp_node_vec.push_back(temp_sp_Node);

		if (our_model.method >= 2 && our_model.method <= 5)
		{
			our_model.sp_node_vec_node_a.push_back(temp_sp_Node);
			our_model.sp_node_vec_node_b.push_back(temp_sp_Node);
		}
	}

	network_file.close();
}

void add_lixels_for_edge(int e, model& our_model)
{
	double cur_dist = 0;
	double middle_dist;
	double next_dist;
	double length = our_model.edge_set[e].length;
	Lixel l;
	
	while (cur_dist < length)
	{
		next_dist = cur_dist + our_model.lixel_reg_length;
		if (next_dist > length)
			next_dist = length;

		middle_dist = (cur_dist + next_dist) / 2.0;
		l.dist_n1 = middle_dist;
		l.dist_n2 = length - middle_dist;
		l.edge_index = e;
		l.KDE_value = -100;
		our_model.lixel_set.push_back(l);

		cur_dist += our_model.lixel_reg_length;
	}
}

void obtain_lixel_set(model& our_model)
{
	for (int e = 0; e < our_model.m; e++)
		add_lixels_for_edge(e, our_model);

	//check the number of lixels
	//cout << our_model.lixel_set.size() << endl;
}

void init_parameters(int argc, char** argv, model& our_model)
{
	//debug
	/*our_model.network_fileName = (char*)"../../../Datasets/Testing/Testing_network";
	our_model.out_NKDV_fileName = (char*)"Results/Testing_M3_K1";
	our_model.method = 3;
	our_model.lixel_reg_length = 1;
	our_model.k_type = 1;
	our_model.bandwidth = 8;*/
	our_model.network_fileName = argv[1];
	our_model.out_NKDV_fileName = argv[2];
	our_model.method = atoi(argv[3]);
	our_model.lixel_reg_length = atoi(argv[4]);
	our_model.k_type = atoi(argv[5]);
	our_model.bandwidth = atof(argv[6]);
	our_model.gamma = 1;

	//Gaussian kernel
	if (our_model.k_type == 0)
		our_model.gamma = 9 / (2 * our_model.bandwidth * our_model.bandwidth);
	//Triangular kernel
	if (our_model.k_type == 1)
		our_model.gamma = 1 / our_model.bandwidth;
	//Epanechnikov and Quartic kernels
	if (our_model.k_type == 2 || our_model.k_type == 3)
		our_model.gamma = 1 / (our_model.bandwidth*our_model.bandwidth);
}

void output_Visual(model& our_model)
{
	int edge_index;
	double dist_n1, dist_n2;
	double KDE_value;
	fstream out_NKDV_file;

	out_NKDV_file.open(our_model.out_NKDV_fileName, ios::in | ios::out | ios::trunc);
	if (out_NKDV_file.is_open() == false)
	{
		cout << "Cannot open output file!" << endl;
		exit(0);
	}

	out_NKDV_file << our_model.lixel_set.size() << endl;
	for (int l = 0; l < (int)our_model.lixel_set.size(); l++)
	{
		edge_index = our_model.lixel_set[l].edge_index;
		dist_n1 = our_model.lixel_set[l].dist_n1;
		dist_n2 = our_model.lixel_set[l].dist_n2;
		KDE_value = our_model.lixel_set[l].KDE_value;
		out_NKDV_file << edge_index << " " << dist_n1 << " " << dist_n2 << " " << KDE_value << endl;
	}

	out_NKDV_file.close();
}