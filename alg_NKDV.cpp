#include "alg_NKDV.h"

void NKDV_algorithm(model& our_model)
{
	#ifdef PROFILING
	double SP_time = 0;
	double edge_KAQ_time = 0;
	#endif

	//Used in method = 2
	int prev_edge_index = -1;
	int num_edge_index = 0;
	int edge_index;

	double run_time;
	//Different algorithms
	auto start_s = chrono::high_resolution_clock::now();
	if (our_model.method == 1) //better version of NKDV baseline method
	{
		for (int l = 0; l < (int)our_model.lixel_set.size(); l++)
		{
			our_model.cur_l = our_model.lixel_set[l];
			#ifdef PROFILING
				auto start_SP_time = chrono::high_resolution_clock::now();
			#endif

			dijkstra(our_model);

			#ifdef PROFILING
				auto end_SP_time = chrono::high_resolution_clock::now();
			#endif

			NKDV_basic(our_model);

			#ifdef PROFILING
				auto end_edge_KAQ_time = chrono::high_resolution_clock::now();
				SP_time += (end_SP_time - start_SP_time).count();
				edge_KAQ_time += (end_edge_KAQ_time - end_SP_time).count();
			#endif

			our_model.lixel_set[l] = our_model.cur_l;
		}
	}
	if (our_model.method >= 2 && our_model.method <= 5) //Save the number of calls in dijstra algorithm
	{
		if (our_model.method == 3 || our_model.method == 5)
			augment_preprocess(our_model);
		if (our_model.method == 4 || our_model.method == 5)
			augment_interval_preprocess(our_model);

		for (int l = 0; l < (int)our_model.lixel_set.size(); l++)
		{
			our_model.cur_l = our_model.lixel_set[l];
			edge_index = our_model.lixel_set[l].edge_index;

			#ifdef PROFILING
				auto start_SP_time = chrono::high_resolution_clock::now();
			#endif

			if (edge_index != prev_edge_index)
			{
				prev_edge_index = edge_index;
				our_model.sel_node_index = our_model.edge_set[edge_index].n1;
				our_model.node_index_pi_a = our_model.sel_node_index;
				dijkstra(our_model);
				copy_sp_info(our_model, true);

				our_model.sel_node_index = our_model.edge_set[edge_index].n2;
				our_model.node_index_pi_b = our_model.sel_node_index;
				dijkstra(our_model);
				copy_sp_info(our_model, false);
			}

			#ifdef PROFILING
				auto end_SP_time = chrono::high_resolution_clock::now();
			#endif

			NKDV_basic(our_model);

			#ifdef PROFILING
				auto end_edge_KAQ_time = chrono::high_resolution_clock::now();
				SP_time += (end_SP_time - start_SP_time).count();
				edge_KAQ_time += (end_edge_KAQ_time - end_SP_time).count();
			#endif

			our_model.lixel_set[l] = our_model.cur_l;
		}
	}

	auto end_s = chrono::high_resolution_clock::now();

	run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;

	std::cout << "method " << our_model.method << ":" << run_time << endl;
	#ifdef PROFILING
	cout << "SP time: " << SP_time / 1000000000.0 << endl;
	cout << "edge KAQ time: " << edge_KAQ_time / 1000000000.0 << endl;
	#endif
}