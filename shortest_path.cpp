#include "shortest_path.h"

inline void clearHeap(PQ& pq)
{
	int heapSize = (int)pq.size();
	for (int h = 0; h < heapSize; h++)
		pq.pop();
}

void init_dijkstra(model& our_model, PQ& pq)
{
	int n1, n2;
	Lixel& cur_l = our_model.cur_l;

	for (int v = 0; v < our_model.n; v++)
	{
		our_model.sp_node_vec[v].node_index = v;
		our_model.sp_node_vec[v].cur_sp_value = inf;
		our_model.sp_node_vec[v].is_opt = false;
	}

	if (our_model.method == 1) //basic method
	{
		n1 = our_model.edge_set[cur_l.edge_index].n1;
		n2 = our_model.edge_set[cur_l.edge_index].n2;
		our_model.sp_node_vec[n1].cur_sp_value = cur_l.dist_n1;
		our_model.sp_node_vec[n2].cur_sp_value = cur_l.dist_n2;

		pq.push(our_model.sp_node_vec[n1]);
		pq.push(our_model.sp_node_vec[n2]);
	}
	if (our_model.method >= 2 && our_model.method <= 5) //share SPD value
	{
		our_model.sp_node_vec[our_model.sel_node_index].cur_sp_value = 0;
		pq.push(our_model.sp_node_vec[our_model.sel_node_index]);
	}
}

void dijkstra(model& our_model)
{
	static PQ pq;
	int node_index, node_index_neighbor;
	int edge_index;
	double length;
	sp_node SP_node;

	init_dijkstra(our_model, pq);

	while (pq.size() > 0)
	{
		SP_node = pq.top();
		pq.pop();

		//discard outdated node
		if (our_model.sp_node_vec[SP_node.node_index].cur_sp_value < SP_node.cur_sp_value - eps)
			continue;

		//discard node with value more than the bandwidth
		if (our_model.sp_node_vec[SP_node.node_index].cur_sp_value > our_model.bandwidth)
			continue;

		our_model.sp_node_vec[SP_node.node_index].is_opt = true;
		node_index = SP_node.node_index;

		for (int e = 0; e < (int)our_model.Network[node_index].size(); e++)
		{
			//Relax (code here)
			edge_index = our_model.Network[node_index][e];
			if (node_index != our_model.edge_set[edge_index].n1)
				node_index_neighbor = our_model.edge_set[edge_index].n1;
			else
				node_index_neighbor = our_model.edge_set[edge_index].n2;

			if (our_model.sp_node_vec[node_index_neighbor].is_opt == true)
				continue;

			length = our_model.edge_set[edge_index].length;
			if (our_model.sp_node_vec[node_index_neighbor].cur_sp_value > our_model.sp_node_vec[node_index].cur_sp_value + length)
			{
				our_model.sp_node_vec[node_index_neighbor].cur_sp_value = our_model.sp_node_vec[node_index].cur_sp_value + length;
				pq.push(our_model.sp_node_vec[node_index_neighbor]);
			}
		}
	}
}

void copy_sp_info(model& our_model, bool isFirst)
{
	for (int v = 0; v < our_model.n; v++)
	{
		if (isFirst == true)
		{
			our_model.sp_node_vec_node_a[v].node_index = our_model.sp_node_vec[v].node_index;
			our_model.sp_node_vec_node_a[v].is_opt = our_model.sp_node_vec[v].is_opt;
			our_model.sp_node_vec_node_a[v].cur_sp_value = our_model.sp_node_vec[v].cur_sp_value;
		}
		else
		{
			our_model.sp_node_vec_node_b[v].node_index = our_model.sp_node_vec[v].node_index;
			our_model.sp_node_vec_node_b[v].is_opt = our_model.sp_node_vec[v].is_opt;
			our_model.sp_node_vec_node_b[v].cur_sp_value = our_model.sp_node_vec[v].cur_sp_value;
		}
	}
}