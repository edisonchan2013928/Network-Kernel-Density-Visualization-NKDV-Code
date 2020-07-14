#include "KAF.h"

double kernel_value(model& our_model, double dist)
{
	double sq_dist = dist * dist;

	if (our_model.k_type == 0) //Gaussian kernel
		return exp(-our_model.gamma*sq_dist);

	//Triangular, Epanechnikov and Quartic kernels
	if (our_model.k_type >= 1 && our_model.k_type <= 3)
	{
		if (dist >= our_model.bandwidth)
			return 0;

		if (our_model.k_type == 1) //Triangular kernel
			return (1 - our_model.gamma*dist);
		if (our_model.k_type == 2) //Epanechnikov kernel
			return (1 - our_model.gamma*sq_dist);
		if (our_model.k_type == 3) //Quartic kernel
			return (1 - our_model.gamma*sq_dist)*(1 - our_model.gamma*sq_dist);
	}

	return 0;
}

//find the correct position in the array (Used in ADA)
int find_position(double value, Edge& e, int which_vector)
{
	vector<Point>& PS = e.PS;
	int position;
	vector<double>::iterator position_iterator;

	if (which_vector == 0)
	{
		position_iterator = lower_bound(e.dist_n1_vec.begin(), e.dist_n1_vec.end(), value);
		position = (int)(position_iterator - e.dist_n1_vec.begin()) - 1;
	}
	if (which_vector == 1)
	{
		position_iterator = lower_bound(e.dist_n2_vec.begin(), e.dist_n2_vec.end(), value, greater<double>());
		position = (int)(position_iterator - e.dist_n2_vec.begin());
	}
	if (which_vector == 2)
	{
		position_iterator = lower_bound(e.aug_dist_diff_vec.begin(), e.aug_dist_diff_vec.end(), value);
		position = (int)(position_iterator - e.aug_dist_diff_vec.begin()) - 1;
	}

	return position;
}

//Used in IA
void find_interval_same_edge(model& our_model, int edge_index, int& left_interval, int& q_interval, int& right_interval, double& edge_KA_value)
{
	double left_value, right_value, q_value;
	double dist_q_p;
	int N_e = our_model.edge_set[edge_index].num_intervals;
	Edge& edge = our_model.edge_set[edge_index];
	Lixel& cur_l = our_model.cur_l;
	bool is_added = true;
	int left_interval_ori, q_interval_ori, right_interval_ori;

	left_value = cur_l.dist_n1 - our_model.bandwidth;
	right_value = cur_l.dist_n1 + our_model.bandwidth;
	q_value = cur_l.dist_n1;

	left_interval = (int)max(floor(left_value / edge.min_dist_diff), 0.0);
	q_interval = (int)floor(q_value / edge.min_dist_diff);
	right_interval = (int)min(floor(right_value / edge.min_dist_diff), (double)N_e - 1.0);

	left_interval_ori = left_interval;
	q_interval_ori = q_interval;
	right_interval_ori = right_interval;

	if (edge.weight_vec[left_interval] > 0 && left_value < edge.interval_point_vec[left_interval])
	{
		dist_q_p = cur_l.dist_n1 - edge.interval_point_vec[left_interval];
		edge_KA_value += edge.weight_vec[left_interval] * (1 - our_model.gamma*dist_q_p*dist_q_p);
		left_interval++;
	}

	if (edge.weight_vec[q_interval] > 0 && left_interval_ori < q_interval_ori)
	{
		dist_q_p = edge.interval_point_vec[q_interval] - cur_l.dist_n1;
		if (dist_q_p < our_model.bandwidth)
			edge_KA_value += edge.weight_vec[q_interval] * (1 - our_model.gamma*dist_q_p*dist_q_p);
	}

	if (edge.weight_vec[right_interval] > 0 && right_value > edge.interval_point_vec[right_interval])
	{
		if (q_interval_ori < right_interval_ori)
		{
			dist_q_p = edge.interval_point_vec[right_interval] - cur_l.dist_n1;
			edge_KA_value += edge.weight_vec[right_interval] * (1 - our_model.gamma*dist_q_p*dist_q_p);
		}
		
		right_interval--;
	}
}

double edge_KAF(model& our_model, int edge_index)
{
	int node_index_pi_c = our_model.edge_set[edge_index].n1;
	int node_index_pi_d = our_model.edge_set[edge_index].n2;
	int node_index_pi_a = our_model.node_index_pi_a;
	int node_index_pi_b = our_model.node_index_pi_b;
	double dist_q_p;
	Lixel& cur_l = our_model.cur_l;

	double edge_KA_value = 0;
	Edge& e = our_model.edge_set[edge_index];
	vector<Point>& PS = e.PS;
	if (PS.size() == 0)
		return 0;

	//Interval Augmentation (IA)
	if (our_model.method == 4)
	{
		int left_interval, q_interval, right_interval;
		double a_I_left_deg_1, a_I_left_deg_2;
		double a_I_right_deg_1, a_I_right_deg_2;
		double weight_c, weight_d;
		double dist_q_c, dist_q_d;
		double left_value, right_value;
		int N_e = our_model.edge_set[edge_index].num_intervals;
		int v;
		Edge& edge = our_model.edge_set[edge_index];

		//q is in the same edge
		if (cur_l.edge_index == edge_index)
		{
			find_interval_same_edge(our_model, edge_index, left_interval, q_interval, right_interval, edge_KA_value);

			//No need to consider the case when left_interval >= q_interval
			if (left_interval < q_interval)
			{
				a_I_left_deg_1 = edge.aug_sum_dist_d_I_deg_1_vec[left_interval] - edge.aug_sum_dist_d_I_deg_1_vec[q_interval];
				a_I_left_deg_2 = edge.aug_sum_dist_d_I_deg_2_vec[left_interval] - edge.aug_sum_dist_d_I_deg_2_vec[q_interval];
				weight_d = edge.agg_weight_d_vec[left_interval] - edge.agg_weight_d_vec[q_interval];

				edge_KA_value += (1 - our_model.gamma*cur_l.dist_n2*cur_l.dist_n2)*weight_d + 2 * our_model.gamma*cur_l.dist_n2*a_I_left_deg_1 - our_model.gamma*a_I_left_deg_2;
			}
				
			if (right_interval > q_interval)
			{
				a_I_right_deg_1 = edge.aug_sum_dist_c_I_deg_1_vec[right_interval] - edge.aug_sum_dist_c_I_deg_1_vec[q_interval];
				a_I_right_deg_2 = edge.aug_sum_dist_c_I_deg_2_vec[right_interval] - edge.aug_sum_dist_c_I_deg_2_vec[q_interval];
				weight_c = edge.agg_weight_c_vec[right_interval] - edge.agg_weight_c_vec[q_interval];
				
				edge_KA_value += (1 - our_model.gamma*cur_l.dist_n1*cur_l.dist_n1)*weight_c + 2 * our_model.gamma*cur_l.dist_n1*a_I_right_deg_1 - our_model.gamma*a_I_right_deg_2;
			}

			return edge_KA_value;
		}

		//q is in different edge
		//compute the distances between q to c and q to d
		dist_q_c = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_c].cur_sp_value,
			cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_c].cur_sp_value);
		dist_q_d = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_d].cur_sp_value,
			cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_d].cur_sp_value);

		//Case 1:
		if (dist_q_c > our_model.bandwidth && dist_q_d > our_model.bandwidth)
			return 0;
		//Case 2:
		if (dist_q_c <= our_model.bandwidth && dist_q_d > our_model.bandwidth)
		{
			left_value = our_model.bandwidth - dist_q_c;
			left_interval = (int)floor(left_value / edge.min_dist_diff);
			if (left_interval > 0)
				edge_KA_value += (1 - our_model.gamma*dist_q_c*dist_q_c)*edge.agg_weight_c_vec[left_interval - 1]
					- 2 * our_model.gamma*dist_q_c*edge.aug_sum_dist_c_I_deg_1_vec[left_interval - 1]
					- our_model.gamma*edge.aug_sum_dist_c_I_deg_2_vec[left_interval - 1];
				

			if (edge.weight_vec[left_interval] > 0)
			{
				dist_q_p = dist_q_c + edge.interval_point_vec[left_interval];
				if (dist_q_p <= our_model.bandwidth)
					edge_KA_value += edge.weight_vec[left_interval] * (1 - our_model.gamma*dist_q_p*dist_q_p);
			}

			return edge_KA_value;
		}
		//Case 3:
		if (dist_q_c > our_model.bandwidth && dist_q_d <= our_model.bandwidth)
		{
			right_value = our_model.bandwidth - dist_q_d;
			left_value = edge.length - right_value; //change to left value
			right_interval = (int)floor(left_value / edge.min_dist_diff);
			if (right_interval < edge.num_intervals - 1)
			{
				edge_KA_value += (1 - our_model.gamma*dist_q_d*dist_q_d)*edge.agg_weight_d_vec[right_interval + 1]
					- 2 * our_model.gamma*dist_q_d*edge.aug_sum_dist_d_I_deg_1_vec[right_interval + 1]
					- our_model.gamma*edge.aug_sum_dist_d_I_deg_2_vec[right_interval + 1];
			}
				

			if (edge.weight_vec[right_interval] > 0)
			{
				dist_q_p = dist_q_d + (edge.length - edge.interval_point_vec[right_interval]);
				if (dist_q_p <= our_model.bandwidth)
					edge_KA_value += edge.weight_vec[right_interval] * (1 - our_model.gamma*dist_q_p*dist_q_p);
			}

			return edge_KA_value;
		}

		//Case 4:
		if (dist_q_c <= our_model.bandwidth && dist_q_d <= our_model.bandwidth)
		{
			left_value = our_model.bandwidth - dist_q_c;
			right_value = our_model.bandwidth - dist_q_d;
			left_interval = (int)floor(left_value / edge.min_dist_diff);
			right_interval = (int)floor((edge.length - right_value) / edge.min_dist_diff);

			if (left_interval < right_interval)
			{
				if (left_interval > 0)
					edge_KA_value += (1 - our_model.gamma*dist_q_c*dist_q_c)*edge.agg_weight_c_vec[left_interval - 1]
					- 2 * our_model.gamma*dist_q_c*edge.aug_sum_dist_c_I_deg_1_vec[left_interval - 1]
					- our_model.gamma*edge.aug_sum_dist_c_I_deg_2_vec[left_interval - 1];

				if (edge.weight_vec[left_interval] > 0)
				{
					dist_q_p = dist_q_c + edge.interval_point_vec[left_interval];
					if (dist_q_p <= our_model.bandwidth)
						edge_KA_value += edge.weight_vec[left_interval] * (1 - our_model.gamma*dist_q_p*dist_q_p);
				}

				if (right_interval < edge.num_intervals - 1)
					edge_KA_value += (1 - our_model.gamma*dist_q_d*dist_q_d)*edge.agg_weight_d_vec[right_interval + 1]
					- 2 * our_model.gamma*dist_q_d*edge.aug_sum_dist_d_I_deg_1_vec[right_interval + 1]
					- our_model.gamma*edge.aug_sum_dist_d_I_deg_2_vec[right_interval + 1];

				if (edge.weight_vec[right_interval] > 0)
				{
					dist_q_p = dist_q_d + (edge.length - edge.interval_point_vec[right_interval]);
					if (dist_q_p <= our_model.bandwidth)
						edge_KA_value += edge.weight_vec[right_interval] * (1 - our_model.gamma*dist_q_p*dist_q_p);
				}

				return edge_KA_value;
			}

			if (left_interval >= right_interval)
			{
				v = (int)floor((edge.length - dist_q_c + dist_q_d) / (2.0*edge.min_dist_diff));

				if (v > 0 && v < edge.num_intervals - 1)
				{
					//left part
					edge_KA_value += (1 - our_model.gamma*dist_q_c*dist_q_c)*edge.agg_weight_c_vec[v - 1]
						- 2 * our_model.gamma*dist_q_c*edge.aug_sum_dist_c_I_deg_1_vec[v - 1]
						- our_model.gamma*edge.aug_sum_dist_c_I_deg_2_vec[v - 1];

					//right part
					edge_KA_value += (1 - our_model.gamma*dist_q_d*dist_q_d)*edge.agg_weight_d_vec[v + 1]
						- 2 * our_model.gamma*dist_q_d*edge.aug_sum_dist_d_I_deg_1_vec[v + 1]
						- our_model.gamma*edge.aug_sum_dist_d_I_deg_2_vec[v + 1];

					//middle part
					if (edge.weight_vec[v] > 0)
					{
						dist_q_p = min(dist_q_c + edge.interval_point_vec[v], dist_q_d + edge.length - edge.interval_point_vec[v]);
						if (dist_q_p <= our_model.bandwidth)
							edge_KA_value += edge.weight_vec[v] * (1 - our_model.gamma*dist_q_p*dist_q_p);
					}
				}

				if (v < 0)
					edge_KA_value += (1 - our_model.gamma*dist_q_d*dist_q_d)*edge.agg_weight_d_vec[0]
					- 2 * our_model.gamma*dist_q_d*edge.aug_sum_dist_d_I_deg_1_vec[0]
					- our_model.gamma*edge.aug_sum_dist_d_I_deg_2_vec[0];
				if (v > edge.num_intervals - 1)
					edge_KA_value += (1 - our_model.gamma*dist_q_c*dist_q_c)*edge.agg_weight_c_vec[edge.num_intervals - 1]
					- 2 * our_model.gamma*dist_q_c*edge.aug_sum_dist_c_I_deg_1_vec[edge.num_intervals - 1]
					- our_model.gamma*edge.aug_sum_dist_c_I_deg_2_vec[edge.num_intervals - 1];

				if (v == 0 || v == edge.num_intervals - 1)
				{
					if (edge.weight_vec[v] > 0)
					{
						dist_q_p = min(dist_q_c + edge.interval_point_vec[v], dist_q_d + edge.length - edge.interval_point_vec[v]);
						if (dist_q_p <= our_model.bandwidth)
							edge_KA_value += edge.weight_vec[v] * (1 - our_model.gamma*dist_q_p*dist_q_p);
					}

					if (v == 0 && edge.num_intervals != 1)
					{
						edge_KA_value += (1 - our_model.gamma*dist_q_d*dist_q_d)*edge.agg_weight_d_vec[1]
							- 2 * our_model.gamma*dist_q_d*edge.aug_sum_dist_d_I_deg_1_vec[1]
							- our_model.gamma*edge.aug_sum_dist_d_I_deg_2_vec[1];
					}

					if (v == edge.num_intervals - 1 && v != 0)
					{
						edge_KA_value += (1 - our_model.gamma*dist_q_c*dist_q_c)*edge.agg_weight_c_vec[edge.num_intervals - 2]
							- 2 * our_model.gamma*dist_q_c*edge.aug_sum_dist_c_I_deg_1_vec[edge.num_intervals - 2]
							- our_model.gamma*edge.aug_sum_dist_c_I_deg_2_vec[edge.num_intervals - 2];
					}
				}
			}
		}

		return edge_KA_value;
	}

	//Aggregate Distance Augmentation (ADA)
	if (our_model.method == 3)
	{
		int left_position, right_position, middle_position;
		double left_value, right_value, middle_value;
		double dist_q_c, dist_q_d;

		//q is in the same edge
		if (cur_l.edge_index == edge_index)
		{
			int q_position;
			double q_value;
			double a_P_left_deg_1, a_P_left_deg_2;
			double a_P_right_deg_1, a_P_right_deg_2;

			left_value = cur_l.dist_n1 - our_model.bandwidth;
			right_value = cur_l.dist_n1 + our_model.bandwidth;
			q_value = cur_l.dist_n1;

			left_position = find_position(left_value, e, 0);
			right_position = find_position(right_value, e, 0);
			q_position = find_position(q_value, e, 0);

			if (right_position == -1) //left, q and right positions are all -1
				return 0;
			
			if (q_position == -1) //left and q positions are all -1
			{
				a_P_right_deg_1 = PS[right_position].aug_sum_dist_c_p_deg_1;
				a_P_right_deg_2 = PS[right_position].aug_sum_dist_c_p_deg_2;

				edge_KA_value = (1 - our_model.gamma*cur_l.dist_n1*cur_l.dist_n1)*(right_position + 1)
					+ 2 * our_model.gamma*cur_l.dist_n1*a_P_right_deg_1 - our_model.gamma*a_P_right_deg_2;

				return edge_KA_value;
			}
				
			if (left_position == -1) //left position is -1
			{
				a_P_left_deg_1 = PS[q_position].aug_sum_dist_c_p_deg_1;
				a_P_left_deg_2 = PS[q_position].aug_sum_dist_c_p_deg_2;
				a_P_right_deg_1 = PS[right_position].aug_sum_dist_c_p_deg_1 - PS[q_position].aug_sum_dist_c_p_deg_1;
				a_P_right_deg_2 = PS[right_position].aug_sum_dist_c_p_deg_2 - PS[q_position].aug_sum_dist_c_p_deg_2;

				edge_KA_value = (1 - our_model.gamma*cur_l.dist_n1*cur_l.dist_n1)*(q_position + 1)
					+ 2 * our_model.gamma*cur_l.dist_n1*a_P_left_deg_1 - our_model.gamma*a_P_left_deg_2
					+ (1 - our_model.gamma*cur_l.dist_n1*cur_l.dist_n1)*(right_position - q_position)
					+ 2 * our_model.gamma*cur_l.dist_n1*a_P_right_deg_1 - our_model.gamma*a_P_right_deg_2;

				return edge_KA_value;
			}

			//left position >= 0
			a_P_left_deg_1 = PS[q_position].aug_sum_dist_c_p_deg_1 - PS[left_position].aug_sum_dist_c_p_deg_1;
			a_P_left_deg_2 = PS[q_position].aug_sum_dist_c_p_deg_2 - PS[left_position].aug_sum_dist_c_p_deg_2;
			a_P_right_deg_1 = PS[right_position].aug_sum_dist_c_p_deg_1 - PS[q_position].aug_sum_dist_c_p_deg_1;
			a_P_right_deg_2 = PS[right_position].aug_sum_dist_c_p_deg_2 - PS[q_position].aug_sum_dist_c_p_deg_2;

			edge_KA_value = (1 - our_model.gamma*cur_l.dist_n1*cur_l.dist_n1)*(q_position - left_position)
				+ 2 * our_model.gamma*cur_l.dist_n1*a_P_left_deg_1 - our_model.gamma*a_P_left_deg_2
				+ (1 - our_model.gamma*cur_l.dist_n1*cur_l.dist_n1)*(right_position - q_position)
				+ 2 * our_model.gamma*cur_l.dist_n1*a_P_right_deg_1 - our_model.gamma*a_P_right_deg_2;

			return edge_KA_value;
		}

		//q is in different edge
		//compute the distances between q to c and q to d 
		dist_q_c = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_c].cur_sp_value,
			cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_c].cur_sp_value);
		dist_q_d = min(cur_l.dist_n1 + our_model.sp_node_vec_node_a[node_index_pi_d].cur_sp_value,
			cur_l.dist_n2 + our_model.sp_node_vec_node_b[node_index_pi_d].cur_sp_value);
		
		//Case 1:
		if (dist_q_c > our_model.bandwidth && dist_q_d > our_model.bandwidth)
			return 0;
		//Case 2:
		if (dist_q_c <= our_model.bandwidth && dist_q_d > our_model.bandwidth)
		{
			left_value = our_model.bandwidth - dist_q_c;
			left_position = find_position(left_value, e, 0);

			if (left_position == -1) //left position is -1
				return 0;
			else
				edge_KA_value = (1 - our_model.gamma*dist_q_c*dist_q_c)*(left_position + 1)
					- 2 * our_model.gamma*dist_q_c*PS[left_position].aug_sum_dist_c_p_deg_1 - our_model.gamma*PS[left_position].aug_sum_dist_c_p_deg_2;
			
			return edge_KA_value;
		}
		//Case 3:
		if (dist_q_c > our_model.bandwidth && dist_q_d <= our_model.bandwidth)
		{
			right_value = our_model.bandwidth - dist_q_d;
			right_position = find_position(right_value, e, 1);
			
			//if (right_position == PS.size() - 1)
			if (right_position == PS.size())
				return 0;
			else
				edge_KA_value = (1 - our_model.gamma*dist_q_d*dist_q_d)*(PS.size() - right_position)
					- 2 * our_model.gamma*dist_q_d*PS[right_position].aug_sum_dist_d_p_deg_1 - our_model.gamma*PS[right_position].aug_sum_dist_d_p_deg_2;

			return edge_KA_value;
		}
		//Case 4:
		if (dist_q_c <= our_model.bandwidth && dist_q_d <= our_model.bandwidth)
		{
			left_value = our_model.bandwidth - dist_q_c;
			left_position = find_position(left_value, e, 0);

			right_value = our_model.bandwidth - dist_q_d;
			right_position = find_position(right_value, e, 1);

			if (left_position == -1 || right_position == PS.size())
			{
				if (left_position == -1 && right_position == PS.size())
					return 0;
				if (right_position == PS.size())
					edge_KA_value = (1 - our_model.gamma*dist_q_c*dist_q_c)*(left_position + 1)
						- 2 * our_model.gamma*dist_q_c*PS[left_position].aug_sum_dist_c_p_deg_1 - our_model.gamma*PS[left_position].aug_sum_dist_c_p_deg_2;
					
				if (left_position == -1)
					edge_KA_value = (1 - our_model.gamma*dist_q_d*dist_q_d)*(PS.size() - right_position)
						- 2 * our_model.gamma*dist_q_d*PS[right_position].aug_sum_dist_d_p_deg_1 - our_model.gamma*PS[right_position].aug_sum_dist_d_p_deg_2;
			}
			else
			{
				if (left_position < right_position)
					edge_KA_value = (1 - our_model.gamma*dist_q_c*dist_q_c)*(left_position + 1)
						- 2 * our_model.gamma*dist_q_c*PS[left_position].aug_sum_dist_c_p_deg_1 - our_model.gamma*PS[left_position].aug_sum_dist_c_p_deg_2
						+ (1 - our_model.gamma*dist_q_d*dist_q_d)*(PS.size() - right_position)
						- 2 * our_model.gamma*dist_q_d*PS[right_position].aug_sum_dist_d_p_deg_1 - our_model.gamma*PS[right_position].aug_sum_dist_d_p_deg_2;
				else
				{
					middle_value = dist_q_d - dist_q_c;
					middle_position = find_position(middle_value, e, 2);

					if (middle_position > -1 && middle_position < (int)PS.size() - 1)
						edge_KA_value = (1 - our_model.gamma*dist_q_c*dist_q_c)*(middle_position + 1)
							- 2 * our_model.gamma*dist_q_c*PS[middle_position].aug_sum_dist_c_p_deg_1 - our_model.gamma*PS[middle_position].aug_sum_dist_c_p_deg_2
							+ (1 - our_model.gamma*dist_q_d*dist_q_d)*(PS.size() - 1 - middle_position)
							- 2 * our_model.gamma*dist_q_d*PS[middle_position + 1].aug_sum_dist_d_p_deg_1 - our_model.gamma*PS[middle_position + 1].aug_sum_dist_d_p_deg_2;
					else
					{
						if (middle_position == -1)
							edge_KA_value = (1 - our_model.gamma*dist_q_d*dist_q_d)*e.PS.size()
								- 2 * our_model.gamma*dist_q_d*PS[0].aug_sum_dist_d_p_deg_1 - our_model.gamma*PS[0].aug_sum_dist_d_p_deg_2;
							
						if (middle_position == PS.size() - 1)
							edge_KA_value = (1 - our_model.gamma*dist_q_c*dist_q_c)*e.PS.size()
								- 2 * our_model.gamma*dist_q_c*PS[PS.size() - 1].aug_sum_dist_c_p_deg_1 - our_model.gamma*PS[PS.size() - 1].aug_sum_dist_c_p_deg_2;							
					}
				}
			}

			return edge_KA_value;
		}
	}

	//basic algorithm for edge_KAF
	if (our_model.method == 1 || our_model.method == 2)
	{
		if (cur_l.edge_index == edge_index)
		{
			for (int i = 0; i < (int)our_model.edge_set[edge_index].PS.size(); i++)
			{
				dist_q_p = fabs(PS[i].dist_n1 - cur_l.dist_n1);
				edge_KA_value += kernel_value(our_model, dist_q_p);
			}

			return edge_KA_value;
		}

		for (int i = 0; i < (int)our_model.edge_set[edge_index].PS.size(); i++)
		{
			if (our_model.method == 1)
				dist_q_p = min(our_model.sp_node_vec[node_index_pi_c].cur_sp_value + PS[i].dist_n1,
					our_model.sp_node_vec[node_index_pi_d].cur_sp_value + PS[i].dist_n2);

			if (our_model.method == 2)
				dist_q_p = min(cur_l.dist_n1 + min(our_model.sp_node_vec_node_a[node_index_pi_c].cur_sp_value + PS[i].dist_n1, our_model.sp_node_vec_node_a[node_index_pi_d].cur_sp_value + PS[i].dist_n2),
					cur_l.dist_n2 + min(our_model.sp_node_vec_node_b[node_index_pi_c].cur_sp_value + PS[i].dist_n1, our_model.sp_node_vec_node_b[node_index_pi_d].cur_sp_value + PS[i].dist_n2));

			edge_KA_value += kernel_value(our_model, dist_q_p);
		}
	}
	
	return edge_KA_value;
}

void NKDV_basic(model& our_model)
{
	Lixel& cur_l = our_model.cur_l;

	cur_l.KDE_value = 0;
	double edge_KA_value;
	int prev_method = -1;
	vector<Point>& PS = our_model.edge_set[our_model.cur_l.edge_index].PS;

	for (int e = 0; e < our_model.m; e++)
	{
		if ((our_model.method == 4 || our_model.method == 5) && our_model.edge_set[e].is_IA == false)
		{
			prev_method = our_model.method;
			if (our_model.method == 4)
				our_model.method = 2;
			if (our_model.method == 5)
				our_model.method = 3;
		}

		edge_KA_value = edge_KAF(our_model, e);
		cur_l.KDE_value += edge_KA_value;

		if (prev_method != -1)
		{
			our_model.method = prev_method;
			prev_method = our_model.method;
		}
	}
}

void augment_preprocess(model& our_model)
{
	double incr_dist_c_p_deg_1;
	double incr_dist_d_p_deg_1;
	double incr_dist_c_p_deg_2;
	double incr_dist_d_p_deg_2;
	int node_c, node_d;
	int size;

	for (int e = 0; e < our_model.m; e++)
	{
		incr_dist_c_p_deg_1 = 0;
		incr_dist_d_p_deg_1 = 0;
		incr_dist_c_p_deg_2 = 0;
		incr_dist_d_p_deg_2 = 0;

		vector<Point>& PS = our_model.edge_set[e].PS;
		node_c = our_model.edge_set[e].n1;
		node_d = our_model.edge_set[e].n2;
		size = PS.size();

		for (int i = 0; i < size; i++)
		{
			incr_dist_c_p_deg_1 += PS[i].dist_n1;
			incr_dist_c_p_deg_2 += PS[i].dist_n1*PS[i].dist_n1;
			PS[i].aug_sum_dist_c_p_deg_1 = incr_dist_c_p_deg_1;
			PS[i].aug_sum_dist_c_p_deg_2 = incr_dist_c_p_deg_2;
			our_model.edge_set[e].dist_n1_vec[i] = PS[i].dist_n1;

			incr_dist_d_p_deg_1 += PS[size - 1 - i].dist_n2;
			incr_dist_d_p_deg_2 += PS[size - 1 - i].dist_n2*PS[size - 1 - i].dist_n2;
			PS[size - 1 - i].aug_sum_dist_d_p_deg_1 = incr_dist_d_p_deg_1;
			PS[size - 1 - i].aug_sum_dist_d_p_deg_2 = incr_dist_d_p_deg_2;
			our_model.edge_set[e].dist_n2_vec[size - 1 - i] = PS[size - 1 - i].dist_n2;
		}

		//obtain dist_diff
		for (int i = 0; i < size; i++)
			our_model.edge_set[e].aug_dist_diff_vec[i]= PS[i].dist_n1 - PS[i].dist_n2;
	}
}

void augment_interval_preprocess(model& our_model)
{
	double prev_dist_n1;
	double min_dist_diff = inf;
	double cur_dist_diff;
	int N_e;
	double i_min, i_max;
	int cur_p;

	#ifdef STATISTICS
		our_model.num_of_edges_w_points = 0;
		our_model.num_of_augmentation = 0;
		double IA_prepro_time;
		auto start_IA_prepro_time = chrono::high_resolution_clock::now();
	#endif

	for (int e = 0; e < our_model.m; e++)
	{
		min_dist_diff = inf;
		our_model.edge_set[e].is_IA = true;
		vector<Point>& PS = our_model.edge_set[e].PS;

		if (PS.size() == 0)
			continue;

		#ifdef STATISTICS
			our_model.num_of_edges_w_points++;
		#endif

		prev_dist_n1 = PS[0].dist_n1;
		for (int p = 1; p < (int)PS.size(); p++)
		{
			cur_dist_diff = PS[p].dist_n1 - prev_dist_n1;
			if (cur_dist_diff < eps) //same value
				continue;

			if (cur_dist_diff < min_dist_diff)
				min_dist_diff = cur_dist_diff;

			prev_dist_n1 = PS[p].dist_n1;
			if (min_dist_diff < const_diff_eps)
			{
				our_model.edge_set[e].is_IA = false;
				break;
			}
		}

		//We do not make interval augmentation for this edge (skip it)
		if (our_model.edge_set[e].is_IA == false)
			continue;

		//number of intervals in edge e using the method IA
		N_e = (int)ceil(our_model.edge_set[e].length / min_dist_diff);
		
		//The method HA (check whether this edge is suitable to use IA)
		if (our_model.method == 5 && our_model.edge_set[e].is_IA == true)
		{
			double cost_ADA;
			double cost_IA;
			double L_e = ceil(our_model.edge_set[e].length / our_model.lixel_reg_length);

			cost_ADA = L_e * log(PS.size());
			cost_IA = L_e + N_e;

			if (cost_ADA < cost_IA)
			{
				our_model.edge_set[e].is_IA = false;
				continue;
			}
		}

		#ifdef STATISTICS
			our_model.num_of_augmentation++;
		#endif

		our_model.edge_set[e].num_intervals = N_e;
		our_model.edge_set[e].min_dist_diff = min_dist_diff;
		our_model.edge_set[e].last_interval_size = our_model.edge_set[e].length - (N_e - 1)*min_dist_diff;
		cur_p = 0;
		for (int i = 0; i < N_e; i++)
		{
			our_model.edge_set[e].aug_sum_dist_c_I_deg_1_vec.push_back(0);
			our_model.edge_set[e].aug_sum_dist_c_I_deg_2_vec.push_back(0);
			our_model.edge_set[e].aug_sum_dist_d_I_deg_1_vec.push_back(0);
			our_model.edge_set[e].aug_sum_dist_d_I_deg_2_vec.push_back(0);
			our_model.edge_set[e].weight_vec.push_back(0);
			our_model.edge_set[e].agg_weight_c_vec.push_back(0);
			our_model.edge_set[e].agg_weight_d_vec.push_back(0);
			our_model.edge_set[e].interval_point_vec.push_back(0);

			i_min = i * min_dist_diff;
			if (i < N_e - 1)
				i_max = i_min + min_dist_diff;
			else
				i_max = our_model.edge_set[e].length;

			if (i > 0)
			{
				our_model.edge_set[e].aug_sum_dist_c_I_deg_1_vec[i] = our_model.edge_set[e].aug_sum_dist_c_I_deg_1_vec[i - 1];
				our_model.edge_set[e].aug_sum_dist_c_I_deg_2_vec[i] = our_model.edge_set[e].aug_sum_dist_c_I_deg_2_vec[i - 1];
				our_model.edge_set[e].agg_weight_c_vec[i] = our_model.edge_set[e].agg_weight_c_vec[i - 1];
			}

			while (cur_p < (int)PS.size())
			{
				if (PS[cur_p].dist_n1 < i_max)
				{
					our_model.edge_set[e].interval_point_vec[i] = PS[cur_p].dist_n1;
					our_model.edge_set[e].aug_sum_dist_c_I_deg_1_vec[i] += PS[cur_p].dist_n1;
					our_model.edge_set[e].aug_sum_dist_c_I_deg_2_vec[i] += PS[cur_p].dist_n1*PS[cur_p].dist_n1;
					our_model.edge_set[e].weight_vec[i]++;
					our_model.edge_set[e].agg_weight_c_vec[i]++;
				}
				else
					break;

				cur_p++;
			}
		}

		cur_p = (int)PS.size() - 1;
		for (int i = N_e - 1; i >= 0; i--)
		{
			i_min = i * min_dist_diff;

			if (i < N_e - 1)
			{
				our_model.edge_set[e].aug_sum_dist_d_I_deg_1_vec[i] = our_model.edge_set[e].aug_sum_dist_d_I_deg_1_vec[i + 1];
				our_model.edge_set[e].aug_sum_dist_d_I_deg_2_vec[i] = our_model.edge_set[e].aug_sum_dist_d_I_deg_2_vec[i + 1];
				our_model.edge_set[e].agg_weight_d_vec[i] = our_model.edge_set[e].agg_weight_d_vec[i + 1];
			}

			while (cur_p >= 0)
			{
				if (PS[cur_p].dist_n2 <= our_model.edge_set[e].length - i_min)
				{
					our_model.edge_set[e].aug_sum_dist_d_I_deg_1_vec[i] += PS[cur_p].dist_n2;
					our_model.edge_set[e].aug_sum_dist_d_I_deg_2_vec[i] += PS[cur_p].dist_n2*PS[cur_p].dist_n2;
					our_model.edge_set[e].agg_weight_d_vec[i]++;
				}
				else
					break;

				cur_p--;
			}
		}

	}

	#ifdef STATISTICS
		auto end_IA_prepro_time = chrono::high_resolution_clock::now();
		IA_prepro_time = (chrono::duration_cast<chrono::nanoseconds>(end_IA_prepro_time - start_IA_prepro_time).count()) / 1000000000.0;
		cout << "IA_preprocessing time: " << IA_prepro_time << endl;
		cout << "number of edges: " << our_model.m << endl;
		cout << "number of edges with points: " << our_model.num_of_edges_w_points << endl;
		cout << "number of augmented edge: " << our_model.num_of_augmentation << endl;
		cout << "number of non-augmented edge: " << our_model.num_of_edges_w_points - our_model.num_of_augmentation << endl;
		exit(0);
	#endif
}