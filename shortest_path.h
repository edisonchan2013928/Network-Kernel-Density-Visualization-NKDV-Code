#pragma once
#ifndef SHORTEST_PATH_H
#define SHORTEST_PATH_H

#include "init.h"

typedef priority_queue<sp_node, vector<sp_node>, comparePriority> PQ;

void init_dijkstra(model& our_model, PQ& pq);
void dijkstra(model& our_model);
void copy_sp_info(model& our_model, bool isFirst);

#endif