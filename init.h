#pragma once
#ifndef INIT_H
#define INIT_H

#include "Network.h"

const double inf = 999999999999999;
const double eps = 0.00000000001;
const double const_diff_eps = 0.00001; //Used in IA

void load_network(model& our_model);
void obtain_lixel_set(model& our_model);
void init_parameters(int argc, char** argv, model& our_model);
void output_Visual(model& our_model);

#endif