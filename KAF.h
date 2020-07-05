#pragma once
#ifndef KAF_H
#define KAF_H

#include "init.h"

double kernel_value(model& our_model, double dist);
double edge_KAF(model& our_model, int edge_index);
void NKDV_basic(model& our_model);

//Used in aggregate distance augmentation (ADA)
void augment_preprocess(model& our_model);
//Used in interval augmentation (IA)
void augment_interval_preprocess(model& our_model);

#endif