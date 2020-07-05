#include "init.h"
#include "alg_NKDV.h"

int main(int argc, char**argv)
{
	model our_model;
	init_parameters(argc, argv, our_model);
	load_network(our_model);
	obtain_lixel_set(our_model);
	NKDV_algorithm(our_model);
	output_Visual(our_model);
}