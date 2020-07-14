#Compile the program
g++ -c init.cpp -w -o init.o -std=c++11
g++ -c shortest_path.cpp -w -o shortest_path.o -std=c++11
g++ -c KAF.cpp -w -o KAF.o -std=c++11
g++ -c alg_NKDV.cpp -w -o alg_NKDV.o -std=c++11
g++ main.cpp -O3 -o main init.o shortest_path.o KAF.o alg_NKDV.o

exit

#Run the code (Unfortunately, the datasets are not provided, due to the space limitation of Github, we do not provide the datasets)
#our_model.network_fileName = argv[1];
#our_model.out_NKDV_fileName = argv[2];
#our_model.method = atoi(argv[3]);
#our_model.lixel_reg_length = atoi(argv[4]);
#our_model.k_type = atoi(argv[5]);
#our_model.bandwidth = atof(argv[6]);
#our_model.gamma = 1;

time_log_file="./Results/time_log_file.txt"
dir="./datasets/"
result_dir="./Results/"
k_type=2 #Epanechnikov kernel
lixel_reg_length=50 #lixel length (50 means: 50m)
bandwidth=1000 #(1 / gamma = 1000, i.e., gamma = 0.001)

#Seattle
dataset="Seattle"
method=3 #method = 1: Baseline method = 2: SPS method = 3: ADA method = 4: IA method = 5: HA
./main $dir$dataset"/"$dataset"_network" $result_dir$dataset"_m"$method"_k"$k_type $method $lixel_reg_length $k_type $bandwidth | tee -a $time_log_file

#Atlanta
dataset="Atlanta"
method=3
./main $dir$dataset"/"$dataset"_network" $result_dir$dataset"_m"$method"_k"$k_type $method $lixel_reg_length $k_type $bandwidth | tee -a $time_log_file

#San_Francisco
dataset="San_Francisco"
method=3
./main $dir$dataset"/"$dataset"_network" $result_dir$dataset"_m"$method"_k"$k_type $method $lixel_reg_length $k_type $bandwidth | tee -a $time_log_file

#New_York
dataset="New_York"
method=3
./main $dir$dataset"/"$dataset"_network" $result_dir$dataset"_m"$method"_k"$k_type $method $lixel_reg_length $k_type $bandwidth | tee -a $time_log_file