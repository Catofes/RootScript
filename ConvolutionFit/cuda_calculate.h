#ifndef CUDA_CALCULADE
#define CUDA_CALCULADE

double sub_cuda_normal_calculate(int bins, double min, double max, double x, double mean, double width, double f_min,
                                 double f_max);

void sub_cuda_gaus_prepare(std::vector<double> &x, std::vector<double> &w, int bins);

double sub_cuda_gaus_calculate(int bins, double min, double max, double x, double mean, double width, double f_min,
                               double f_max);

#endif
