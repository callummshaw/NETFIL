#include "params.h"
#include "rng.h"

using namespace std;

double random_real(){
    uniform_real_distribution<> distribution(0.0,1.0);

    return distribution(gen);
}

int poisson(double rate){
    poisson_distribution<int> distribution(rate);

    return distribution(gen);
}

double normal(double mean, double stddev){
    normal_distribution<double> distribution(mean, stddev);

    return distribution(gen);
}

double bite_gamma(double shape, double scale){
    gamma_distribution<double> distribution(shape, scale);

    return distribution(gen);
}

double init_beta(double a, double b){
    
    gamma_distribution<> X(a, 1.0);
    gamma_distribution<> Y(b, 1.0);

    double x = X(gen);
    double y = Y(gen);

    return  x / (x + y);
}

void partial_shuffle(vector<double>& vec, int start, int end){
   shuffle(vec.begin() + start, vec.begin() + end, gen);
}