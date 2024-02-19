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
