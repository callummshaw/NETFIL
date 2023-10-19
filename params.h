#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <set>
#include <map>

using namespace std;


#define work_bite_rate              0.2 //proportion of bites in working hours
#define offwork_bite_rate           0.8 //proportion of bites in the rest of the day

#define immature_period_mean        36 //mean immature period in weeks
#define immature_period_mean_std    2 //STD dev of immature period in weeks

#define mature_period_mean          312 //mean mature period in weeks
#define mature_period_mean_std      10 //STD dev of mature period in weeks

#define proportion_male             0.5 //proportion of worms that are male

//defining prob functions that are used
double random_real();
double normal(double mean, double stddev);
int poisson(double rate);
