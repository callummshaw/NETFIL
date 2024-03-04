#ifndef headers_h
#define headers_h

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

#define immature_period_mean        30*3 //mean immature period  
#define immature_period_mean_std    30*0.1 //STD dev of immature period 

#define mature_period_mean          364*5 //mean mature period 
#define mature_period_mean_std      364*0.51 //STD dev of mature period 

#define proportion_male_worm        0.5 //proportion of worms that are male
#define proportion_male_agent       0.5 //proportion of agents that are male 

#define n_age_groups                16 //number of age brackets (for seeding pop)

#define sim_years                   7//defining prob functions that are used

#define max_init_age                80 //maximum age of agent upon init

#define init_prev_min               2.75 //minimum initial antigen prev
#define init_prev_max               3.75 //maximum initial antigen prev

#define init_ratio_min              0.15
#define init_ratio_max              0.18

#define ant_0                       0.0325 //initial antigen prev
#define sigma_g                     1.1311//Household standard dev
#define beta_0                      -3.9515//beta_0
       
#define start_year                  2010 //model starting year

#define commuting_prop              0.5 //proportion of group that commut daily (over 5 years old)
#define DailyProbLoseAntigen        0.992327946   //set so the half-life is 90 days i.e. pow(0.5,1/90)

#define ABC_fitting                 true
#define ABC_fitting_init            false

double random_real();
double normal(double mean, double stddev);
int poisson(double rate);
double bite_gamma(double shape, double scale);

#if ABC_fitting || ABC_fitting_init
    #define datadir "LF2/data/"
    #define outdir ""
    #define config "LF2/$config/"
    #define config_pop "LF2/$config/pop/"
    #define Tran_param "TranParams-temp"
#else
    #define datadir "../data/"
    #define outdir "../output/"
    #define config "../$config/"
    #define config_pop "../$config/pop/"
    #define Tran_param "TranParams.csv"
#endif

/*
#define datadir                     "../data/"
#define outdir                      "../output/"
#define config                      "../$config/"
#define config_pop                  "../$config/pop/"
#define Tran_param                  "TranParams.csv"
*/

#define group_name                  "group_names.csv"
#define group_locations             "group_locations.csv"
#define group_populations           "group_populations.csv"
#define exposure_age                "exposure_age.csv"

#define birth_file                  "birth_rates.csv"
#define mortality_file              "mortality_rates.csv"
#define age_brackets                "pop_age_dist.csv"

#define crow_distance               "euc_dist.csv"
#define car_distance                "road_dist.csv"

#define MDA_params                  "MDAParams.csv"
#endif /* headers_h */
