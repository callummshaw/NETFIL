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

#define n_age_groups                16 //number of age brackets (for seeding pop)

#define sim_years                   26//defining prob functions that are used

#define max_init_age                80*52 //maximum age of agent upon init

#define init_prev_min               2.75 //minimum initial antigen prev
#define init_prev_max                3.75 //maximum initial antigen prev

#define sigma_g                     1.31//Household standard dev
#define beta_0                      -3.9515//beta_0

#define start_year                  2010 //model starting year

#define commuting_prop              0.5 //proportion of group that commut daily (over 5 years old)
#define DailyProbLoseAntigen        0.992327946   //set so the half-life is 90 days i.e. pow(0.5,1/90)
double random_real();
double normal(double mean, double stddev);
int poisson(double rate);

#define datadir                     "../data/"
#define outdir                      "../output/"
#define config                      "../$config/"
#define config_pop                  "../$config/pop/"

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
