#ifndef network_hpp
#define network_hpp
#include <string>

#include "mda.h"
#include "agent.h"

using namespace std;

class group;                           //groups of people akin to villages
class region;                          //region which is comprised of the groups!

class group{
public:

    int gid;                           //Group ID
    
    double day_strength;            //strength of infection during the day
    double night_strength;           //strength of infection during the day
    
    double lat, lon;                    //latitude & longitude
    region *rgn;                       //region!
    double sum_mf;                      //NEED TO DEFINE

    map<int, agent*> group_pop;       //group population (out of work hours)

    //communiting data
    struct c_node{ //used to store distances to all other groups from current group
        int gid; //other group idea
        double dis; //the distance!
        c_node(int gid, double dis): gid(gid), dis(dis) {}
    };

    vector<c_node*> commuting_dist; //storing the distances
    map<int, double> commuting_pop; //the prop of commuters from each location
    map<int,double> commuting_cumsum; //cumsum of commuters from each location
    map<int, agent*> day_population;    //the commuters to current group! and agents from group that did not commute!

    void add_member(agent *p);
    void rmv_member(agent *p);
    
    void bld_group_pop();  //build initial population
  

    group(int gid, region *rgn, double lat, double lon);
    ~group();
};

class region{
public:
    int rid;                           //region ID 
    string rname;                      //region name                   
    double init_prev;              // initial prevalence
    double init_ratio;
    int rpop;                          //region population
    int next_aid;                      //agent ID tracker for births
    bool init;                         // Has the population been built before?    
    
    double theta1;                      //transmission parameters for the different mf maturation scalings!
    double theta2;
    double theta3;

    double immature_to_antigen;     //at init the ratio of people with just immature worms to antigen positive people (fitted) 
    double immature_and_ant;       //at init the ratio of people with people who are ant pos but mf negative with immature worms

    double worktonot = 0;               //where to a majaortiy of bites occur?

    double mf_to_ant_2014; //used to save fitting data
    
    double agg_param;
    double agg_scale;
    double agg_param_init;

    double age_dist[n_age_groups];     //container for the age distribution
    int age_dist_lower[n_age_groups];
    int age_dist_upper[n_age_groups];
    //used to keep track of total population for easy analysis
    map<int, agent*> pre_indiv;        //collection of immautre worms individuals
    map<int, agent*> inf_indiv;        //collection ofinfectious individuals
    map<int, agent*> uninf_indiv;      //collection of peple with adult worms but are uninfectious individuals (single gender or sterile)

    vector<agent*> pvec[n_age_groups]; //storing all people of certain age group

    //now all the information about the groups
    int next_gid, group_blocks;
    map<int, group*> groups;            //storing all groups in region
    map<string, int> group_names;       //each group assigned name to index
    map<int, string> group_numbers;     //each group assigned number to index
    map<int, double*> group_coords;     //coords of each group

    // For fitting
    double mf_2014 = 0;
    double ant_2014 = 0;

    double mf_2016 = 0;
    double ant_2016 = 0;

    //distances 
    double *euclid_dst;                 //euclidean (L2) distance between groups
    double *road_dst;                   //road (L1) distance between groups

    map<int, int> group_pops;           //pop in each group
 
    double mortality_rate[n_age_groups];                //mortatlity rates
    double birth_rate[n_age_groups];
    double exposure_by_age[16];

    double achieved_coverage[sim_years]; // the actual drug coverage achieved each year (for each year of the simulation). Will be zero for most years.
    int number_treated[sim_years];

    region(int rid, string rname);

    //Functions that run on region
    void sim(int year, mda_strat strategy);                     //wrapper to run simulation
    void handl_commute(int year);                               // generate commuter network and assign
    void remove_agent(agent *p);                                   //remove dead people from population
    void radt_model(char m);                                    //radiation model for daily trips (work/school)
    //void hndl_migrt(int day);                                //TODO long term migration between groups (to help avoid groups that have died out)
    void renew_pop(int year, int day, int dt);
    void hndl_birth(int year, int day, int dt);                         // handle new births
    void calc_risk();         //find prevalence in each village
    void calc_risk_single();
    void update_epi_status(int year, int day, int dt);                  //update agent's epi status
    void seed_lf();                                             //seed LF in population
    void seed_lf_single();

    double mf_functional_form(char form, double worm_strength);            //converts worm strength to mf load

    void implement_MDA(int year, mda_strat strat);           //MDA!
    
    bool pop_reload();
    void read_groups();                                 //read input data
    void bld_groups();                                  //build the model groups 
    void bld_region_population();//build the population of the region
    void read_parameters();

    void reset_population();
    void reset_prev();
    void output_epidemics(int year, int day, mda_strat strategy);    //output outbreak data
    void output_abc_epidemics(int year);
    void output_abc_epidemics_single(int year);
    void output_abc_epidemics_init();
    int factorial(int n);
    vector<double> prob_worms(double prev);
    int number_worms(vector<double> cum_sum_prob, double prob);
};
#endif /* network_hpp */
