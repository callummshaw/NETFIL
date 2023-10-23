#ifndef network_hpp
#define network_hpp
#include <string>
#include "agent.h"
#include "mda.h"
using namespace std;

class group;                           //groups of people akin to villages
class region;                          //region which is comprised of the groups!

class group{
public:
    int gid;                           //Group ID
    int MDA = 0;                       //Time since MDA
    int days_before_MDA = 0;           //Time till MDA
    double lat, lon;                    //latitude & longitude
    region *rgn;                       //region!
    double sum_mf;                      //NEED TO DEFINE

    map<int, agent*> group_pop;       //group population
 

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
    double init_prev = 0;              // initial prevalence
    int rpop;                          //region population
    int next_aid;                      //agent ID tracker for births
    bool init;                         // Has the population been built before?
    
    double age_dist[n_age_groups];     //container for the age distribution

    int age_dist_lower[n_age_groups];
    int age_dist_upper[n_age_groups];
    //used to keep track of total population for easy analysis
    map<int, agent*> pre_indiv;        //collection of latent individuals
    map<int, agent*> inf_indiv;        //collection ofinfectious individuals
    map<int, agent*> uninf_indiv;      //collection of peple with adult worms but are uninfectious individuals (single gender or sterile)

    vector<agent*> pvec[n_age_groups]; //storing all people of certain age group

    //now all the information about the groups
    int next_gid, group_blocks;
    map<int, group*> groups;            //storing all groups in region
    map<string, int> group_names;       //each group assigned name to index
    map<int, string> group_numbers;     //each group assigned number to index
    map<int, double*> group_coords;     //coords of each group

    //distances 
    double *euclid_dst;                 //euclidean (L2) distance between groups
    double *road_dst;                   //road (L1) distance between groups

    map<int, int> group_pops;           //pop in each group
 
    double mortality_rate[n_age_groups];                //mortatlity rates
    double birth_rate[n_age_groups];
    double exposure_by_age[16];

    region(int rid, string rname);

    //Functions that run on region
    void sim(int year, mda_strat strategy);              //wrapper to run simulation
    void rmv_agent(agent *p);                           //remove dead people from population
    void radt_model(char m);                            //radiation model for daily trips (work/school)
    void hndl_migrt(int week);                          //long term migration between regions
    void hndl_birth(int week);                          // handle new births
    void calc_risk(int week, mda_strat strat);          //find prevalence in each village
    void update_epi_status(int week);                   //update agent's epi status
    void seed_clustered_epidemics();                    //seed LF in population
    void implement_MDA(int week, mda_strat strat);      //MDA!
    
    double achieved_coverage[sim_years];                //MDA coverage!
    bool pop_reload();
    void read_groups();                                 //read input data
    void bld_groups();                                  //build the model groups 
    void bld_region_population();//build the population of the region
    void read_parameters();

    void reset_population();
    void reset_prevalence();

};
#endif /* network_hpp */
