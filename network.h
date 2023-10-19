#ifndef network_hpp
#define network_hpp
#include <string>
#include "agent.h"

using namespace std;

class group; //groups of people akin to villages
class region; //region which is comprised of the groups!

class region{
public:
    int rid;                //region ID 
    string rname;           //region name                   
    double init_prev = 0;   // initial prevalence
    int rpop;               //region population
    int next_aid;           //agent ID tracker for births
    
    //used to keep track of total population for easy analysis
    map<int, agent*> pre_indiv;        //collection of latent individuals
    map<int, agent*> inf_indiv;        //collection ofinfectious individuals
    map<int, agent*> uninf_indiv;      //collection of peple with adult worms but are uninfectious individuals (single gender or sterile)

    //now all the information about the groups
    map<int, group*> groups;            //storing all groups in region
    map<string, int> group_numbers;     //each group assigned number to index
    map<int, string> group_names;       //storage of all group names
    map<int, double*> group_coords;     //coords of each group

    //distances 
    double *euclid_dst;                 //euclidean (L2) distance between groups
    double *road_dst;                   //road (L1) distance between groups

    map<int, int> group_mpops;          //male pop in each group
    map<int, int> group_fpops;          //female pop in each group
};
#endif /* network_hpp */
