#include "network.h"
#include <stdio.h>
#include <string.h>
#include <cstring>

using namespace std;

//constructer of region
region::region(int rid, string rname){
    this->rid = rid; //region id
    this->rname = rname;

    //Trackers for IDs
    next_aid = 1;
    next_gid = 1;
    group_blocks = 0;

    for(int i = 0; i < n_age_groups; i++){
        age_dist_lower[i] = i*5;
        age_dist_upper[i] = i*5 + 4;
    };

    //Distances
    euclid_dst = NULL;
    road_dst = NULL;
    
    //Clearing counts region level infection status
    pre_indiv.clear();
    uninf_indiv.clear();
    inf_indiv.clear();

    //recreate population when there multiple simulations
    init = pop_reload();

    if (!init){ //first simulation and we havent built population before
        read_groups();
        bld_groups();
        bld_region_population();
    }

}

bool region::pop_reload(){

    //checking if we have already generated the input!
    string file = config;   file = file + rname;    file = file + ".init";
    ifstream in(file.c_str());
    
    if(!in) return false;
}

void region::read_groups(){
    //function to read in group info!

    ifstream in;
    string line, file;
    
    //firstly reading in group names! 
    file = datadir; file = file + group_name;
    in.open(file.c_str());

    getline(in, line); //header

     while(getline(in, line)){
        if(group_names.find(line) == group_names.end()){
            group_names.insert(pair<string, int>(line, next_gid));
            group_numbers.insert(pair<int, string>(next_gid++, line));
            
        }
    }
    in.close();

    //reading in age distribution
    file = datadir;    file = file + age_brackets;
    in.open(file.c_str());
    
    while(getline(in, line)){
        if(line[0] == '*') continue;
        if(line.length() <= 1) continue;  //empty line with carriage return
        break;
    }
    int ii = 0;

    while(getline(in, line)){
        age_dist[ii] = atof(line.c_str());
        ii++;
    }
    in.close();

    char *str;
    char *p = std::strtok(str, " ,");

    //reading in village populations
    file = datadir;    file = file + group_populations;
    in.open(file.c_str());
    
    getline(in, line);          //skip header

    while(getline(in,line)){
        str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        p = std::strtok(str, ",");  //village name may have space
        
        //deal with gender record
        int gid = group_names[p];
        
        int pop = 0;

        p = std::strtok(NULL, ", ");    pop = atoi(p);
        
        group_pops.insert(pair<int, int>(gid, pop));
    
        delete []str;
    }
    in.close();


    //now reading in village locations
    file = datadir;     file = file + group_locations;
    in.open(file.c_str());
    
    getline(in, line); //header
    
    while(getline(in, line)){
        str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        
        p = std::strtok(str, ",");  //village name may have space
        map<string, int>::iterator k = group_names.find(p);
        if(k == group_names.end()) continue;
        
        int mid = k->second;
        p = std::strtok(NULL, ",");     double lat = atof(p);
        p = std::strtok(NULL, ",");     double log = atof(p);
        
        double *r = new double[2];
        r[0] = lat;     r[1] = log;
        group_coords.insert(pair<int, double*>(mid, r));
    }
    in.close();
    
    if(group_coords.size() < group_names.size()){
        cout << "Group coordinates are missing" << endl;
        exit(1);
    }

    //calculate cpop
    for(map<int, int>::iterator j = group_pops.begin(); j != group_pops.end(); ++j){
        int gid = j->first;
        rpop += group_pops[gid];
    }
}

void region::bld_groups(){
    for(map<int, double*>::iterator j = group_coords.begin(); j != group_coords.end(); ++j){
        int gid = j->first;
        double lat = j->second[0], log = j->second[1];
        groups.insert(pair<int, group*>(gid, new group(gid, this, lat, log)));
    }
};

void region::bld_region_population(){

    //going village by village
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        group *grp = j->second;
        grp->bld_group_pop();
    }

};


//constructer of groups
group::group(int gid, region *rgn, double lat, double lon){
    this->gid = gid;
    this->rgn = rgn;
    this->lat = lat;
    this->lon = lon;

    this->sum_mf = 0;
};

group::~group(){
    rgn = NULL;

    for(map<int, agent*>::iterator j = group_pop.begin(); j != group_pop.end(); ++j)
        delete j->second;
    group_pop.clear();
};

//build individual group populations!
void group::bld_group_pop(){

    int group_population = rgn->group_pops[gid];
    double* age_dist = rgn->age_dist;
    
    for(int i = 0; i< n_age_groups; i++){
        if (age_dist[i] == 0) continue;

        int lower_bound = rgn->age_dist_lower[i]; //lower bound of our age bracket
        int upper_bound = rgn->age_dist_upper[i]; //upper_band of our age bracket

        int pp = group_population*age_dist[i];

        while(pp-- > 0){
            int id = rgn->next_aid++;
            int age = 52*(lower_bound + (upper_bound - lower_bound)*random_real()); // age in weeks
            agent *p = new agent(id,age); //creating new agent of correct age!
            add_member(p);
        }

    };
};

void group::add_member(agent *p){
    group_pop.insert(pair<int, agent*>(p->aid, p));
};

