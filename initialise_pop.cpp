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

    //reading in village populations
    char *str;
    char *p = std::strtok(str, " ,");


    file = datadir;    file = file + group_populations;
    in.open(file.c_str());
    
    getline(in, line);          //skip header

    while(getline(in,line)){
        str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        p = std::strtok(str, ",");  //village name may have space
        
        //deal with gender record
        int mid = group_names[p];
        
        int pop = 0;

        p = std::strtok(NULL, ", ");    pop = atoi(p);
        
        group_pops.insert(pair<int, int>(mid, pop));
    
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

void group::bld_group_pop(){
    int pp = rgn->group_pops[gid];
    
};

