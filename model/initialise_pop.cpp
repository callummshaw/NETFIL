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

        //now saving all our config files (to use on other runs!)

        //Saving group meta data
        string file = config;   file = file + rname;    file = file + ".init";
        ofstream out(file.c_str());

        out << std::setprecision(2) << std::setiosflags(std::ios::fixed);
        
        out << "TRUE" << endl;
        out << rpop << endl; //region population
        out << next_aid << endl; //next agent id
        out << next_gid << endl; //next group id
        out << group_blocks << endl; //number of groups

        //saving village names/numbers and locations
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << j->first << "," << group_numbers[j->first] << "," << j->second->lon << "," << j->second->lat << endl;
        }
        out.close();

        //now saving group data, group by group! 
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            group *grp = j->second; //individual group data

            string grp_str = group_numbers[grp->gid]; //group name

            //file name 
            file = config_pop;  file = file + grp_str;  file = file + "_pop.init";

            out.open(file.c_str());
            out << "ID,age" << endl;

            //iterating over agents!
            for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                agent *cur = k->second;
                out << cur->aid << "," << cur->age << endl;
            }
            out.close();
        }
    }

    read_parameters();

}

bool region::pop_reload(){

    //checking if we have already generated the input!
    string file = config;   file = file + rname;    file = file + ".init";
    ifstream in(file.c_str());
    
    if(!in) return false;

    string line;
    getline(in, line);
    if(line != "TRUE") return false;
    
    getline(in, line);      rpop = atoi(line.c_str());
    getline(in, line);      next_aid = atoi(line.c_str());
    getline(in, line);      next_gid = atoi(line.c_str());
    getline(in, line);      group_blocks = atoi(line.c_str());
    
    while(getline(in, line)){
        char *str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        
        char *p = std::strtok(str, ",");        int id = atoi(p);
        p = std::strtok(NULL, ",");             string grp = p;
        p = std::strtok(NULL, ",");             double lon = atof(p);
        p = std::strtok(NULL, ",");             double lat = atof(p);
        
        group_names.insert(pair<string, int>(grp, id));
        group_numbers.insert(pair<int, string>(id, grp));
        groups.insert(pair<int, group*>(id, new group(id, this, lat, lon)));
        
        delete []str;
    }
    in.close();
    
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        string grp_str = group_numbers[j->first];
        
        file = config_pop;   file = file + grp_str;   file = file + "_pop.init";
        in.open(file.c_str());
        
        string line;
        getline(in, line);
        while(getline(in, line)){
            char *str = new char[line.size()+1];
            std::strcpy(str, line.c_str());
            
            char *p = std::strtok(str, ",");        int id = atoi(p);
            p = std::strtok(NULL, ",");             int age = atoi(p);
            
            agent *pp = new agent(id, age);

            j->second->add_member(pp);
            
            delete []str;
        }
        in.close();
    }

    return true;
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
    
    group_blocks = (int)group_names.size();

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
    char *p = NULL;

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
}

void region::bld_region_population(){

    //going village by village
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        group *grp = j->second;
        grp->bld_group_pop();
    }
}

void region::read_parameters(){
    
    ifstream in;
    string line, file;

    //reading in birth rate data 
    file = datadir;     file = file + birth_file;
    in.open(file.c_str());

    //skip the description
    while(getline(in, line)){
        if(line[0] == '*') continue;
        if(line.length() <= 1) continue;  //empty line with carriage return
        break;
    }

    int ii = 0;

    while(getline(in, line)){
        birth_rate[ii] = atof(line.c_str());
        ii++;
    }

    in.close();


     //reading in mortality rate data 
    file = datadir;     file = file + mortality_file;
    in.open(file.c_str());

    //skip the description
    while(getline(in, line)){
        if(line[0] == '*') continue;
        if(line.length() <= 1) continue;  //empty line with carriage return
        break;
    }

    ii = 0;
    while(getline(in, line)){
        mortality_rate[ii] = atof(line.c_str());
        ii++;
    }
    in.close();

    //read exposure by age
    file = datadir;    file = file + exposure_age;
    in.open(file.c_str());
    if(!in){
        cout << "open " << file << " failed" << endl;
        exit(1);
    }
    getline(in, line);
    while(getline(in,line)){
        char *str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        char *p = NULL;
        p = std::strtok(str, ",");      int age = atoi(p);
        p = std::strtok(NULL, ",");
        p = std::strtok(NULL, ",");     double expo = atof(p);
        
        exposure_by_age[age] = expo;
        delete []str;
    }
    in.close();

    //reading road_dst & euclid_dst
    int len = group_blocks*(group_blocks-1)/2;
    road_dst = new double[len];     memset(road_dst, 0, sizeof(double)*len);
    euclid_dst = new double[len];   memset(euclid_dst, 0, sizeof(double)*len);
    
    file = datadir;     file = file + car_distance;
    in.open(file.c_str());
    
    getline(in, line);
    vector<string> grp_vec;
    {
        char *str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        
        char *p = std::strtok(str, ",");
        while(p != NULL){
            grp_vec.push_back(p);
            p = std::strtok(NULL, ",");
        }
        delete []str;
        
        while(getline(in, line)){
            str = new char[line.size()+1];
            std::strcpy(str, line.c_str());
            
            p = std::strtok(str, ",");
            string src = p;
            int src_id = group_names[src];
            
            int index = 0;
            p = std::strtok(NULL, ",");
            while(p != NULL){
                string tag = grp_vec[index++];
                int tag_id = group_names[tag];
                
                if(tag_id > src_id){
                    double dd = atof(p);
                    
                    int ii = (src_id-1)*(group_blocks*2-src_id)/2 + tag_id-src_id - 1;
                    road_dst[ii] = dd;
                }
                
                p = std::strtok(NULL, ",");
            }
            delete []str;
        }
    }
    grp_vec.clear();
    grp_vec.shrink_to_fit();
    in.close();
    
    file = datadir;     file = file + crow_distance;
    in.open(file.c_str());
    
    getline(in, line);
    {
        char *str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        
        char *p = std::strtok(str, ",");
        while(p != NULL){
            grp_vec.push_back(p);
            p = std::strtok(NULL, ",");
        }
        delete []str;
        
        while(getline(in, line)){
            str = new char[line.size()+1];
            std::strcpy(str, line.c_str());
            
            p = std::strtok(str, ",");
            string src = p;
            int src_id = group_names[src];
            
            int index = 0;
            p = std::strtok(NULL, ",");
            while(p != NULL){
                string tag = grp_vec[index++];
                int tag_id = group_names[tag];
                
                if(tag_id > src_id){
                    double dd = atof(p);
                    
                    int ii = (src_id-1)*(group_blocks*2-src_id)/2 + tag_id-src_id - 1;
                    euclid_dst[ii] = dd;
                }
                
                p = std::strtok(NULL, ",");
            }
            delete []str;
        }
    }
    in.close();
    


}

void region::reset_population(){
   
    //resetting population
    pre_indiv.clear();
    inf_indiv.clear();
    uninf_indiv.clear();

    for(map<int, group*>::iterator j = groups.begin();  j != groups.end(); ++j){ //iterating through groups
        delete j->second;
    }
    groups.clear();

    for(map<int, double*>::iterator j = group_coords.begin();  j != group_coords.end(); ++j){ //iterating through groups
        delete [] j->second;
    }
    group_coords.clear();

    group_names.clear();
    group_numbers.clear();

    group_pops.clear();

    rpop = 0;
    next_aid = 1;
    next_gid = 1;
    group_blocks = 0;

    if(!pop_reload()){
        cout << "reload pop err" << endl;
        exit(1);
    }
    read_parameters();

}

//constructer of groups
group::group(int gid, region *rgn, double lat, double lon){
    this->gid = gid;
    this->rgn = rgn;
    this->lat = lat;
    this->lon = lon;

    this->sum_mf = 0;
}

group::~group(){
    rgn = NULL;

    for(map<int, agent*>::iterator j = group_pop.begin(); j != group_pop.end(); ++j)
        delete j->second;
    group_pop.clear();
}

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
            int age = 365*(lower_bound + (upper_bound - lower_bound)*random_real()); // age in weeks
            agent *p = new agent(id,age); //creating new agent of correct age!
            add_member(p);
        }

    }
}

void group::add_member(agent *p){
    group_pop.insert(pair<int, agent*>(p->aid, p));
}
