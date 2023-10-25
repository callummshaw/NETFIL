#include "network.h"
#include <stdio.h>
#include <cstdlib>
#include <algorithm>

void region::implement_MDA(int year, mda_strat strat){
    int n_pop = 0;
    int n_treated = 0;
    int n_under_min = 0;

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //for every group
        
        group *grp = j->second;
        
        n_pop += grp->group_pop.size();
        
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //for all people 
            
            double age = k->second->age/365.0; // agents age!

            if(age<strat.min_age) ++n_under_min; 

        }
    }

    double target_prop = 1 - n_under_min /(double)n_pop;

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //for every group
        
        group *grp = j->second;
        
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //for all people 
            
            double age = k->second->age/365.0; // agents age!

            if(random_real() <= strat.Coverage/(double)target_prop){
                ++n_treated;
                k->second->mda(strat.drug);
            }
        }
    }

    number_treated[year] = n_treated;
    achieved_coverage[year] = n_treated/(double)n_pop;
}

void region::handl_commute(int year){
    
    int recalc_years = 5; //how often we want to recalc commuter network
    char distance_type = 'r'; // r for road distance, e for euclidean 

    //firstly need to clear previous storage
    
    if (year % recalc_years == 0) radt_model(distance_type); //generating commuting network

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //now using the network
        group *grp = j->second;
        int no_commute_id = grp->gid;

        //now iterating over all group members
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            agent *cur = k->second; //our agent
            
            if(random_real() > commuting_prop){ //Will not commute!
                grp->day_population.insert(pair<int, agent*>(cur->aid, cur)); //storing them in current group for day population
            } 
            else{ //person will commute
                int commute_id;
                //but commute where?
               for(map<int, double>::iterator j = grp->commuting_pop.begin(); j != grp->commuting_pop.end(); ++j){
                    if(j->second < 1) continue;//every has commuted that we need to this village so we will skip it
                    else{ //person will go here
                        commute_id = j->first;
                        j->second -= 1; //remove member
                        break;
                    }
                }
                groups[commute_id]->day_population.insert(pair<int, agent*>(cur->aid, cur)); //storing them in commuting group for day population
            }
        }
    }
}


