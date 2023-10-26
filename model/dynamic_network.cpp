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
    
    int recalc_years = 2; //how often we want to recalc commuters
    char distance_type = 'r'; // r for road distance, e for euclidean 

    //firstly need to clear previous storage

    if (year % recalc_years == 0){    
        radt_model(distance_type); //generating commuting network

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
                for(map<int, double>::iterator i = grp->commuting_pop.begin(); i != grp->commuting_pop.end(); ++i){
                        if(i->second < 1) continue;//every has commuted that we need to this village so we will skip it
                        else{ //person will go here
                            commute_id = i->first;
                            cur->dgp = groups[commute_id]; //assigning agent to day group
                            i->second -= 1; //remove member
                            break;
                        }
                    }
                    groups[commute_id]->day_population.insert(pair<int, agent*>(cur->aid, cur)); //storing them in commuting group for day population
                }
            }
        }
    }
}

void region::calc_risk(int year, int day, mda_strat strat){
    
    char form = 'l'; //l for limitation, f for facilation, or anything else for linear 

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //setting the force of transmission in each group to 0
        group *grp = j->second;
        grp->day_strength = 0;
        grp->night_strength = 0;
    }
    
    //now looping over all infected agents
    for(map<int, agent*>::iterator j = inf_indiv.begin(); j != inf_indiv.end(); ++j){
        group *dgrp = j->second->dgp; //infected agents daytime group
        group *ngrp = j->second->ngp; //infected agents nightime group

        dgrp->day_strength += (mf_functional_form(form, j->second->worm_strength);
        ngrp->night_strength += mf_functional_form(form, j->second->worm_strength);
    }
}

double region::mf_functional_form(char form, double worm_strength){
    if(form == 'l'){ // limitation
        //TODO write limitation function 
    }
    else if(form == 'f'){ //facilitation
        //TODO write facilitation function 
    }
    else{ //asumme linear
        return worm_strength;
    }
}

