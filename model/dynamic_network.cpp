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
                double cum_sum_floor = 0;
                if(random_real() > commuting_prop){ //Will not commute!
                    grp->day_population.insert(pair<int, agent*>(cur->aid, cur)); //storing them in current group for day population
                    cur->dgp = groups[no_commute_id];
                } 
                else{ //person will commute
                    int commute_id;
                    double commute_dest = random_real();
                    //but commute where?
                    for(map<int, double>::iterator i = grp->commuting_cumsum.begin(); i != grp->commuting_pop.end(); ++i){
                        
                        if ((cum_sum_floor < commute_dest) && (commute_dest <= i->second)){
                            commute_id = i->first;
                            cur->dgp = groups[commute_id]; //assigning agent to day group
                            groups[commute_id]->day_population.insert(pair<int, agent*>(cur->aid, cur)); //storing them in commuting group for day population
                            goto found_commute;
                        } 
                        else {
                            cum_sum_floor = i->second;
                        }   
                    }
                }
                found_commute:;
            }
        }
    }
    rpop = 0;
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        group *grp = j->second;   
        rpop += grp->group_pop.size();
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
        
        agent *cur =j->second;
        group *dgrp = cur->dgp; //infected agents daytime group
        group *ngrp = cur->ngp; //infected agents nightime group

        int age = int(cur->age / 365);
        double c = 1;

        if(age <= 15) c = exposure_by_age[age];
             
        dgrp->day_strength += mf_functional_form(form, j->second->worm_strength) / (double) dgrp->day_population.size();
        ngrp->night_strength += mf_functional_form(form, j->second->worm_strength) / (double) ngrp->day_population.size();
    }

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //looping over groups
        group *grp = j->second;
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //night infections
            agent *cur = k->second; //our person
            char prev_status = cur ->status;
            double c = 1; //
            int age = int(cur->age / 365);
            if (age <= 15) c = exposure_by_age[age];
            
            cur->sim_bites(grp->night_strength,'N', c); // simulating the bites!

            if(cur->status == 'E' && prev_status == 'S'){
                pre_indiv.insert(pair<int, agent *>(cur->aid, cur));
            }
        }

        for(map<int, agent*>::iterator k = grp->day_population.begin(); k != grp->day_population.end(); ++k){ //day infections

            agent *cur = k->second; //our person
            char prev_status = cur ->status;
            double c = 1; //
            int age = int(cur->age / 365);
            if (age <= 15) c = exposure_by_age[age];
            
            cur->sim_bites(grp->day_strength,'D', c); // simulating the bites!

            if(cur->status == 'E' && prev_status == 'S'){
                pre_indiv.insert(pair<int, agent *>(cur->aid, cur));
            }

        }
    }
}

double region::mf_functional_form(char form, double worm_strength){
    if(form == 'l'){ // limitation
        //TODO write limitation function 
        return worm_strength;
    }
    else if(form == 'f'){ //facilitation
        //TODO write facilitation function 
        return worm_strength;
    }
    else{ //asumme linear
        return worm_strength;
    }
}

void region::update_epi_status(int year, int day){

    for(map<int, agent*>::iterator j = pre_indiv.begin(); j != pre_indiv.end();){ //looking at all agents with mature worms but not a set!
        
        agent *p = j->second;

        p->update(year,day);

        if(p->status != 'E'){ //agent has left this stage of infection
            p->ChangedEpiToday = true;
            pre_indiv.erase(j++);

            if(p->status == 'I'){ //infective now!
                inf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
            else if(p->status == 'U'){ //have immature worms only npw
                uninf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
        }
        else ++j;
    }

    for(map<int, agent*>::iterator j = uninf_indiv.begin(); j != pre_indiv.end();){ //looking at all agents with only immature worms!
        
        agent *p = j->second;

        if(!p->ChangedEpiToday) p->update(year, day);
        else p->ChangedEpiToday = false;
    
        if(p->status != 'I'){
            uninf_indiv.erase(j++);
            if(p->status == 'I'){
                p->ChangedEpiToday = true;
                inf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
            else if(p->status == 'E'){
                pre_indiv.insert(pair<int, agent*>(p->aid, p));
            }
        }
        else ++j;
    }

    for(map<int, agent*>::iterator j = inf_indiv.begin(); j != inf_indiv.end();){
        agent *p = j->second;
        if(!p->ChangedEpiToday) p->update(year, day);
        else p->ChangedEpiToday = false;
        if(p->status != 'I'){
            inf_indiv.erase(j++);
            if(p->status == 'E'){
                pre_indiv.insert(pair<int, agent*>(p->aid, p));
            }
            else if(p->status == 'U'){
                uninf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
        }
        else ++j;
    }

}

void region::renew_pop(int year, int day){
    //handleing deaths!
    vector<agent*> deaths;

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
        group *grp = j->second;

        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //going through group members
            agent *cur = k->second;
            int index = int(int(cur->age/365)/5);
            if(index > 15) index = 15; //all 75+ the same

            double prob = 1 - exp(-mortality_rate[index]);
            if(random_real() < prob) deaths.push_back(cur); //seeing if agent dies depending on age
            else ++cur->age; //increase everyones age
        }
    }

    while(deaths.size() > 0){ //now removing agents that have died
        agent *cur = deaths.back();
        remove_agent(cur);
        deaths.pop_back();
    }       
}

void region::remove_agent(agent *p){

    //removing agent from lists of infected
   
    if(p->status == 'E') pre_indiv.erase(p->aid);
    else if(p->status == 'I') inf_indiv.erase(p->aid);
    else if(p->status == 'U') uninf_indiv.erase(p->aid);
    
    //daytime group
    group *dgrp = p->dgp;
    //nightime group
    group *ngrp = p->ngp;

    ngrp->group_pop.erase(p->aid);
    dgrp->day_population.erase(p->aid);
    
   
    
    delete p;
}

void region::hndl_birth(int year, int day){ //deal with births
    
    int total_births  = 0;

    for(map<int,group*>::iterator j = groups.begin(); j != groups.end(); j++){//looping over groups
        group *grp = j->second;
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){//over agents
            agent *cur = k->second;
            if(cur->age >= 18*365 && cur->age < 60*365){
                if(random_real() < female_prop){ // excluding male population
                    int index = int((int(cur->age/365)-15)/5);
                    double prob = 1 - exp(-birth_rate[index]);
                    if(random_real() < prob) ++total_births;
                }
            }
        }
        //now assigning births
        while (total_births > 0) {
            agent *bb = new agent(next_aid++, 0); //have birth!
            bb->ngp = grp;
            bb->dgp = grp;
            grp->add_member(bb); //assigning baby to group
            grp->day_population.insert(pair<int, agent*>(bb->aid, bb));//assume baby stays within group during day
            
            --total_births;
        }
    }
}