#include "network.h"
#include "mda.h"

void region::sim(int year, mda_strat strat){
    
    //if the first year, must seed LF in the population
    if(year == 0){
        init_prev = 0;

        while((init_prev < init_prev_min) || (init_prev > init_prev_max)){
            seed_lf();
        }
    }

    handl_commute(year);
    achieved_coverage[year] = 0;

    if(strat.is_mda_year(year+start_year)){
        //implement_MDA(year,strat);
    }
    output_epidemics(year, strat); 
    
    int dt = 1; 

    for(int day = 0; day < 364; ++day){
        //if (day % dt == 0){
            if(!(inf_indiv.empty() & pre_indiv.empty() & uninf_indiv.empty())) { //If disease has not been eliminated
            
                calc_risk(year, day, dt, strat); //Determine who gets infected with new worms today - doesn't update epi status
            
                update_epi_status(year, day, dt); //update everyone's LF epi status (including the status of each of their worms)
                
          //  }
        //if (day % dt == 0){
            renew_pop(year, day, dt); //deaths
            hndl_birth(year, day); //births
        }
    }
}

void region::seed_lf(){
    reset_prev();
    double ant_pos = 0;
    
    //iterating over groups
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        
        group *grp = j->second;

        double group_effect = normal(0.0,sigma_g);
        double group_log_odds = group_effect + beta_0;
        double group_prev = 1/(1+exp(-group_log_odds));

        //iterating over group members!
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){

            agent *cur = k->second;
            cur->ngp = grp; //assigning night group
           
            if(random_real() <= group_prev){ //agent is antigen positive!
                ++ant_pos;

                if(random_real() <= matedtonot){//has mated worms!
                //eg if we assume everyone has 1 pair!
                    cur->status='I';
                    cur->worm_strength = 1.0;
                    double mature_period = normal(mature_period_mean,mature_period_mean_std);
                    cur->wvec.push_back(new worm('M', 0, mature_period ,'M'));
                    cur->wvec.push_back(new worm('M', 0, mature_period ,'F'));
                    inf_indiv.insert(pair<int, agent*>(cur->aid, cur));

                } 
                else{ //antigen postive but does not have mated worms!
                    double rr = random_real();
                    if (rr <= 0.666){ // has one mature worm
                        cur->status='U';
                        double mature_period = normal(mature_period_mean,mature_period_mean_std);
                        if(random_real() < proportion_male_agent){ //has one male mature
                            cur->wvec.push_back(new worm('M', 0, mature_period ,'M'));
                            uninf_indiv.insert(pair<int, agent*>(cur->aid, cur));
                        }
                        else{ //has one female mature
                            cur->wvec.push_back(new worm('M', 0, mature_period ,'F'));
                            uninf_indiv.insert(pair<int, agent*>(cur->aid, cur));
                        }
                    }
                    else{ //has immature worms! for now will assume only 1 of either gender
                        cur->status='E';
                        double immature_period = normal(immature_period_mean,immature_period_mean_std);
                        double mature_period = normal(mature_period_mean,mature_period_mean_std);
                        if(random_real() < proportion_male_worm){
                            cur->wvec.push_back(new worm('P', immature_period, mature_period ,'M'));
                            pre_indiv.insert(pair<int, agent*>(cur->aid, cur));
                        }
                        else{
                            cur->wvec.push_back(new worm('P', immature_period, mature_period ,'F'));
                            pre_indiv.insert(pair<int, agent*>(cur->aid, cur));
                        }
                    }
                }
            
            //TODO finish worm seeding  

            }
        }
    }

    init_prev = 100 * ant_pos / rpop;
    
    cout << "Init prev: " << init_prev << "%" << endl;
}