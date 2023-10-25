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
        implement_MDA(year,strat);
    }

}

void region::seed_lf(){
    
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

            double r = random_real();

            if(r <= group_prev){ //agent is antigen positive!
            ++ant_pos;
            
            //NEED TO WRITE                
            }
        }
    }

    init_prev = 100 * ant_pos / rpop;
    
    cout << "Init prev: " << init_prev << "%" << endl;
}