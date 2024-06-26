#include <math.h>
#include "network.h"
#include "mda.h"

extern string prv_out_loc;

void region::sim(int year, mda_strat strat){
    
    bool debug_fit = false; //prints out yearly data


    //if the first year, must seed LF in the population
    if(year == 0){
        init_prev = 0;
        init_ratio  = 0;
        while(((init_prev < init_prev_min) || (init_prev > init_prev_max)) || ((init_ratio < init_ratio_min) || (init_ratio > init_ratio_max))){
            
            seed_lf();
        }
        
        cout << "Init prev: " << init_prev << "%" << endl;
        cout << "MF to Ant: " << init_ratio << endl;
        
    }

    if (prv_out_loc == "print"){
        debug_fit = true;
    }
    
    else{
        handl_commute(year);
        achieved_coverage[year] = 0;

        if (ABC_fitting){
            if (groups.size() == 1){
                output_abc_epidemics_single(year);
            }
            else{
                output_abc_epidemics(year);
            }
        }

        
        int epi_dt = 7; 
        int population_dt = 28;

        for(int day = 0; day < 364; ++day){
            
            if (day % epi_dt == 0){
                if(!(inf_indiv.empty() & pre_indiv.empty() & uninf_indiv.empty())) { //If disease has not been eliminated
                    
                    calc_risk();
                    update_epi_status(year, day, epi_dt); //update everyone's LF epi status (including the status of each of their worms)
                }
            }   
        
            if (day % population_dt == 0){
                renew_pop(year, day, population_dt); //deaths
                hndl_birth(year, day, population_dt); //births
            }

            if((strat.is_mda_year(year+start_year)) && (day == 28)){
                implement_MDA(year,strat);
            } 
        
            if ((day % population_dt == 0) && (!ABC_fitting)){
                output_epidemics(year, day, strat); 
            }
                
        }
    }
}

void region::seed_lf(){
    reset_prev();
    double ant_pos = 0;

    string filename = datadir;
    filename  = filename + "initaggs.csv";
    ifstream file(filename);
    string line;
    vector<double> Init_prev, Init_k;

    while (getline(file, line)) {
        stringstream ss(line);
        string valueA, valueB;
        getline(ss, valueA, ',');
        getline(ss, valueB, ',');
        Init_prev.push_back(stod(valueA));
        Init_k.push_back(stod(valueB));
    }

    file.close();

    //iterating over groups
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        
        double minDifference = 1;
        int closestIndex = -1;

        double ki; //init agg to get correct ant to mf prevalence
        double in_agg; //inverse init agg used in calcs oftern
        double Mean_load; //mean worm burden in group
        double group_prev; //antigen prev in our group
        
        group *grp = j->second;

        if (groups.size() > 1){ //for multiple groups we need to determine antigen prev which is clustered at the group level
            double group_effect = normal(0.0,sigma_g);
            double group_log_odds = group_effect + beta_0;
            group_prev = 1/(1+exp(-group_log_odds));
        
            if (group_prev > 0.48){ //ensuring there is a maximum group prev (rarely hits this but is important to include)
                group_prev = 0.48;
            };
        }else{
            group_prev = ant_0; //for single groups we already know the prev!
        }

        for (int i = 0; i < Init_prev.size(); ++i) { //now given the prev we are finding the init aggregation from values we have calculated previously (finding the roots here was too intensive so I found them previously in mathematica and included them here)
            double difference = abs(Init_prev[i] - group_prev);
            if (difference < minDifference) {
                minDifference = difference;
                closestIndex = i;
            }
        }

        ki = Init_k[closestIndex];
        in_agg = 1.0 / ki;
        Mean_load = ki*(1-pow((1-group_prev),in_agg))/pow((1-group_prev),in_agg);

        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            agent *cur = k->second;
            cur->ngp = grp; //assigning night group
            double agent_running_agg;
            double agent_init_agg; //draw an agents initial bite_scale
            
            if (agg_param < ki){ //region has been seeded with a high init prevalence then we base init on running
                agent_running_agg = bite_gamma(agg_param,1);
                agent_init_agg = agent_running_agg + bite_gamma((ki - agg_param),1);
                agent_running_agg *= agg_scale;
            }else{ //we base running on init!
                agent_init_agg = bite_gamma(ki,1);
                agent_running_agg =  agg_scale*(agent_init_agg + bite_gamma((agg_param-ki),1)); //drawing the agents running scale which is determined by input
            }
        
            cur->bite_scale = agent_running_agg; 

            int worm_count = poisson(in_agg*agent_init_agg*Mean_load); //the number of worms
            
            int wm = 0; //male worms
            int wf = 0; //female worms
            
            if (worm_count > 0){ // person has adult worms
                
                ++ant_pos;

                for(int i = 1; i <= worm_count; ++ i){
                    double mature_period;
                  
                    mature_period = (1-init_beta(1,init_beta_b))*normal(mature_period_mean,mature_period_mean_std);
                    if(random_real() < proportion_male_agent){ //has one male mature
                        cur->wvec.push_back(new worm('M', 0, mature_period ,'M'));
                        ++wm;
                    }
                    else{ //has one female mature
                        cur->wvec.push_back(new worm('M', 0, mature_period ,'F'));
                        ++wf;
                    }
                }
                
                if ((wf > 0) && (wm >0)){ //agent has breeding pair of worms!
                    cur->status='I';
                    cur->worm_strength = wf;
                    inf_indiv.insert(pair<int, agent*>(cur->aid, cur)); //storing the person as infected!
                } 
                else{
                    cur->status='U';
                    uninf_indiv.insert(pair<int, agent*>(cur->aid, cur));
                }

                if(random_real() <= immature_and_ant){ //assigning immature worms to ant positive 
                    int n_worms;
                    n_worms = max(poisson(init_poisson),1);
                
                    for(int i = 1; i <= n_worms; ++ i){ 
                        double immature_period = (1-init_beta(1,init_beta_b))*normal(immature_period_mean,immature_period_mean_std);
                        double mature_period = normal(mature_period_mean,mature_period_mean_std);

                        if(random_real() < proportion_male_worm){
                            cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period ,'M'));
                        }
                        else{
                            cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period ,'F'));
                        }
                    }
                }

            }
            else if (random_real() <= group_prev*immature_to_antigen){
                cur->status='E';
                pre_indiv.insert(pair<int, agent*>(cur->aid, cur));

                int n_worms;
                n_worms = max(poisson(init_poisson),1);
                
                for(int i = 1; i <= n_worms; ++ i){ 
                    double immature_period = (1-init_beta(1,init_beta_b))*normal(immature_period_mean,immature_period_mean_std);
                    double mature_period = normal(mature_period_mean,mature_period_mean_std);

                    if(random_real() < proportion_male_worm){
                        cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period ,'M'));
                    }
                    else{
                        cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period ,'F'));
                    }
                }
            
            }else{
                no_worms_indiv.insert(pair<int, agent*>(cur->aid, cur));
            }
            
        }
    }
    
    init_ratio = (double)inf_indiv.size() / (double)ant_pos;
    init_prev = 100 * (double)ant_pos / (double)rpop;
    
}

int region::factorial(int n) {
    if (n == 0) return 1;
    return n * factorial(n - 1);
}
