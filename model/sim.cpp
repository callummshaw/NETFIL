#include <math.h>
#include "network.h"
#include "mda.h"

extern string prv_out_loc;

void region::sim(int year, mda_strat strat){
    
    bool debug_fit = false; //prints out yearly data

    if (ABC_fitting && ABC_fitting_init){
        cout << "Warning both stages of fitting turned on" << endl;
    }
    //if the first year, must seed LF in the population
    if(year == 0){
        init_prev = 0;
        init_ratio  = 0;
        while(((init_prev < init_prev_min) || (init_prev > init_prev_max)) || ((init_ratio < init_ratio_min) || (init_ratio > init_ratio_max))){
            if (groups.size() > 1){
                seed_lf();
                if( ABC_fitting_init) break;
            }
            else{
                seed_lf_single();
                if( ABC_fitting_init) break;
            }
        }
        if (!ABC_fitting_init){
            cout << "Init prev: " << init_prev << "%" << endl;
            cout << "MF to Ant: " << init_ratio << endl;
        }
    }

    if (prv_out_loc == "print"){
        debug_fit = true;
    }
    
    if(ABC_fitting_init){
        output_abc_epidemics_init();
    }
    else{
        handl_commute(year);
        achieved_coverage[year] = 0;

        if(strat.is_mda_year(year+start_year)){
            implement_MDA(year,strat);
        } 
        if ((ABC_fitting && debug_fit) || (!ABC_fitting && !ABC_fitting_init)){
            output_epidemics(year, strat); 
        }
        else{
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
            
            if (day % epi_dt == theta2*0){
                if(!(inf_indiv.empty() & pre_indiv.empty() & uninf_indiv.empty())) { //If disease has not been eliminated
                    if (groups.size() > 1){
                        calc_risk();
                    } 
                    else{
                        calc_risk_single(); // can speed up process if only one group!
                    }
                    update_epi_status(year, day, epi_dt); //update everyone's LF epi status (including the status of each of their worms)
                }
            }   
        
            if (day % population_dt == 0){
                renew_pop(year, day, population_dt); //deaths
                hndl_birth(year, day, population_dt); //births
            }
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
        double group_prev = 1/(1+exp(-group_log_odds)); //prev in our group
     
        vector<double> cum_sum_prob = prob_worms(group_prev); //found the prob of different worm burdens in village assuming negative bionmial
        //iterating over group members!                        uninf_indiv.insert(pair<int, agent*>(cur->aid, cur));

        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            agent *cur = k->second;
            cur->ngp = grp; //assigning night group

            double inf_status = random_real(); //used to determine number of worms
            int worm_count = number_worms(cum_sum_prob, inf_status); //the number of worms
            
            int wm = 0; //male worms
            int wf = 0; //female worms
            
            if (worm_count > 0){ // person has adult worms
                
                ++ant_pos;

                for(int i = 1; i <= worm_count; ++ i){
                    double mature_period = random_real()*normal(mature_period_mean,mature_period_mean_std);
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
                    double immature_period = random_real()*normal(immature_period_mean,immature_period_mean_std);
                    double mature_period = normal(mature_period_mean,mature_period_mean_std);
                    if(random_real() < proportion_male_agent){ //has one male mature
                        cur->wvec.push_back(new worm('M', immature_period, mature_period ,'M'));
                    }
                    else{ //has one female mature
                        cur->wvec.push_back(new worm('M', immature_period, mature_period ,'F'));
                    }
                }

            }
            else if (random_real() <= group_prev*immature_to_antigen){
                cur->status='E';
                pre_indiv.insert(pair<int, agent*>(cur->aid, cur));
                double immature_period = random_real()*normal(immature_period_mean,immature_period_mean_std);
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
    init_ratio = (double)inf_indiv.size() / (double)ant_pos;
    init_prev = 100 * (double)ant_pos / (double)rpop;
    
}

void region::seed_lf_single(){
    reset_prev();
    double ant_pos = 0;
    //iterating over groups
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        
        group *grp = j->second;
       
        vector<double> cum_sum_prob = prob_worms(ant_0); //found the prob of different worm burdens in village assuming negative bionmial
        //iterating over group members!                        uninf_indiv.insert(pair<int, agent*>(cur->aid, cur));

        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            agent *cur = k->second;
            cur->ngp = grp; //assigning night group

            double inf_status = random_real(); //used to determine number of worms
            int worm_count = number_worms(cum_sum_prob, inf_status); //the number of worms
            
            int wm = 0; //male worms
            int wf = 0; //female worms
            
            if (worm_count > 0){ // person has adult worms
                
                ++ant_pos;

                for(int i = 1; i <= worm_count; ++ i){
                    double mature_period = random_real()*normal(mature_period_mean,mature_period_mean_std);
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
                    double immature_period = random_real()*normal(immature_period_mean,immature_period_mean_std);
                    double mature_period = normal(mature_period_mean,mature_period_mean_std);
                    if(random_real() < proportion_male_agent){ //has one male mature
                        cur->wvec.push_back(new worm('M', immature_period, mature_period ,'M'));
                    }
                    else{ //has one female mature
                        cur->wvec.push_back(new worm('M', immature_period, mature_period ,'F'));
                    }
                }

            }
            else if (random_real() <= ant_0*immature_to_antigen){
                cur->status='E';
                pre_indiv.insert(pair<int, agent*>(cur->aid, cur));
                double immature_period = random_real()*normal(immature_period_mean,immature_period_mean_std);
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
    init_ratio = (double)inf_indiv.size() / (double)ant_pos;
    init_prev = 100 * (double)ant_pos / (double)rpop;
}

vector<double> region::prob_worms(double prev){
    
    int n_worms  = 10; //number of worms we want to consider for init
   
    double alpha = 1 / agg_param_init;
    double worm_mean = agg_param_init*(1 - pow(( 1 - prev),alpha))/pow(( 1 - prev),alpha); //mean worm burden in group (assuming negative binomial)
   
    double worm_prob = 0.0;
    
    vector<double> cum_sum_prob;

    for (int i = 0; i <= n_worms; ++i){ //now working out prob of each worm burden assuming negative binomial
        worm_prob += (tgamma(i + agg_param_init) / (tgamma(agg_param_init) * factorial(i))) * pow((worm_mean / (agg_param_init + worm_mean)),i) * pow((agg_param_init / (agg_param_init + worm_mean)),agg_param_init);

        cum_sum_prob.push_back(worm_prob); // storing
    
    } 
    //Now calculating relative prob of each non-zero worm burden
    return cum_sum_prob;

}

int region::factorial(int n) {
    if (n == 0) return 1;
    return n * factorial(n - 1);
}

int region::number_worms(vector<double> cum_sum_prob, double prob){
    int n_worms = cum_sum_prob.size();

    for ( int i = 0; i <= n_worms; ++i){
        if (prob <= cum_sum_prob[i]) return i;
    }
    
    //cout << "Over 10" << endl;

    return n_worms;
}