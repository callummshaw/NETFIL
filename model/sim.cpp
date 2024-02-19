#include <math.h>
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

    if (!ABC_fitting){
        output_epidemics(year, strat); 
    }
    else if(ABC_fitting){
        output_abc_epidemics(year);
    }

    int epi_dt = 7; 
    int population_dt = 28;

    for(int day = 0; day < 364; ++day){
        
        if (day % epi_dt == 0){
            if(!(inf_indiv.empty() & pre_indiv.empty() & uninf_indiv.empty())) { //If disease has not been eliminated
            
                calc_risk(); //Determine who gets infected with new worms today - doesn't update epi status
            
                update_epi_status(year, day, epi_dt); //update everyone's LF epi status (including the status of each of their worms)
                
            }
        }   
        
        if (day % population_dt == 0){
            renew_pop(year, day, population_dt); //deaths
            hndl_birth(year, day, population_dt); //births
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
     
        int n_worms = 5; //max  number of worms we want to consider

        vector<double> cum_sum_prob = prob_worms(group_prev, n_worms); //found the prob of different worm burdens in village assuming negative bionmial
        //iterating over group members!                        uninf_indiv.insert(pair<int, agent*>(cur->aid, cur));

        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            agent *cur = k->second;
            cur->ngp = grp; //assigning night group
          
            if(random_real() <= group_prev){ //agent is antigen positive!
                ++ant_pos;
                
                if(random_real() <= matedtonot){//has mated worms! How many mated worms?
                    cur->status='I'; 
                   
                    double inf_status = random_real(); //used to determine number of worms
                    double mature_period = normal(mature_period_mean,mature_period_mean_std); //how long worms will live for (all the same)
                    int worm_count = number_worms(cum_sum_prob, inf_status); //the number of worms
                    cur->worm_strength = worm_count; //the strength of worms (equal to number of femeales)
                    inf_indiv.insert(pair<int, agent*>(cur->aid, cur)); //storing the person as infected!
                  
                    for(int i = 1; i <= worm_count; ++ i){ //assigning worms
                        double period = random_real()*mature_period;

                        cur->wvec.push_back(new worm('M', 0, period ,'M'));
                        cur->wvec.push_back(new worm('M', 0, period ,'F'));
                    }
                } 
                else{ //antigen postive but does not have mated worms!
                     // has one mature worm. can be of either gender
                    cur->status='U';
                    uninf_indiv.insert(pair<int, agent*>(cur->aid, cur));
                    double mature_period = normal(mature_period_mean,mature_period_mean_std);
                    if(random_real() < proportion_male_agent){ //has one male mature
                        cur->wvec.push_back(new worm('M', 0, random_real()*mature_period ,'M'));
                    }
                    else{ //has one female mature
                        cur->wvec.push_back(new worm('M', 0, random_real()*mature_period ,'F'));
                    }

                    if(random_real() <= immature_and_unif){ //agent also have immature worm!
                        double immature_period = normal(immature_period_mean,immature_period_mean_std);
                        if(random_real() < proportion_male_agent){ //has one male mature
                            cur->wvec.push_back(new worm('M', random_real()*immature_period, mature_period ,'M'));
                        }
                        else{ //has one female mature
                            cur->wvec.push_back(new worm('M', random_real()*immature_period, mature_period ,'F'));
                        }
                    }
                }
            }
            else {
               if(random_real() <= group_prev*immature_to_antigen){//agent only has immature worm!
                    cur->status='E';
                    pre_indiv.insert(pair<int, agent*>(cur->aid, cur));
                    double immature_period = normal(immature_period_mean,immature_period_mean_std);
                    double mature_period = normal(mature_period_mean,mature_period_mean_std);
                    
                    if (random_real() > twotoone){ //one worm
                        if(random_real() < proportion_male_worm){
                            cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period ,'M'));
                        }
                        else{
                            cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period ,'F'));
                            
                        }
                    }
                    else{ //Two worms
                        for (int i = 0; i < 2; i++) {
                            if (random_real() < proportion_male_worm){ // worm is male!
                                cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period, 'M'));
                            }
                            else{ // worm is female!
                                cur->wvec.push_back(new worm('P', random_real()*immature_period, mature_period, 'F'));
                            }
                        }
                    }
                }
            }
        }
    }

    init_prev = 100 * ant_pos / rpop;
    cout << "Init prev: " << init_prev << "%" << endl;
}

vector<double> region::prob_worms(double prev, int n_worms){
    
    
    double no_worm_prev = 1.0 - matedtonot * prev; //proportion of pop without any worms
    double worm_mean = agg_param * (pow((1/no_worm_prev), (1/agg_param)) - 1.0); //mean worm burden in group (assuming negative binomial)
    double worm_prob = 0.0;
    
    
    vector<double> cum_sum_prob;

    for (int i = 1; i <= n_worms; ++i){ //now working out prob of each worm burden assuming negative binomial
        worm_prob +=  (tgamma(i + agg_param) / (tgamma(agg_param) * factorial(i))) * pow((worm_mean / (agg_param + worm_mean)),i) * pow((agg_param / (agg_param + worm_mean)),agg_param);  
        cum_sum_prob.push_back(worm_prob); // storing
    } 
    //Now calculating relative prob of each non-zero worm burden

    for (auto& n : cum_sum_prob) n = n/cum_sum_prob.back(); //relative prob of worms
    return cum_sum_prob;

}

int region::factorial(int n) {
    if (n == 0)
       return 1;
    return n * factorial(n - 1);
}

int region::number_worms(vector<double> cum_sum_prob, double prob){
    int n_worms = cum_sum_prob.size();

    for ( int i = 1; i <= n_worms; ++i){
        if (prob <= cum_sum_prob[i]) return i;
    }
    
    cout << "ERROR" << endl;

    return n_worms;
}