#include "network.h"
#include "rng.h"
#include <cstring>

extern string prv_out_loc;
extern int SimulationNumber;

void region::output_epidemics(int year, int day, mda_strat strategy){
    
    //total pop
    double pop_total = 0;
    double inf_total = 0;
    double ant_total = 0;

    //worm burdens
    double immature_worm_only = 0;
    double non_mated_adult = 0;
    double one_mated_adult = 0;
    double two_mated_adult = 0;
    double three_mated_adult = 0;
    double four_mated_adult = 0;
    double five_mated_adult = 0;
    double six_mated_adult = 0;
    double seven_mated_adult = 0;
    double eight_mated_adult = 0;
    double nine_mated_adult = 0;
    double tenplus_mated_adult = 0;
    
    vector<double> inf_groups;
    inf_groups.resize(groups.size());

    vector<double> antigen_pos_groups;
    antigen_pos_groups.resize(groups.size());

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
        group *grp = j->second;
        pop_total += grp->group_pop.size();

        //now over people
        for(map<int,agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            agent *a = k-> second;
            int age = int(a->age/365);

            if(a->status == 'I'){//person is infectious
                ++inf_groups[j->first - 1];
                double ws = a->worm_strength;
                ++inf_total;
                if (ws <= 1) ++one_mated_adult;
                if (ws > 1 && ws <= 2) ++two_mated_adult;
                if (ws > 2 && ws <= 3) ++three_mated_adult;
                if (ws > 3 && ws <= 4) ++four_mated_adult;
                if (ws > 4 && ws <= 5) ++five_mated_adult;
                if (ws > 5 && ws <= 6) ++six_mated_adult;
                if (ws > 6 && ws <= 7) ++seven_mated_adult;
                if (ws > 7 && ws <= 8) ++eight_mated_adult;
                if (ws > 8 && ws <= 9) ++nine_mated_adult;
                if (ws > 9 ) ++tenplus_mated_adult;

            }
            if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, (year*365 +day) - a->lastwormtime) ){ //all people infected with any number of mature worms or who still have lingering antibodies are counted
                
                ++antigen_pos_groups[j->first - 1];
                ++ant_total;
            }
            if (a->status == 'U') ++non_mated_adult;
            if (a->status == 'E') ++immature_worm_only;
        }

    }
    if (day == 0){
    cout << endl;
    
    cout << year+start_year << ": " << "prepatent = " << pre_indiv.size() << " uninfectious = " << uninf_indiv.size() << " infectious = " << inf_indiv.size() << " antigen positive = " << ant_total << endl;
    cout << "overall mf prevalence = " << fixed << setprecision(2) << inf_indiv.size()/(double)rpop*100 << "%" << endl;
    cout<< "overall ant prevalence = " << fixed << setprecision(2) << ant_total/(double)rpop*100 << "%" << endl;
    cout<< "overall ratio prevalence = " << fixed << setprecision(2) << ant_total/inf_total << endl;
    }
    string prv_dat = outdir;    prv_dat = prv_dat + prv_out_loc; 
    ofstream out;   ifstream in;
    in.open(prv_dat.c_str()); // try opening the target for output
    if(!in){ // if it doesn't exist write a heading
        out.open(prv_dat.c_str());
        out << "SimulationNumber,";
        out << "Year,";
        out << "Day,";
        out << "Agg,";
        out << "Theta1,";
        out << "Theta2,";
        out << "WorktoNot,";
        out << "AntandImmature,";
        out << "OnlyImmature,"; 
        out << "MDACoverageAttempted,";
        out << "MDAKillProb,";
        out << "MDAFullSterProb,";
        out << "MDAPartSterProb,";
        out << "MDASterDur,";
        out << "MDAPartSterMagnitude";
        out << "MDAMinAge,";
        out << "MDAStartYear,";
        out << "MDANumRounds,";
        out << "MDAYearsBetweenRounds,";
        out << "AchievedMDACoverage,";
        out << "SimYears,";
        out << "Pop_total,";
        out << "inf_total,"; 
        out << "ant_total,";
        out << "treated,";
        out << "immature_worm,";
        out << "non_mated,";
        out << "one_mated,";
        out << "two_mated,";
        out << "three_mated,";
        out << "four_mated,";
        out << "five_mated,";
        out << "six_mated,";
        out << "seven_mated,";
        out << "eight_mated,";
        out << "nine_mated,";
        out << "tenplus_mated,";
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << "pop_" << group_numbers[j -> second -> gid] << ","; 
        }
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << "mf_" << group_numbers[j -> second -> gid] << ","; 
        }
        out << endl;
        out.close();
    }
    else in.close();

    //write the prevalence for whole populations, by gender, by age group and for each village
    out.open(prv_dat.c_str(), ios::app);
    
    out << SimulationNumber << ",";
    out << year + start_year << ",";
    out << day << ",";
    out << agg_param << ",";
    out << theta1 << ",";
    out << theta2 << ",";
    out << worktonot << ",";
    out << immature_and_ant  << ",";
    out << immature_to_antigen << ",";
    out << strategy.Coverage << ",";
    out << strategy.drug.KillProb << ",";
    out << strategy.drug.FullSterProb << ",";
    out << strategy.drug.PartSterProb << ",";
    out << strategy.drug.SterDur << ",";
    out << strategy.drug.PartSterMagnitude << ",";
    out << strategy.StartYear << ",";
    out << strategy.NumRounds << ",";
    out << strategy.YearsBetweenRounds << ",";
    out << achieved_coverage[year] << ",";
    out << sim_years << ",";
    out << pop_total << ",";
    out << inf_total << ",";
    out << ant_total << ",";
    out << number_treated[year] << ",";
    out << immature_worm_only << ",";
    out << non_mated_adult << ",";
    out << one_mated_adult << ",";
    out << two_mated_adult<< ",";
    out << three_mated_adult << ",";
    out << four_mated_adult << ",";
    out << five_mated_adult<< ",";
    out << six_mated_adult << ",";
    out << seven_mated_adult<< ",";
    out << eight_mated_adult << ",";
    out << nine_mated_adult << ",";
    out << tenplus_mated_adult<< ",";
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        double n_village = (j -> second -> group_pop).size();
        if(n_village==0) out << "NA,"; // there's a chance that populations in small villages might drop to zero - this is to avoid crashes in that situation
        else out << n_village << ",";
    }
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        double n_village = (j -> second -> group_pop).size();
        if(n_village==0) out << "NA,"; // there's a chance that populations in small villages might drop to zero - this is to avoid crashes in that situation
        else out <<  inf_groups[j -> first - 1] << ",";
    }
    out << endl;
    out.close();
}


void region::output_abc_epidemics(int year){

    int icc_2016 = 0;
    int all_2016 = 1000; //number of people we need for ICC calc

    vector<int> group_number;
    vector<float> agev;
    vector<int> inf;
    vector<unsigned> keys; // Vector that is used to choose the groups that we want (keys for the map of groups) 
    int ppg = 50; //max people per group
    int tot_groups = groups.size();

   
    if (year + start_year == 2014) {
        
        int ant_count = 0;
        
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
            group *grp = j->second;
            //now over people
            for(map<int,agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                agent *a = k-> second;
                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime) ){ 
                    ++ant_count; 
                }
            }

        }
        mf_to_ant_2014 =  (double)inf_indiv.size()/(double)ant_count;
    }

    if (year + start_year == 2016) {
        double mf_2016;
        double ant_2016;

        int ant_count = 0;
        
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
            group *grp = j->second;
            //now over people
            for(map<int,agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                agent *a = k-> second;
                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime) ){ 
                    ++ant_count; 
                }
            }

        }

        int ppg = 100; //max people per group
        int tot_groups = groups.size();
        vector<unsigned> keys; // Vector that is used to choose the groups that we want (keys for the map of groups)

        for (unsigned int i = 0; i < tot_groups; ++i){
            keys.push_back(i + 1);
        }

        shuffle(keys.begin(),keys.end(),gen);

        for (auto const &i: keys) {//iterating through groups
            int npeople = 0;
            for (auto  const& k: groups[i]->group_pop){//iterating through people selected group
                if( drand48()>0.5){ //inducing some randomness into the process, choosing a person
                    
                    agent *a = k.second;
                    int age = int(a->age/365);

                    icc_2016 += 1;
                    npeople += 1;

                    if (icc_2016 < all_2016){
                        group_number.push_back(i);
                        agev.push_back(age); 
                        if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime))
                        {
                            inf.push_back(1);
                        } else{
                            inf.push_back(0);
                        }
                    }
        
                    
                }

                if (npeople == ppg) goto next_group2016; 
            }
            next_group2016:;
        }
        
        mf_2016 =  (double)inf_indiv.size()/(double)rpop*100;
        ant_2016 = (double)ant_count/(double)rpop*100;
        int number_rows = inf.size();
        
        string survey_out = outdir;    survey_out = survey_out + "survey_" + prv_out_loc;
        string icc_out = outdir;   icc_out = icc_out + "icc_" + prv_out_loc; 

        ofstream out;   ifstream in;
        out.open(survey_out.c_str());
        out <<  "Ratio_2014 ";
        out <<  "Antigen_2016 ";
        out <<  "MF_2016"; 
        out << endl;
        out << mf_to_ant_2014;
        out << ' ' << ant_2016;
        out << ' ' << mf_2016;
        out << endl;
        out.close();

        out.open(icc_out.c_str());
        out << "Village,";
        out << "Status,";
        out << "Age,";
        out << endl;
        for (int i = 0; i < number_rows; i++){
            out << group_number[i]  << "," << inf[i] << "," << agev[i]  << endl;
        }
        out << endl;
        out.close();
    }
}

void region::output_abc_epidemics_single(int year){

    if (year + start_year == 2014) {
        
        int ant_count = 0;
        
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
            group *grp = j->second;
            //now over people
            for(map<int,agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                agent *a = k-> second;
                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime) ){ 
                    ++ant_count; 
                }
            }

        }
        mf_to_ant_2014 =  (double)inf_indiv.size()/(double)ant_count;
    }

    if (year + start_year == 2016) {
        double mf_2016;
        double ant_2016;

        int ant_count = 0;
        
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
            group *grp = j->second;
            //now over people
            for(map<int,agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                agent *a = k-> second;
                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime) ){ 
                    ++ant_count; 
                }
            }

        }

        mf_2016 =  (double)inf_indiv.size()/(double)rpop*100;
        ant_2016 = (double)ant_count/(double)rpop*100;
        
        string init_out = "summary_stats_temp.txt";

        ofstream out;   ifstream in;
        out.open(init_out.c_str());
        out <<  "Ratio_2014 ";
        out <<  "Antigen_2016 ";
        out <<  "MF_2016"; 
        out << endl;
        out << mf_to_ant_2014;
        out << ' ' << ant_2016;
        out << ' ' << mf_2016;
        out << endl;
        out.close();
    }
}