#include "network.h"
#include <cstring>


extern string prv_out_loc;
extern int SimulationNumber;

void region::output_epidemics(int year, mda_strat strategy){
    
    // community based
    double pop_0_9 = 0;
    double inf_0_9 = 0;
    double ant_0_9 = 0;
    
    double pop_10_19 = 0;
    double inf_10_19 = 0;
    double ant_10_19 = 0;
    
    double pop_20_29 = 0;
    double inf_20_29 = 0;
    double ant_20_29 = 0;

    double pop_30_39 = 0;
    double inf_30_39 = 0;
    double ant_30_39 = 0;

    double pop_40_49 = 0;
    double inf_40_49 = 0;
    double ant_40_49 = 0;

    double pop_50_59 = 0;
    double inf_50_59 = 0;
    double ant_50_59 = 0;

    double pop_60_69 = 0;
    double inf_60_69 = 0;
    double ant_60_69 = 0;

    double pop_70_79 = 0;
    double inf_70_79 = 0;
    double ant_70_79 = 0;

    double pop_80_plus = 0;
    double inf_80_plus = 0;
    double ant_80_plus = 0;

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


            // add to appropriate tally of total population by age
            if(age <= 9)++pop_0_9;
            if(age >=10 && age <=19)++pop_10_19;
            if(age >=20 && age <=29)++pop_20_29;
            if(age >=30 && age <=39)++pop_30_39;
            if(age >=40 && age <=49)++pop_40_49;
            if(age >=50 && age <=59)++pop_50_59;
            if(age >=60 && age <=69)++pop_60_69;
            if(age >=70 && age <=79)++pop_70_79;
            if(age >=80)++pop_80_plus;

            if(a->status == 'I'){//person is infectious

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

                if(age <= 9)++inf_0_9;
                if(age >=10 && age <=19)++inf_10_19;
                if(age >=20 && age <=29)++inf_20_29;
                if(age >=30 && age <=39)++inf_30_39;
                if(age >=40 && age <=49)++inf_40_49;
                if(age >=50 && age <=59)++inf_50_59;
                if(age >=60 && age <=69)++inf_60_69;
                if(age >=70 && age <=79)++inf_70_79;
                if(age >=80)++inf_80_plus;
            }
            if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime) ){ //all people infected with any number of mature worms or who still have lingering antibodies are counted
                
                ++antigen_pos_groups[j->first - 1];
                ++ant_total;
                
                if(age <= 9)++ant_0_9;
                if(age >=10 && age <=19)++ant_10_19;
                if(age >=20 && age <=29)++ant_20_29;
                if(age >=30 && age <=39)++ant_30_39;
                if(age >=40 && age <=49)++ant_40_49;
                if(age >=50 && age <=59)++ant_50_59;
                if(age >=60 && age <=69)++ant_60_69;
                if(age >=70 && age <=79)++ant_70_79;
                if(age >=80)++ant_80_plus;
            }
            if (a->status == 'U') ++non_mated_adult;
            if (a->status == 'E') ++immature_worm_only;
        }

    }

    cout << endl;
    cout << year+start_year << ": " << "prepatent = " << pre_indiv.size() << " uninfectious = " << uninf_indiv.size() << " infectious = " << inf_indiv.size() << " antigen positive = " << ant_total << endl;
    cout << "overall mf prevalence = " << fixed << setprecision(2) << inf_indiv.size()/(double)rpop*100 << "%" << endl;
    cout<< "overall ant prevalence = " << fixed << setprecision(2) << ant_total/(double)rpop*100 << "%" << endl;
    cout<< "overall ratio prevalence = " << fixed << setprecision(2) << ant_total/inf_total << endl;

    string prv_dat = outdir;    prv_dat = prv_dat + prv_out_loc; 

    ofstream out;   ifstream in;
    in.open(prv_dat.c_str()); // try opening the target for output
    if(!in){ // if it doesn't exist write a heading
        out.open(prv_dat.c_str());
        out << "SimulationNumber,";
        out << "Year,";
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
        out << "Pop_0_9,";
        out << "Pop_10_19,";
        out << "Pop_20_29,";
        out << "Pop_30_39,";
        out << "Pop_40_49,";
        out << "Pop_50_59,";
        out << "Pop_60_69,";
        out << "Pop_70_79,";
        out << "Pop_80_plus,";
        out << "inf_total,"; 
        out << "inf_0_9,";
        out << "inf_10_19,";
        out << "inf_20_29,";
        out << "inf_30_39,";
        out << "inf_40_49,";
        out << "inf_50_59,";
        out << "inf_60_69,";
        out << "inf_70_79,";
        out << "inf_80_plus,";
        out << "ant_total,";
        out << "ant_0_9,";
        out << "ant_10_19,";
        out << "ant_20_29,";
        out << "ant_30_39,";
        out << "ant_40_49,";
        out << "ant_50_59,";
        out << "ant_60_69,";
        out << "ant_70_79,";
        out << "ant_80_plus,";
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
            out << "Prevalence" << group_numbers[j -> second -> gid] << ","; //For each village prints "Prevalence<Village Name>,"
        }
        for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << "Prev_Antigen" << group_numbers[j -> second -> gid] << ","; //For each village prints "Prev_Antigen<Village Name>,"
        }
        out << endl;
        out.close();
    }
    else in.close();

    //write the prevalence for whole populations, by gender, by age group and for each village
    out.open(prv_dat.c_str(), ios::app);
    
    out << SimulationNumber << ",";
    out << year + start_year << ",";
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
    out << pop_0_9<< ",";
    out << pop_10_19 << ",";
    out << pop_20_29 << ",";
    out << pop_30_39 << ",";
    out << pop_40_49 << ",";
    out << pop_50_59 << ",";
    out << pop_60_69 << ",";
    out << pop_70_79 << ",";
    out << pop_80_plus << ",";
    out << inf_total << ",";
    out << inf_0_9<< ",";
    out << inf_10_19 << ",";
    out << inf_20_29 << ",";
    out << inf_30_39 << ",";
    out << inf_40_49 << ",";
    out << inf_50_59 << ",";
    out << inf_60_69 << ",";
    out << inf_70_79 << ",";
    out << inf_80_plus << ",";
    out << ant_total << ",";
    out << ant_0_9<< ",";
    out << ant_10_19 << ",";
    out << ant_20_29 << ",";
    out << ant_30_39 << ",";
    out << ant_40_49 << ",";
    out << ant_50_59 << ",";
    out << ant_60_69 << ",";
    out << ant_70_79 << ",";
    out << ant_80_plus << ",";
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
        else out <<  inf_groups[j -> first - 1]/(double)n_village << ",";
    }
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        double n_village = (j -> second -> group_pop).size();
        if(n_village==0) out << "NA,"; // there's a chance that populations in small villages might drop to zero - this it to avoid crashes in that situation
        else out <<  antigen_pos_groups[j -> first - 1]/(double)n_village << ",";
    }

    out << endl;
    out.close();
}
void region::output_abc_epidemics(int year){

    vector<int> group_number;
    vector<float> agev;
    vector<int> inf;
    
    //Survey Data
    //2014
    int under_2014_real = 425;
    int over_2014_real = 653;
    int mf_2014_real = 20;

    //2016
    int under_2016_real = 879;
    int over_2016_real = 1617;
    int mf_2016_real = 86;


    // Model Data
    //2016
    int under_2016 = 0;
    int over_2016 = 0;

    int over_2016_ant = 0;
    int under_2016_ant = 0;

    int mf_2016 = 0;
    int mf_2016_positive = 0;

    int icc_2016 = 0;
    //2014
    int under_2014 = 0;
    int over_2014 = 0;

    int over_2014_ant = 0;
    int under_2014_ant = 0;

    int mf_2014 = 0;
    int mf_2014_positive = 0;

    int ppg = 50; //max people per group
    int tot_groups = groups.size();

    vector<unsigned> keys; // Vector that is used to choose the groups that we want (keys for the map of groups)

    for (unsigned int i = 0; i < tot_groups; ++i) //Filling vector with numbers from 1-number of groups
    {
        keys.push_back(i + 1);
    }
    random_shuffle(keys.begin(),keys.end());

    if (year + start_year == 2014) {
        while (under_2014 < under_2014_real || over_2014 < over_2014_real) {
            for (auto const &i: keys) {//iterating through groups
                int npeople = 0;
                for (auto  const& k: groups[i]->group_pop){//iterating through people selected group
                    if( drand48()>0.5) { //inducing some randomness into the process, choosing a person
                        npeople += 1;
                        
                        agent *a = k.second;
                        int age = int(a->age/365);

                        if( age < 25){
                            if (under_2014 < under_2014_real) {
                                under_2014 += 1; //tested
                                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime) ){ //ant pos
                                    under_2014_ant += 1;
                                    if (mf_2014 < mf_2014_real) {
                                        mf_2014 += 1;
                                        if(a->status == 'I'){
                                            mf_2014_positive += 1;
                                        }
                                    }
                                }
                            }
                        }
                        else{
                            if (over_2014 < over_2014_real) {
                                over_2014 += 1; //test
                                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime) ){ //ant pos
                                    over_2014_ant += 1;
                                    if (mf_2014 < mf_2014_real) {
                                        mf_2014 += 1;
                                        if(a->status == 'I'){
                                            mf_2014_positive += 1;
                                        }
                                    }
                                }
                            }
                        }
                        
                    }
                    if (npeople == ppg) goto next_group; 
                }
                next_group:;
            }
        }

        data_ant_u = under_2014_ant;
        data_ant_o = over_2014_ant;
        data_mf = mf_2014_positive;
        
    }

    if (year + start_year == 2016) {
        while (under_2016 < under_2016_real || over_2016 < over_2016_real) {
            for (auto const &i: keys) {//iterating through groups
                int npeople = 0;
                for (auto  const& k: groups[i]->group_pop){//iterating through people selected group
                    if( drand48()>0.5) { //inducing some randomness into the process, choosing a person
                        npeople += 1;

                        agent *a = k.second;
                        int age = int(a->age/365);

                        icc_2016 += 1;
                        if (icc_2016 < (under_2016_real+over_2016_real)){
                            group_number.push_back(i);
                            agev.push_back(age); 
                            if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime))
                            {
                                inf.push_back(1);
                            } else{
                                inf.push_back(0);
                            }
                        }
                        
                        if( age < 20){
                            if (under_2016 < under_2016_real) {
                                under_2016 += 1; //tested
                                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime)){ //ant pos
                                    under_2016_ant += 1;
                                    if (mf_2016 < mf_2016_real) {
                                        mf_2016 += 1;
                                        if(a->status == 'I'){
                                            mf_2016_positive += 1;
                                        }
                                    }

                                }
                            }
                        }
                        else{
                            if (over_2016 < over_2016_real) {
                                over_2016 += 1; //test
                                if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DailyProbLoseAntigen, year*365 - a->lastwormtime)){ //ant pos
                                    over_2016_ant += 1;
                                    if (mf_2016 < mf_2016_real) {
                                        mf_2016 += 1;
                                        if(a->status == 'I'){
                                            mf_2016_positive += 1;
                                        }
                                    }
                                }
                            }
                        }
                        
                    }

                    if (npeople == ppg) goto next_group2016; 
                }
                next_group2016:;
            }
        }
    }

    string survey_out = outdir;    survey_out = survey_out + "survey_" + prv_out_loc;
    string icc_out = outdir;   icc_out = icc_out + "icc_" + prv_out_loc; 

    int number_rows = inf.size();

    if (year+start_year == 2016){
        ofstream out;   ifstream in;
        out.open(survey_out.c_str());
        out <<  "Antigen_2014_Under,";
        out <<  "Antigen_2014_Over,";
        out <<  "MF_2014_Total,";
        out <<  "Antigen_2016_Under,";
        out <<  "Antigen_2016_Over,";
        out <<  "MF_2016_Total";
        out << endl;
        out << data_ant_u;
        out << ',' << data_ant_o;
        out << ',' << data_mf;
        out << ',' << under_2016_ant;
        out << ',' << over_2016_ant;
        out << ',' << mf_2016_positive;
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

void region::output_abc_epidemics(int year){

    int icc_2016 = 0;
    int all_2016 = 2500; //number of people we need for ICC calc

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

        random_shuffle(keys.begin(),keys.end());

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