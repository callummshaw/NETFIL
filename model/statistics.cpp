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