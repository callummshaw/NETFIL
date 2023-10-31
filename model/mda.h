#ifndef mda_h
#include "params.h"
#define mda_h
#include <vector>
#include <algorithm>

using namespace std;

class mda_strat;
class drugs;


int count_mda_scenarios(string file);
mda_strat get_mda_strat(string filename, int N);


class drugs{ //Drug profile for MDA
public:
    double KillProb; // prob will kill worms
    double FullSterProb; // prob will fully sterilise
    double PartSterProb; // prob will partially sterilise
    double SterDur; // duration of sterilisation
    double PartSterMagnitude; //magnitude of partial sterilisation

    
    drugs(double K=0.0, double FSP=0.0, double PSP=0.0, double SD=0.0, double PSM=0.0){
        KillProb = K;
        FullSterProb = FSP;
        PartSterProb = PSP;
        SterDur = SD;
        PartSterMagnitude = PSM;
    }
    
    
    void print_drugs(){ //print drug profiles!
        cout << "KillProb: " << KillProb << endl;
        cout << "FullSterProb: " << FullSterProb << endl;
        cout << "PartSterProb: " << PartSterProb << endl;
        cout << "SterDur: " << SterDur << endl;
        cout << "PartSterMagnitude: " << PartSterMagnitude << endl;
    }
    
};

class mda_strat{
public:
    double Coverage; //mda coverage
    drugs drug; //drug profile
    int min_age; //min age to take MDA
    int StartYear; //mda starty year
    int NumRounds; //number of rounds
    int YearsBetweenRounds; //time between roundss
    vector<int> MDAYears; //vector storing all mda years
    int NumSims; //number of simulations

    mda_strat(double C, drugs D, int MN, int S, int N, int Y, int NS){
        Coverage = C;
        drug = D;
        min_age = MN;
        StartYear = S;
        NumRounds = N;
        YearsBetweenRounds = Y;
        NumSims = NS;

        MDAYears.resize(NumRounds);
        for(int i = 0; i<NumRounds; ++i){
            MDAYears[i] = StartYear + i * YearsBetweenRounds;
            cout << MDAYears[i] << " is a MDA year" << endl;
        }   
    }

    void print_mda_strat(){
        cout << endl;
        cout << "Coverage: " << Coverage * 100 << "%" << endl;
        drug.print_drugs();
        cout << "Start Year: " << StartYear << endl;
        cout << "Number of Rounds: " << NumRounds << endl;
        cout << "Number of Simulations: " << NumSims << endl;
        cout << endl;
    }
    
    bool is_mda_year(int Year){ // Returns true if given year is an MDA year for the strategy
        return find(MDAYears.begin(), MDAYears.end(), Year) != MDAYears.end();
    }

};

#endif /* mda_h */