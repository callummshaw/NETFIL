#include <iostream>
#include <ctime>
#include <unistd.h>

#include "main.h"
#include "mda.h"

using namespace std;

 
int SimulationNumber;

string prv_out_loc;

int main(int argc, const char * argv[]){
    
    prv_out_loc = argv[1];
   
    region *rgn = new region(region_id, region_name);
    
    string mda_data = datadir; mda_data = mda_data + MDA_params;

    //Counting the number of different simulations we will perform
    int MDAScenario_count = count_mda_scenarios(mda_data);
   // cout << "There are " << MDAScenario_count << " scenarios" << endl;

    //now looping over scenarios
    for (int scenario_count = 0; scenario_count < MDAScenario_count; ++scenario_count){

        //generating mda strategy!
        mda_strat strategy = get_mda_strat(mda_data, scenario_count + 1);

        //Now looping over simulations
         for (int i = 0; i < strategy.NumSims; ++i){
            
            //resetting the populations from previous simulation
            rgn->reset_population();
           
            //run run the simulation year by year
            for(int year = 0; year < sim_years; ++year){
                
                rgn->sim(year, strategy);
              
            }
        }
    }

    return 0;
    
}