#include "agent.h"
#include "mda.h"
#include "network.h"
#include <limits>

//Constructer of agent
agent::agent(int aid, int age = -1){
    this->aid = aid;
    this->age = age;
    
   
    status = 's';

    worm_strength = 0;

    lastwormweek = - std::numeric_limits<double>::infinity();
}

agent::~agent(){

    g_p = NULL;

    for(int i = 0; i < wvec.size(); ++i){
        delete wvec[i];
    }
    wvec.clear();

}

void agent::sim_bites(double prev, char time, double c, double theta){
    
    double pos_inf_bite_rate = c*prev*theta;
    
    if (time == 'D') pos_inf_bite_rate *= work_bite_rate; //bites during working hours
    else pos_inf_bite_rate *= offwork_bite_rate;  //bites outside working hours

    int InfectiveBites = poisson(pos_inf_bite_rate); // the number of bites agent will recieve

    for(int i = 0; i < InfectiveBites; ++i){ //looping through bites
        int immature_period = normal(immature_period_mean, immature_period_mean_std); //immature period of worm
        int mature_period = normal(mature_period_mean, mature_period_mean_std); //mature period of worm

        if (random_real() < proportion_male){ // worm is male!
            wvec.push_back(new worm('I', immature_period, mature_period, 'M'));
        }
        else{ // worm is female!
            wvec.push_back(new worm('I', immature_period, mature_period, 'F'));
        }
    }
}

void agent::mda(drugs drug){
    if(wvec.size() > 0){ //if person has worms
        double rr = random_real(); //same thing will occur to all worms!
        for(int i = 0; i < wvec.size(); i++){ // looping through worms
            wvec[i]->age_mda = drug.SterDur;
            if (rr <= drug.KillProb){
                wvec[i]->status = 'D';
            }
            else if(rr <= drug.KillProb + drug.FullSterProb){ // sterilise worms with probability FullSterProb
                wvec[i]->mda_sterile = 0.0; //worm is sterile!
            }
            else if (rr <= drug.KillProb + drug.FullSterProb + drug.PartSterProb){ //Partially sterilise with probability PartSterProb
                wvec[i]->mda_sterile = min(wvec[i]->mda_sterile, 1 - drug.PartSterMagnitude); // worm is partially sterile, we also ensure that mda does NOT increase infectivity of an already sterilised worm
            }
        }
    }
}

//update people!
void agent::update(int week){

    //Firstly update status of all worms in the body
    if(wvec.size() > 0){ //Now will update each worm
        for(int i = 0; i < wvec.size();){ 
            wvec[i]->update();
            
            //Removing dead worms!
            if(wvec[i]->status == 'd'){ 
                delete wvec[i];
                wvec.erase(wvec.begin() + i);
            }
            else ++i; // iterator increment is placed here because if a worm is removed from the list the next worm will now have the index of deleted worm
        }
    }

    char prevstatus = status; // agents previous infection status

    bool mature_worm = false; //has mature worms of either sex (fertile does not matter)

    double worm_strength_female = 0.0; //counter to track strength of females
    double worm_strength_male = 0.0; // counter to track strength of males

    //update agent's epi status
    if(wvec.size() == 0){ //no Worms!
        status = 'S'; // therfore in S class!
        worm_strength = 0.0; //agent has no worm strength
    } 
    else{ //person has worms!
        for(int i =0; i < wvec.size(); ++i){ //iterating over worms
            if(wvec[i]->status == 'M'){ // if worm is mature
                mature_worm = true;
                if (wvec[i]->sex == 'M'){//male worm
                    if (wvec[i]->mda_sterile > 0){
                        worm_strength_male += wvec[i]->mda_sterile;
                    }
                }
                else{ // female worm
                    if (wvec[i]->mda_sterile > 0){
                        worm_strength_female += wvec[i]->mda_sterile;
                    }
                } 
            } 

        }
        if((worm_strength_female > 0) && (worm_strength_male > 0)){ //agent is infectious!
            status = 'I';// person is infectious

            if(worm_strength_male > 1.0) worm_strength_male = 1.0; //We assume polygamous worms that relies on females

            worm_strength = worm_strength_male * worm_strength_female; // total worm strength

        } 
        else{ 
            worm_strength = 0.0;
            if(mature_worm) status = 'U';//agent has worm(s) but not a mature fertile set
            else{
                status = 'E'; //agent has only immature worm(s)
            }
        }
    }

    //Record if worm has died!
    if((prevstatus == 'U' || prevstatus == 'I') && (status == 'S' || status == 'E')){ // all mature worms have died!
        lastwormweek = week;
    }
    
}

//update worms!
void worm::update(){
    //Three statuses: I-immature, M-mature, D-Dead 

    if(status == 'I'){
        if(age_immature == 0) status = 'M'; //worm is now mature!
        else --age_immature; // countdown
    }

    if(status == 'M'){
        if(age_mature == 0) status = 'D'; //worm is dead!
        else --age_mature;
    }

    if (age_mda > 0) --age_mda; //now looking at sterilisation status of worm!
    else mda_sterile = 1.0; //sterile period now over!
}