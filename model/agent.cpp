#include "agent.h"
#include "mda.h"
#include "network.h"
#include <limits>

//Constructer of agent
agent::agent(int aid,  double bite_shape, int age){
    this->aid = aid;
    this->age = age;
   
    ChangedEpiToday = false;
    status = 'S';

    worm_strength = 0;

    lastwormtime = - std::numeric_limits<double>::infinity();

    bite_scale = bite_gamma(bite_shape, 1/bite_shape);
}

agent::~agent(){
    dgp = NULL;
    ngp = NULL;
    for(int i = 0; i < wvec.size(); ++i){
        delete wvec[i];
    }
    wvec.clear();

}

void agent::sim_bites(double c, double worktonot, bool single){
    
    int total_bites;

    if(single){
        total_bites = poisson(c *  ngp->night_strength * bite_scale);
    }else{
        int day_bites;
        int night_bites;
        

        day_bites  = poisson(c *  dgp->day_strength * bite_scale * worktonot);
        night_bites = poisson(c * ngp->night_strength * bite_scale * (1.0 - worktonot));

        total_bites = day_bites + night_bites;
    }
    
    for(int i = 0; i < total_bites; ++i){ //looping through infective bites and assigning worms
        int immature_period = normal(immature_period_mean, immature_period_mean_std); //immature period of worm
        int mature_period = normal(mature_period_mean, mature_period_mean_std); //mature period of worm

        if (random_real() < proportion_male_worm){ // worm is male!
            wvec.push_back(new worm('P', immature_period, mature_period, 'M'));
        }
        else{ // worm is female!
            wvec.push_back(new worm('P', immature_period, mature_period, 'F'));
        }
    }

    if(total_bites > 0 && status == 'S') status = 'E';
}

void agent::mda(drugs drug){
    if(wvec.size() > 0){ //if person has worms
        double rr = random_real(); //same thing will occur to all worms!
        for(int i = 0; i < wvec.size(); i++){ // looping through worms
            wvec[i]->age_mda = drug.SterDur*365;
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
void agent::update(int year, int day, int dt){
    //Firstly update status of all worms in the body
    if(wvec.size() > 0){ //Now will update each worm
        for(int i = 0; i < wvec.size();){ 
            wvec[i]->update(dt);
            
            //Removing dead worms!
            if(wvec[i]->status == 'D'){ 
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
        status = 'S'; // therefore in S class!
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
                else if (wvec[i]->sex == 'F'){// female worm
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
        lastwormtime = year * 365 + day*dt;
    }
    
}

//update worms!
void worm::update(int dt){
    //Three statuses: P-immature, M-mature, D-Dead 

    if(status == 'P'){
        if(age_immature <= 0){ 
            status = 'M';
            } //worm is now mature!
        else age_immature -= dt; // countdown
    }

    if(status == 'M'){
        if(age_mature <= 0) status = 'D'; //worm is dead!
        else age_mature -= dt;
    }

    if (age_mda > 0) age_mda -= dt; //now looking at sterilisation status of worm!
    else mda_sterile = 1.0; //sterile period now over!
}