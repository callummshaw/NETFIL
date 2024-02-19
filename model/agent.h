#ifndef agent_hpp
#define agent_hpp

#include "params.h"

#include<iostream>

using namespace std;

//Classes that we need
class group;
class drugs;
class mda_strat;

class agent; //people in the model
class worm; //worms!

class worm{
public:
    char status; // stage the worm is in
    int age_immature; // immature age
    int age_mature; // mature age
    int age_mda; // time since mda
    double mda_sterile; // to track sterilisation from MDA
    char sex; // sex of the worm 

    worm(char s, int ia, int ma, char sx){
        
        status = s;
        age_immature = ia;
        age_mature = ma;
        sex = sx;

        age_mda = 0;
        mda_sterile = 1.0;
    }

    void update(int dt);
};

class agent{
public:
    
    int aid; // agent's id
    int age; // agent's age
   
    char status; // epi status
    
    double lastwormtime; // time since last adult worm
    double worm_strength;//keep track of number and sterility of mature female worms when there is an adult male!

    bool ChangedEpiToday;

    group *dgp; //daytime group
    group *ngp; //nightime group
    
    vector<worm*> wvec;
   
    agent(int aid, int age = -1);

    ~agent();

    void sim_bites(double prev, char time, double c, double biteload, double w2n);
    void update(int day, int year, int dt);
    void mda(drugs drug);

};

#endif /* agent_hpp */