#include "params.h"

#include<iostream>

using namespace std;

//Classes that we need
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

    void update();
};

class agent{
public:
    
    int aid; // agent's id
    int age; // agent's age
   
    char status; // epi status
    
    double lastwormweek; // time since last adult worm
    double worm_strength;//keep track of number and sterility of mature female worms when there is an adult male!

    vector<worm*> wvec;
    group *g_p;

    agent(int aid, int age = -1,  group *g_p = NULL);

    ~agent();

    void sim_bites(double prev, char time, double c, double theta);
    void update(int week);
    void mda(drugs drug);

};