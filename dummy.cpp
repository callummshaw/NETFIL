#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <set>
#include <map>
#include <stdio.h>
#include <string.h>
#include <cstring>

using namespace std;

int main(){
    ifstream in;
    string line, file;

    double drummystore[16]; 
    //reading in age distribution
    file = "data/";    file = file + "pop_age_dist.csv";
    in.open(file.c_str());
    
    //skip the description
    while(getline(in, line)){
        if(line[0] == '*') continue;
        if(line.length() <= 1) continue;  //empty line with carriage return
        break;
    }
    int ii = 0;

    while(getline(in, line)){
        drummystore[ii] = atof(line.c_str());
        ii++;
    }
    cout<<"HERE"<<endl;
    for (int i = 0; i < 16; i++) cout << drummystore[i] <<endl;
return 0;
}
