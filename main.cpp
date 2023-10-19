#include <iostream>
#include <ctime>
#include <unistd.h>
#include <stdlib.h>
#include <random>


using namespace std;

unsigned seed;

random_device r;
seed_seq seeds{r(), r(), r(), r(), r(), r(), r(), r()};
mt19937 gen(seeds);