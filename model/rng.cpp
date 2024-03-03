#include "rng.h"
#include <chrono>
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 gen(seed); // Seed the generator

// Warm up the generator
void warmup() {
    
    for(int i = 0; i < 10000; ++i) {
        gen();
    }
}

int dummy = (warmup(), 0);