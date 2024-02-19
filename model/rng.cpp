#include "rng.h"

unsigned seed = 12345; // Use current time as seed
std::mt19937 gen(seed); // Seed the generator

// Warm up the generator
void warmup() {
    
    for(int i = 0; i < 10000; ++i) {
        gen();
    }
}

int dummy = (warmup(), 0);