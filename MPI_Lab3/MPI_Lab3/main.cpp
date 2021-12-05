#include <stdlib.h>
#include <iostream>
#include "ParallelGaussAlgorithm.h"

int main(int argc, char* argv[]) {
    auto algo = new ParallelGaussAlgorithm();

    algo->execute(argc, argv);

    return 0;
}