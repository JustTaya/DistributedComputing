#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include "BlockSchemaAlgorithm.h"

int main(int argc, char* argv[]) {
    auto algo = new BlockSchemaAlgorithm();

    algo -> execute(argc, argv);

    return 0;
}