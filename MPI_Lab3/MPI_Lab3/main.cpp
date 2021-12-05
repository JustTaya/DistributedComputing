#include <stdlib.h>
#include <iostream>
#include "SerialAlgorithm.h"

int main(int argc, char* argv[]) {
    auto algo = new SerialAlgorithm();

    algo->execute(argc, argv);

    return 0;
}