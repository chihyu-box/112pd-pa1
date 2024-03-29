#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include "partitioner.h"
using namespace std;
using namespace chrono;

int main(int argc, char** argv)
{
    auto start = high_resolution_clock::now();
    fstream input, output;

    if (argc == 3) {
        input.open(argv[1], ios::in);
        output.open(argv[2], ios::out);
        if (!input) {
            cerr << "Cannot open the input file \"" << argv[1]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output) {
            cerr << "Cannot open the output file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
    }
    else {
        cerr << "Usage: ./fm <input file> <output file>" << endl;
        exit(1);
    }

    Partitioner* partitioner = new Partitioner(input);
    partitioner->partition();
    partitioner->printSummary();
    partitioner->writeResult(output);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    cout << "runtime  " << duration.count() << " ms" << endl;

    return 0;
}
