// 1D Ising spin system simulation with N spins with ferromagnetic coupling

// including used headers/libraries
#include <iostream>
#include <fstream>
#include <numeric>
#include <math.h>
#include <random>
#include <array>
#include <algorithm>
#include <string>

// size of system
const int N = 50;

// -------------------------------------------------------------------------------------------

// write to file
auto WriteToFile = [](auto &file, auto const &data) {
    for (auto &e : data)
        file << e << " ";
    file << "\n";
};

// -------------------------------------------------------------------------------------------

// generate uniform random real values from given [low, high] interval
auto RandUniReal = [](auto const &low, auto const &high) {
    // random number generation via std
    std::random_device rd{};
    std::mt19937 gen(rd());
    using R = decltype(low / high);
    std::uniform_real_distribution<R> Uni(low, high);
    return Uni(gen);
};

// -------------------------------------------------------------------------------------------

// generate uniform random integer values from given [low, high] interval
auto RandUniInt = [](auto const &low, auto const &high) {
    // random number generation via std
    std::random_device rd{};
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> Uni(low, high);
    return Uni(gen);
};

// -------------------------------------------------------------------------------------------

// fuction to choose random spin direction
auto pm = []() {
    return RandUniReal(0., 1.) > 0.5 ? 1 : -1;
};

// initialize spins
auto Init = [&]() {
    std::array<int, N> vec;
    std::generate(vec.begin(), vec.end(), pm);
    return vec;
};

// initialize spins (ordered: +/- 1)
auto InitOrder = [&](int const &upOrDown) {
    std::array<int, N> vec;
    vec.fill(upOrDown);
    return vec;
};

// -------------------------------------------------------------------------------------------

// "to flip, or not to flip?"
auto ToFlipOrNotToFlip = [&](int const &index, auto &spinState, auto const &betaJ) {
    int si{spinState[index]}, before_si{0}, after_si{0};
    // periodic boundary conditions
    if (index == 0)
        before_si = spinState[N - 1], after_si = spinState[1];
    else if (index == N - 1)
        before_si = spinState[N - 2], after_si = spinState[0];
    else
        before_si = spinState[index - 1], after_si = spinState[index + 1];

    // deltaE
    auto delta = 2 * si * (before_si + after_si);

    // check sign and act accordingly
    if (delta < 0)
        spinState[index] = -spinState[index];
    else if (delta == 0)
        spinState[index] = RandUniReal(0., 1.) < 0.5 ? -spinState[index] : spinState[index];
    else
        spinState[index] = RandUniReal(0., 1.) < std::exp(-betaJ * delta) ? -spinState[index] : spinState[index];
};

// -------------------------------------------------------------------------------------------

// magnetization
auto Magnetizaion = [&](auto const &spinState) {
    return std::accumulate(spinState.begin(), spinState.end(), 0.) / N;
};

// -------------------------------------------------------------------------------------------

// main function
// argv[1] --> beta * J
// argv[2] --> t
// argv[3] --> file name
int main(int argc, char **argv)
{
    // check argument list
    if (argc < 4)
    {
        std::cout << "ERROR\nNot enough parameter." << std::endl;
        std::exit(-1);
    }

    // declare beta * J parameter
    double betaJ{std::stod(argv[1])};

    // declare integration time
    int t{std::stoi(argv[2])};

    // declare file name
    std::string fileName = argv[3];

    // initialize spins
    //std::array<int, N> spinState(Init());
    std::array<int, N> spinState(InitOrder(1));

    // save data
    std::ofstream file;
    file.open(fileName);
    // simulation
    int i{0};
    while (i < N * t)
    {
        double m = Magnetizaion(spinState);
        file << m << "\n";
        int index{RandUniInt(0, N)};
        ToFlipOrNotToFlip(index, spinState, betaJ);
        //WriteToFile(file, spinState);
        i++;
    }
    file.close();
}
