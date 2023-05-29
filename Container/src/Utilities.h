#include <vector>
#include <random>

class Utilities {
public:
    static std::vector<double> createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase);
    static std::vector<double> createMatrixRandom(int x, int y, int z, double minVal, double maxVal);
};

