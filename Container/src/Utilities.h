#include <vector>
#include <random>

class Utilities {
public:
    static std::vector<float> createMatrixWave(int x, int y, int z, float frequency, float amplitude, float phase);
    static std::vector<float> createMatrixRandom(int x, int y, int z, float minVal, float maxVal);
};

