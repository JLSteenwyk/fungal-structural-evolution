// Export the frozen FastML library's category rates without changing inference.
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "gammaDistribution.h"
int main() {
    std::string id;
    double alpha;
    std::cout << std::setprecision(17);
    while (std::cin >> id >> alpha) {
        if (!std::isfinite(alpha) || alpha <= 0) return 2;
        gammaDistribution distribution(alpha, 4);
        std::cout << id;
        for (int i=0; i<4; ++i) std::cout << '\t' << distribution.rates(i);
        std::cout << '\n';
    }
    return std::cin.eof() ? 0 : 3;
}
