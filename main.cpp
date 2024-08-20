#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include <thread>
#include <vector>

#include "matplotlibcpp.h"  // Include the matplotlib-cpp header
#include "dubins_path.h"    // Include the Dubins path header


using namespace std;
namespace plt = matplotlibcpp;

int main() 
{


    std::vector<std::tuple<double, double, double>> states = {
        {-3, 3, 120}, {10, -7, -30}, {10, 13, 30}, {20, 5, -25},
        {35, 10, 180}, {32, -10, 180}, {5, -12, 90}
    };

    double max_c = 0.25;
    std::vector<double> path_x, path_y, yaw;

    // Create an instance of dubins_path
    dubins_path dp;


    // Calculate the path for each segment
    for (size_t i = 0; i < states.size() - 1; ++i) {
        double s_x = std::get<0>(states[i]);
        double s_y = std::get<1>(states[i]);
        double s_yaw = std::get<2>(states[i]) * M_PI / 180.0;  // Convert degrees to radians
        double g_x = std::get<0>(states[i + 1]);
        double g_y = std::get<1>(states[i + 1]);
        double g_yaw = std::get<2>(states[i + 1]) * M_PI / 180.0;  // Convert degrees to radians

        std::cout << "Processing segment " << i << ": Start (" << s_x << ", " << s_y << ", " << s_yaw << ") -> "
                  << "Goal (" << g_x << ", " << g_y << ", " << g_yaw << ")" << std::endl;

        // Calculate the Dubins path for this segment
        Path path_i = dp.calc_dubins_path(s_x, s_y, s_yaw, g_x, g_y, g_yaw, max_c);

        if (!path_i.x.empty() && !path_i.y.empty() && !path_i.yaw.empty()) {
            path_x.insert(path_x.end(), path_i.x.begin(), path_i.x.end());
            path_y.insert(path_y.end(), path_i.y.begin(), path_i.y.end());
            yaw.insert(yaw.end(), path_i.yaw.begin(), path_i.yaw.end());
        } else {
            std::cerr << "Invalid path generated for segment " << i << std::endl;
        }

        std::cout << "Path segment " << i << " processed. Path length: " << path_i.x.size() << std::endl;
    }

    // Plot the path
    plt::plot(path_x, path_y, "gray");
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::title("Waypoints and Vehicle Path with Orientation");
    plt::grid(true);

    // Display the plot
    plt::show();

    return 0;
}


// Compile and run the code
// g++ -I/usr/include/eigen3 -I/usr/include/python3.10 main.cpp dubins_path.cp -o dubins_path -lpython3.10 -g