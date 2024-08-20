#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <queue>
#include "dubins_path.h" // Include the Dubins path header

class Car {
public:
    static constexpr double maxSteerAngle = 0.6;
    static constexpr int steerPrecision = 10;
    static constexpr double wheelBase = 3.5;
    static constexpr double axleToFront = 4.5;
    static constexpr double axleToBack = 1.0;
    static constexpr double width = 3.0;
};

class Cost {
public:
    static constexpr double reverse = 10.0;
    static constexpr double directionChange = 150.0;
    static constexpr double steerAngle = 1.0;
    static constexpr double steerAngleChange = 5.0;
    static constexpr double hybridCost = 50.0;
};

class Node {
public:
    std::tuple<int, int, int> gridIndex;  // grid block x, y, yaw index
    std::vector<std::tuple<double, double, double>> traj;  // trajectory x, y, yaw of a simulated node
    double steeringAngle;  // steering angle throughout the trajectory
    int direction;  // direction throughout the trajectory
    double cost;  // node cost
    std::tuple<int, int, int> parentIndex;  // parent node index

    // Constructor
    Node(std::tuple<int, int, int> gridIdx, 
         const std::vector<std::tuple<double, double, double>>& trajectory, 
         double steerAngle, 
         int dir, 
         double nodeCost, 
         std::tuple<int, int, int> parentIdx)
        : gridIndex(gridIdx), 
          traj(trajectory), 
          steeringAngle(steerAngle), 
          direction(dir), 
          cost(nodeCost), 
          parentIndex(parentIdx) 
    {}
};



int main() {



    return 0;
}
