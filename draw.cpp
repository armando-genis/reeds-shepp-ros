#include "draw.h"
#include "matplotlibcpp.h"
#include <cmath>
#include <string>
#include <vector>

namespace plt = matplotlibcpp;
namespace draw {

const double PI = 3.14159265358979323846;

Arrow::Arrow(double x, double y, double theta, double L, const std::string& c) {
    double angle = deg2rad(30);
    double d = 0.5 * L;
    double w = 2;

    double x_start = x;
    double y_start = y;
    double x_end = x + L * std::cos(theta);
    double y_end = y + L * std::sin(theta);

    double theta_hat_L = theta + PI - angle;
    double theta_hat_R = theta + PI + angle;

    double x_hat_start = x_end;
    double x_hat_end_L = x_hat_start + d * std::cos(theta_hat_L);
    double x_hat_end_R = x_hat_start + d * std::cos(theta_hat_R);

    double y_hat_start = y_end;
    double y_hat_end_L = y_hat_start + d * std::sin(theta_hat_L);
    double y_hat_end_R = y_hat_start + d * std::sin(theta_hat_R);

    std::map<std::string, std::string> props;
    props["color"] = c;
    // Use a different color for the starting point
    std::map<std::string, std::string> start_point_props;
    start_point_props["color"] = "red"; // Specify the color for the start point

    plt::scatter(std::vector<double>{x_start}, std::vector<double>{y_start}, 10.0, start_point_props);
    plt::scatter(std::vector<double>{x_hat_start, x_hat_end_L}, std::vector<double>{y_hat_start, y_hat_end_L}, 10.0, props);
    plt::scatter(std::vector<double>{x_hat_start, x_hat_end_R}, std::vector<double>{y_hat_start, y_hat_end_R}, 10.0, props);

}

double Arrow::deg2rad(double deg) {
    return deg * PI / 180.0;
}

Car::Car(double x, double y, double yaw, double w, double L) {
    double theta_B = PI + yaw;

    double xB = x + L / 4 * std::cos(theta_B);
    double yB = y + L / 4 * std::sin(theta_B);

    double theta_BL = theta_B + PI / 2;
    double theta_BR = theta_B - PI / 2;

    double x_BL = xB + w / 2 * std::cos(theta_BL);
    double y_BL = yB + w / 2 * std::sin(theta_BL);
    double x_BR = xB + w / 2 * std::cos(theta_BR);
    double y_BR = yB + w / 2 * std::sin(theta_BR);

    double x_FL = x_BL + L * std::cos(yaw);
    double y_FL = y_BL + L * std::sin(yaw);
    double x_FR = x_BR + L * std::cos(yaw);
    double y_FR = y_BR + L * std::sin(yaw);

    // plt::plot({x_BL, x_BR, x_FR, x_FL, x_BL}, {y_BL, y_BR, y_FR, y_FL, y_BL}, "black");

    Arrow(x, y, yaw, L / 2, "black");
}

}  // namespace draw
