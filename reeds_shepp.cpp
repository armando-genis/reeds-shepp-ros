#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>
#include <map>
#include "matplotlibcpp.h"
#include "draw.h"

namespace plt = matplotlibcpp;


const double STEP_SIZE = 0.2;
const double MAX_LENGTH = 1000.0;
const double PI = M_PI;

class PATH {
public:
    std::vector<double> lengths;      // lengths of each part of path
    std::vector<std::string> ctypes;  // type of each part of the path
    double L;                         // total path length
    std::vector<double> x;            // final x positions
    std::vector<double> y;            // final y positions
    std::vector<double> yaw;          // final yaw angles
    std::vector<int> directions;      // forward: 1, backward: -1

    // Constructor
    PATH(const std::vector<double>& lengths,
         const std::vector<std::string>& ctypes,
         double L,
         const std::vector<double>& x,
         const std::vector<double>& y,
         const std::vector<double>& yaw,
         const std::vector<int>& directions)
        : lengths(lengths), ctypes(ctypes), L(L), x(x), y(y), yaw(yaw), directions(directions) {}
};


// Forward declarations
std::vector<PATH> calc_all_paths(double sx, double sy, double syaw,
                                 double gx, double gy, double gyaw,
                                 double maxc, double step_size);

std::vector<PATH> generate_path(const std::array<double, 3>& q0, const std::array<double, 3>& q1, double maxc);

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int>>
generate_local_course(double L, const std::vector<double>& lengths, const std::vector<char>& mode,
                      double maxc, double step_size);

double pi_2_pi(double theta);

std::tuple<double, double> R(double x, double y);

double M(double theta);

std::tuple<bool, double, double, double> SLS(double x, double y, double phi);

std::tuple<bool, double, double, double> LRLRn(double x, double y, double phi);

std::tuple<bool, double, double, double> LRLRp(double x, double y, double phi);

std::vector<PATH> SCS(double x, double y, double phi, std::vector<PATH> paths);

std::vector<PATH> CSC(double x, double y, double phi, std::vector<PATH> paths);

std::vector<PATH> CCC(double x, double y, double phi, std::vector<PATH> paths);

std::vector<PATH> CCCC(double x, double y, double phi, std::vector<PATH> paths);

std::tuple<bool, double, double, double> LSL(double x, double y, double phi);

std::tuple<bool, double, double, double> LSR(double x, double y, double phi);

std::tuple<bool, double, double, double> LRL(double x, double y, double phi);

std::vector<PATH> CCSC(double x, double y, double phi, std::vector<PATH> paths);

std::tuple<bool, double, double, double> LRSR(double x, double y, double phi);

std::tuple<bool, double, double, double> LRSL(double x, double y, double phi);

std::vector<PATH> CCSCC(double x, double y, double phi, std::vector<PATH> paths);

std::tuple<bool, double, double, double> LRSLR(double x, double y, double phi);

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int>> 
interpolate(int ind, double l, char m, double maxc, double ox, double oy, double oyaw,
            std::vector<double>& px, std::vector<double>& py, std::vector<double>& pyaw, std::vector<int>& directions);

// Reeds-Shepp path types

PATH calc_optimal_path(double sx, double sy, double syaw, 
                       double gx, double gy, double gyaw, 
                       double maxc, double step_size = STEP_SIZE) {

    // Print the sx and sy values
    std::cout << "sx: " << sx << std::endl;
    std::cout << "sy: " << sy << std::endl;
    std::cout << "syaw: " << syaw << std::endl;
    
    std::vector<PATH> paths = calc_all_paths(sx, sy, syaw, gx, gy, gyaw, maxc, step_size);

    double minL = paths[0].L;
    std::size_t mini = 0;

    for (std::size_t i = 0; i < paths.size(); ++i) {
        if (paths[i].L <= minL) {
            minL = paths[i].L;
            mini = i;
        }
    }

    return paths[mini];
}


std::vector<PATH> calc_all_paths(double sx, double sy, double syaw, 
                                 double gx, double gy, double gyaw, 
                                 double maxc, double step_size = STEP_SIZE) 
{
    std::array<double, 3> q0 = {sx, sy, syaw};
    std::array<double, 3> q1 = {gx, gy, gyaw};

    std::vector<PATH> paths = generate_path(q0, q1, maxc);

    // cout the path ctype
    std::cout << "path ctypes: ";
    for (const auto& ctype : paths[0].ctypes) {
        std::cout << ctype << " ";
    }
    std::cout << std::endl;

    // print the q1 and q0 values
    std::cout << "q0: " << q0[0] << " " << q0[1] << " " << q0[2] << std::endl;
    std::cout << "q1: " << q1[0] << " " << q1[1] << " " << q1[2] << std::endl;

    for (auto& path : paths) {
        std::vector<double> x, y, yaw;
        std::vector<int> directions;

        // Convert ctypes (strings) to mode (characters)
        std::vector<char> mode;
        for (const auto& ctype : path.ctypes) {
            if (!ctype.empty()) {
                mode.push_back(ctype[0]);  // Take the first character of each string
            }
        }

        std::tie(x, y, yaw, directions) = generate_local_course(path.L, path.lengths, 
                                                                mode, maxc, step_size * maxc);

        // print path length
        std::cout << "path length: " << path.L << std::endl;


        // Convert global coordinates
        path.x.resize(x.size());
        path.y.resize(y.size());
        path.yaw.resize(yaw.size());

        for (size_t i = 0; i < x.size(); ++i) {
            path.x[i] = cos(-q0[2]) * x[i] + sin(-q0[2]) * y[i] + q0[0];
            path.y[i] = -sin(-q0[2]) * x[i] + cos(-q0[2]) * y[i] + q0[1];
            path.yaw[i] = pi_2_pi(yaw[i] + q0[2]);
        }

        path.directions = directions;
        
        for (auto& l : path.lengths) {
            l /= maxc;
        }

        path.L /= maxc;
    }

    // cout final path lenght 
    std::cout << "final path length: " << paths[0].L << std::endl;

    // cout ceparetor
    std::cout << "--------------------------------" << std::endl;

    return paths;
}
    



std::vector<PATH> set_path(std::vector<PATH>& paths, 
                           const std::vector<double>& lengths, 
                           const std::vector<std::string>& ctypes) {
    PATH path({}, ctypes, 0.0, {}, {}, {}, {});
    path.ctypes = ctypes;
    path.lengths = lengths;

    // Check if the same path exists
    for (const auto& path_e : paths) {
        if (path_e.ctypes == path.ctypes) {
            double sum_diff = 0.0;
            for (size_t i = 0; i < path_e.lengths.size(); ++i) {
                sum_diff += std::abs(path_e.lengths[i] - path.lengths[i]);
            }
            if (sum_diff <= 0.01) {
                return paths;  // Not inserting path
            }
        }
    }

    path.L = std::accumulate(lengths.begin(), lengths.end(), 0.0, [](double acc, double val) {
        return acc + std::abs(val);
    });

    if (path.L >= MAX_LENGTH) {
        return paths;
    }

    assert(path.L >= 0.01);
    paths.push_back(path);

    return paths;
}

std::tuple<bool, double, double, double> LSL(double x, double y, double phi) {
    auto [u, t] = R(x - std::sin(phi), y - 1.0 + std::cos(phi));

    if (t >= 0.0) {
        double v = M(phi - t);
        if (v >= 0.0) {
            return {true, t, u, v};
        }
    }

    return {false, 0.0, 0.0, 0.0};
}


std::tuple<bool, double, double, double> LSR(double x, double y, double phi) {
    double u1, t1;
    std::tie(u1, t1) = R(x + std::sin(phi), y - 1.0 - std::cos(phi));
    u1 = u1 * u1;

    if (u1 >= 4.0) {
        double u = std::sqrt(u1 - 4.0);
        double theta = std::atan2(2.0, u);
        double t = M(t1 + theta);
        double v = M(t - phi);

        if (t >= 0.0 && v >= 0.0) {
            return std::make_tuple(true, t, u, v);
        }
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}

std::tuple<bool, double, double, double> LRL(double x, double y, double phi) {
    double u1, t1;
    std::tie(u1, t1) = R(x - std::sin(phi), y - 1.0 + std::cos(phi));

    if (u1 <= 4.0) {
        double u = -2.0 * std::asin(0.25 * u1);
        double t = M(t1 + 0.5 * u + PI);
        double v = M(phi - t + u);

        if (t >= 0.0 && u <= 0.0) {
            return std::make_tuple(true, t, u, v);
        }
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}


std::vector<PATH> SCS(double x, double y, double phi, std::vector<PATH> paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = SLS(x, y, phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"S", "L", "S"});

    std::tie(flag, t, u, v) = SLS(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"S", "R", "S"});

    return paths;
}

std::tuple<bool, double, double, double> SLS(double x, double y, double phi) {
    phi = M(phi);

    if (y > 0.0 && 0.0 < phi && phi < PI * 0.99) {
        double xd = -y / std::tan(phi) + x;
        double t = xd - std::tan(phi / 2.0);
        double u = phi;
        double v = std::sqrt((x - xd) * (x - xd) + y * y) - std::tan(phi / 2.0);
        return std::make_tuple(true, t, u, v);
    } else if (y < 0.0 && 0.0 < phi && phi < PI * 0.99) {
        double xd = -y / std::tan(phi) + x;
        double t = xd - std::tan(phi / 2.0);
        double u = phi;
        double v = -std::sqrt((x - xd) * (x - xd) + y * y) - std::tan(phi / 2.0);
        return std::make_tuple(true, t, u, v);
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}

std::vector<PATH> CSC(double x, double y, double phi, std::vector<PATH> paths) {
    bool flag;
    double t, u, v;

    // LSL
    std::tie(flag, t, u, v) = LSL(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, u, v}, {"L", "S", "L"});
    }

    std::tie(flag, t, u, v) = LSL(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -v}, {"L", "S", "L"});
    }

    std::tie(flag, t, u, v) = LSL(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, u, v}, {"R", "S", "R"});
    }

    std::tie(flag, t, u, v) = LSL(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -v}, {"R", "S", "R"});
    }

    // LSR
    std::tie(flag, t, u, v) = LSR(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, u, v}, {"L", "S", "R"});
    }

    std::tie(flag, t, u, v) = LSR(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -v}, {"L", "S", "R"});
    }

    std::tie(flag, t, u, v) = LSR(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, u, v}, {"R", "S", "L"});
    }

    std::tie(flag, t, u, v) = LSR(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -v}, {"R", "S", "L"});
    }

    return paths;
}


std::vector<PATH> CCC(double x, double y, double phi, std::vector<PATH> paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRL(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, u, v}, {"L", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRL(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -v}, {"L", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRL(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, u, v}, {"R", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRL(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -v}, {"R", "L", "R"});
    }

    // backwards
    double xb = x * std::cos(phi) + y * std::sin(phi);
    double yb = x * std::sin(phi) - y * std::cos(phi);

    std::tie(flag, t, u, v) = LRL(xb, yb, phi);
    if (flag) {
        paths = set_path(paths, {v, u, t}, {"L", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRL(-xb, yb, -phi);
    if (flag) {
        paths = set_path(paths, {-v, -u, -t}, {"L", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRL(xb, -yb, -phi);
    if (flag) {
        paths = set_path(paths, {v, u, t}, {"R", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRL(-xb, -yb, phi);
    if (flag) {
        paths = set_path(paths, {-v, -u, -t}, {"R", "L", "R"});
    }

    return paths;
}


std::tuple<double, double> calc_tauOmega(double u, double v, double xi, double eta, double phi) {
    double delta = M(u - v);
    double A = std::sin(u) - std::sin(delta);
    double B = std::cos(u) - std::cos(delta) - 1.0;

    double t1 = std::atan2(eta * A - xi * B, xi * A + eta * B);
    double t2 = 2.0 * (std::cos(delta) - std::cos(v) - std::cos(u)) + 3.0;

    double tau = t2 < 0 ? M(t1 + PI) : M(t1);
    double omega = M(tau - u + v - phi);

    return std::make_tuple(tau, omega);
}

std::tuple<bool, double, double, double> LRLRn(double x, double y, double phi) {
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    double rho = 0.25 * (2.0 + std::sqrt(xi * xi + eta * eta));

    if (rho <= 1.0) {
        double u = std::acos(rho);
        double t, v;
        std::tie(t, v) = calc_tauOmega(u, -u, xi, eta, phi);

        if (t >= 0.0 && v <= 0.0) {
            return std::make_tuple(true, t, u, v);
        }
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}

std::tuple<bool, double, double, double> LRLRp(double x, double y, double phi) {
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    double rho = (20.0 - xi * xi - eta * eta) / 16.0;

    if (rho >= 0.0 && rho <= 1.0) {
        double u = -std::acos(rho);
        if (u >= -0.5 * PI) {
            double t, v;
            std::tie(t, v) = calc_tauOmega(u, u, xi, eta, phi);

            if (t >= 0.0 && v >= 0.0) {
                return std::make_tuple(true, t, u, v);
            }
        }
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}



std::vector<PATH> CCCC(double x, double y, double phi, std::vector<PATH> paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRLRn(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, u, -u, v}, {"L", "R", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRLRn(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, u, -v}, {"L", "R", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRLRn(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, u, -u, v}, {"R", "L", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRLRn(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, u, -v}, {"R", "L", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRLRp(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, u, u, v}, {"L", "R", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRLRp(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -u, -v}, {"L", "R", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRLRp(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, u, u, v}, {"R", "L", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRLRp(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, -u, -u, -v}, {"R", "L", "R", "L"});
    }

    return paths;
}

std::tuple<bool, double, double, double> LRSR(double x, double y, double phi) {
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    double rho, theta;
    std::tie(rho, theta) = R(-eta, xi);

    if (rho >= 2.0) {
        double t = theta;
        double u = 2.0 - rho;
        double v = M(t + 0.5 * PI - phi);
        if (t >= 0.0 && u <= 0.0 && v <= 0.0) {
            return std::make_tuple(true, t, u, v);
        }
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}

std::tuple<bool, double, double, double> LRSL(double x, double y, double phi) {
    double xi = x - std::sin(phi);
    double eta = y - 1.0 + std::cos(phi);
    double rho, theta;
    std::tie(rho, theta) = R(xi, eta);

    if (rho >= 2.0) {
        double r = std::sqrt(rho * rho - 4.0);
        double u = 2.0 - r;
        double t = M(theta + std::atan2(r, -2.0));
        double v = M(phi - 0.5 * PI - t);
        if (t >= 0.0 && u <= 0.0 && v <= 0.0) {
            return std::make_tuple(true, t, u, v);
        }
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}


std::vector<PATH> CCSC(double x, double y, double phi, std::vector<PATH> paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRSL(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, -0.5 * PI, u, v}, {"L", "R", "S", "L"});
    }

    std::tie(flag, t, u, v) = LRSL(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"L", "R", "S", "L"});
    }

    std::tie(flag, t, u, v) = LRSL(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, -0.5 * PI, u, v}, {"R", "L", "S", "R"});
    }

    std::tie(flag, t, u, v) = LRSL(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"R", "L", "S", "R"});
    }

    std::tie(flag, t, u, v) = LRSR(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, -0.5 * PI, u, v}, {"L", "R", "S", "R"});
    }

    std::tie(flag, t, u, v) = LRSR(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"L", "R", "S", "R"});
    }

    std::tie(flag, t, u, v) = LRSR(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, -0.5 * PI, u, v}, {"R", "L", "S", "L"});
    }

    std::tie(flag, t, u, v) = LRSR(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"R", "L", "S", "L"});
    }

    // Backwards
    double xb = x * std::cos(phi) + y * std::sin(phi);
    double yb = x * std::sin(phi) - y * std::cos(phi);

    std::tie(flag, t, u, v) = LRSL(xb, yb, phi);
    if (flag) {
        paths = set_path(paths, {v, u, -0.5 * PI, t}, {"L", "S", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRSL(-xb, yb, -phi);
    if (flag) {
        paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"L", "S", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRSL(xb, -yb, -phi);
    if (flag) {
        paths = set_path(paths, {v, u, -0.5 * PI, t}, {"R", "S", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRSL(-xb, -yb, phi);
    if (flag) {
        paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"R", "S", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRSR(xb, yb, phi);
    if (flag) {
        paths = set_path(paths, {v, u, -0.5 * PI, t}, {"R", "S", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRSR(-xb, yb, -phi);
    if (flag) {
        paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"R", "S", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRSR(xb, -yb, -phi);
    if (flag) {
        paths = set_path(paths, {v, u, -0.5 * PI, t}, {"L", "S", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRSR(-xb, -yb, phi);
    if (flag) {
        paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"L", "S", "L", "R"});
    }

    return paths;
}

std::tuple<bool, double, double, double> LRSLR(double x, double y, double phi) {
    // formula 8.11 *** TYPO IN PAPER ***
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    double rho, theta;
    std::tie(rho, theta) = R(xi, eta);

    if (rho >= 2.0) {
        double u = 4.0 - std::sqrt(rho * rho - 4.0);
        if (u <= 0.0) {
            double t = M(std::atan2((4.0 - u) * xi - 2.0 * eta, -2.0 * xi + (u - 4.0) * eta));
            double v = M(t - phi);

            if (t >= 0.0 && v >= 0.0) {
                return std::make_tuple(true, t, u, v);
            }
        }
    }

    return std::make_tuple(false, 0.0, 0.0, 0.0);
}

std::vector<PATH> CCSCC(double x, double y, double phi, std::vector<PATH> paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRSLR(x, y, phi);
    if (flag) {
        paths = set_path(paths, {t, -0.5 * PI, u, -0.5 * PI, v}, {"L", "R", "S", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRSLR(-x, y, -phi);
    if (flag) {
        paths = set_path(paths, {-t, 0.5 * PI, -u, 0.5 * PI, -v}, {"L", "R", "S", "L", "R"});
    }

    std::tie(flag, t, u, v) = LRSLR(x, -y, -phi);
    if (flag) {
        paths = set_path(paths, {t, -0.5 * PI, u, -0.5 * PI, v}, {"R", "L", "S", "R", "L"});
    }

    std::tie(flag, t, u, v) = LRSLR(-x, -y, phi);
    if (flag) {
        paths = set_path(paths, {-t, 0.5 * PI, -u, 0.5 * PI, -v}, {"R", "L", "S", "R", "L"});
    }

    return paths;
}


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int>> 
generate_local_course(double L, const std::vector<double>& lengths, const std::vector<char>& mode, 
                      double maxc, double step_size) {
    
    int point_num = static_cast<int>(L / step_size) + lengths.size() + 3;

    std::vector<double> px(point_num, 0.0);
    std::vector<double> py(point_num, 0.0);
    std::vector<double> pyaw(point_num, 0.0);
    std::vector<int> directions(point_num, 0);

    int ind = 1;

    directions[0] = (lengths[0] > 0.0) ? 1 : -1;
    double d = (lengths[0] > 0.0) ? step_size : -step_size;

    double pd = d;
    double ll = 0.0;

    for (size_t i = 0; i < mode.size(); ++i) {
        double l = lengths[i];
        char m = mode[i];

        d = (l > 0.0) ? step_size : -step_size;

        double ox = px[ind], oy = py[ind], oyaw = pyaw[ind];
        --ind;

        if (i >= 1 && (lengths[i - 1] * lengths[i]) > 0) {
            pd = -d - ll;
        } else {
            pd = d - ll;
        }

        while (std::abs(pd) <= std::abs(l)) {
            ++ind;
            std::tie(px, py, pyaw, directions) = 
                interpolate(ind, pd, m, maxc, ox, oy, oyaw, px, py, pyaw, directions);
            pd += d;
        }

        ll = l - pd - d;  // calc remain length

        ++ind;
        std::tie(px, py, pyaw, directions) = 
            interpolate(ind, l, m, maxc, ox, oy, oyaw, px, py, pyaw, directions);
    }

    // remove unused data
    while (!px.empty() && px.back() == 0.0) {
        px.pop_back();
        py.pop_back();
        pyaw.pop_back();
        directions.pop_back();
    }

    return std::make_tuple(px, py, pyaw, directions);
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int>> 
interpolate(int ind, double l, char m, double maxc, double ox, double oy, double oyaw,
            std::vector<double>& px, std::vector<double>& py, std::vector<double>& pyaw, std::vector<int>& directions) {
    
    if (m == 'S') {
        px[ind] = ox + l / maxc * std::cos(oyaw);
        py[ind] = oy + l / maxc * std::sin(oyaw);
        pyaw[ind] = oyaw;
    } else {
        double ldx = std::sin(l) / maxc;
        double ldy;

        if (m == 'L') {
            ldy = (1.0 - std::cos(l)) / maxc;
        } else if (m == 'R') {
            ldy = (1.0 - std::cos(l)) / (-maxc);
        }

        double gdx = std::cos(-oyaw) * ldx + std::sin(-oyaw) * ldy;
        double gdy = -std::sin(-oyaw) * ldx + std::cos(-oyaw) * ldy;
        px[ind] = ox + gdx;
        py[ind] = oy + gdy;
    }

    if (m == 'L') {
        pyaw[ind] = oyaw + l;
    } else if (m == 'R') {
        pyaw[ind] = oyaw - l;
    }

    directions[ind] = (l > 0.0) ? 1 : -1;

    return std::make_tuple(px, py, pyaw, directions);
}

std::vector<PATH> generate_path(const std::array<double, 3>& q0, 
                                const std::array<double, 3>& q1, 
                                double maxc) 
{
    double dx = q1[0] - q0[0];
    double dy = q1[1] - q0[1];
    double dth = q1[2] - q0[2];
    double c = std::cos(q0[2]);
    double s = std::sin(q0[2]);
    double x = (c * dx + s * dy) * maxc;
    double y = (-s * dx + c * dy) * maxc;

    std::vector<PATH> paths;

    paths = SCS(x, y, dth, paths);
    paths = CSC(x, y, dth, paths);
    paths = CCC(x, y, dth, paths);
    paths = CCCC(x, y, dth, paths);
    paths = CCSC(x, y, dth, paths);
    paths = CCSCC(x, y, dth, paths);

    return paths;
}

// Utils
// Function to normalize theta to the range [-PI, PI]
double pi_2_pi(double theta) {
    while (theta > PI) {
        theta -= 2.0 * PI;
    }

    while (theta < -PI) {
        theta += 2.0 * PI;
    }

    return theta;
}

// Function to return the polar coordinates (r, theta) of the point (x, y)
std::tuple<double, double> R(double x, double y) {
    double r = std::hypot(x, y);
    double theta = std::atan2(y, x);

    return std::make_tuple(r, theta);
}

// Function to regulate theta to the range [-PI, PI]
double M(double theta) {
    double phi = std::fmod(theta, 2.0 * PI);

    if (phi < -PI) {
        phi += 2.0 * PI;
    }
    if (phi > PI) {
        phi -= 2.0 * PI;
    }

    return phi;
}

std::string get_label(const PATH& path) {
    std::string label;

    for (std::size_t i = 0; i < path.ctypes.size(); ++i) {
        label += path.ctypes[i];
        if (path.lengths[i] > 0.0) {
            label += "+";
        } else {
            label += "-";
        }
    }

    return label;
}



std::tuple<std::vector<double>, std::vector<double>> calc_curvature(
    const std::vector<double>& x, 
    const std::vector<double>& y, 
    const std::vector<double>& yaw, 
    const std::vector<int>& directions) {

    std::vector<double> c, ds;

    for (std::size_t i = 1; i < x.size() - 1; ++i) {
        double dxn = x[i] - x[i - 1];
        double dxp = x[i + 1] - x[i];
        double dyn = y[i] - y[i - 1];
        double dyp = y[i + 1] - y[i];
        double dn = std::hypot(dxn, dyn);
        double dp = std::hypot(dxp, dyp);

        double dx = 1.0 / (dn + dp) * (dp / dn * dxn + dn / dp * dxp);
        double ddx = 2.0 / (dn + dp) * (dxp / dp - dxn / dn);
        double dy = 1.0 / (dn + dp) * (dp / dn * dyn + dn / dp * dyp);
        double ddy = 2.0 / (dn + dp) * (dyp / dp - dyn / dn);

        double curvature = (ddy * dx - ddx * dy) / (dx * dx + dy * dy);
        double d = (dn + dp) / 2.0;

        if (std::isnan(curvature)) {
            curvature = 0.0;
        }

        if (directions[i] <= 0.0) {
            curvature = -curvature;
        }

        if (c.empty()) {
            ds.push_back(d);
            c.push_back(curvature);
        }

        ds.push_back(d);
        c.push_back(curvature);
    }

    if (!ds.empty()) {
        ds.push_back(ds.back());
        c.push_back(c.back());
    }

    return std::make_tuple(c, ds);
}


void check_path(double sx, double sy, double syaw, 
                double gx, double gy, double gyaw, 
                double maxc) {
    std::vector<PATH> paths = calc_all_paths(sx, sy, syaw, gx, gy, gyaw, maxc);

    assert(!paths.empty());

    for (const auto& path : paths) {
        assert(std::abs(path.x[0] - sx) <= 0.01);
        assert(std::abs(path.y[0] - sy) <= 0.01);
        assert(std::abs(path.yaw[0] - syaw) <= 0.01);
        assert(std::abs(path.x.back() - gx) <= 0.01);
        assert(std::abs(path.y.back() - gy) <= 0.01);
        assert(std::abs(path.yaw.back() - gyaw) <= 0.01);

        std::vector<double> dx(path.x.size() - 1);
        std::vector<double> dy(path.y.size() - 1);

        // Calculate differences between consecutive points
        std::transform(path.x.begin() + 1, path.x.end(), path.x.begin(), dx.begin(), std::minus<double>());
        std::transform(path.y.begin() + 1, path.y.end(), path.y.begin(), dy.begin(), std::minus<double>());

        std::vector<double> d(dx.size());

        for (std::size_t i = 0; i < dx.size(); ++i) {
            d[i] = std::hypot(dx[i], dy[i]);
            assert(std::abs(d[i] - STEP_SIZE) <= 0.001);
        }
    }
}

void main_simulation() 
{


    std::vector<std::tuple<double, double, double>> states = {
        {0, 0, 0},
        {10, 0, 0},
        {20, 13, 30},
        {20, 5, -25},
        {35, 10, 120},
        {32, -10, 180},
        {5, -12, 90}
    };

    double max_c = 0.1;  // max curvature
    std::vector<double> path_x, path_y, yaw;
    


    for (size_t i = 0; i < states.size() - 1; ++i) {
        double s_x = std::get<0>(states[i]);
        double s_y = std::get<1>(states[i]);
        double s_yaw = std::get<2>(states[i]) * M_PI / 180.0;  // Convert to radians

        double g_x = std::get<0>(states[i + 1]);
        double g_y = std::get<1>(states[i + 1]);
        double g_yaw = std::get<2>(states[i + 1]) * M_PI / 180.0;  // Convert to radians



        PATH path_i = calc_optimal_path(s_x, s_y, s_yaw, g_x, g_y, g_yaw, max_c);

        path_x.insert(path_x.end(), path_i.x.begin(), path_i.x.end());
        path_y.insert(path_y.end(), path_i.y.begin(), path_i.y.end());
        yaw.insert(yaw.end(), path_i.yaw.begin(), path_i.yaw.end());
    }

    std::cout << "Path length: " << path_x.size() << std::endl;
    std::cout << "Path length: " << path_y.size() << std::endl;
    std::cout << "Path length: " << yaw.size() << std::endl;
    

    // Plot the states using Arrow objects
    for (const auto& state : states) {
        double x = std::get<0>(state);
        double y = std::get<1>(state);
        double theta = std::get<2>(state) * M_PI / 180.0;  // Convert to radians
        draw::Arrow arrow(x, y, theta, 2, "blueviolet");
    }


    // Plot the path
    plt::plot(path_x, path_y, "gray");
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::title("Waypoints and Vehicle Path with Orientation");
    plt::grid(true);

    // Display the plot
    plt::show();

    // // Begin animation
    // plt::ion();  // Interactive mode on
    // plt::figure(1);

    // for (size_t i = 0; i < path_x.size(); ++i) {
    //     plt::clf();  // Clear the current figure
    //     plt::plot(path_x, path_y, "gray");

    //     for (const auto& state : states) {
    //         double x = std::get<0>(state);
    //         double y = std::get<1>(state);
    //         double theta = std::get<2>(state) * M_PI / 180.0;  // Convert to radians
    //         draw::Arrow arrow(x, y, theta, 2, "blueviolet");
    //     }

    //     draw::Car car(path_x[i], path_y[i], yaw[i], 1.5, 3);

    //     plt::axis("equal");

    //     try {
    //         plt::title("Simulation of Reeds-Shepp Curves");  // Wrap in try-catch block
    //     } catch (const std::runtime_error& e) {
    //         std::cerr << "Runtime error while setting title: " << e.what() << std::endl;
    //     }

    //     plt::axis({-10, 42, -20, 20});
    //     plt::draw();  // Update the figure
    //     plt::pause(0.001);  // Pause for a short time
    // }

    // plt::pause(1);  // Pause before closing the figure
    // plt::show();    // Keep the plot window open

}

int main() {
    main_simulation();
    return 0;
}


// std::vector<PATH> calc_all_paths(double sx, double sy, double syaw,
//                                  double gx, double gy, double gyaw,
//                                  double maxc, double step_size = STEP_SIZE)
// {
    
//     std::array<double, 3> q0 = {sx, sy, syaw};
//     std::array<double, 3> q1 = {gx, gy, gyaw};


//     std::vector<PATH> paths = generate_path(q0, q1, maxc);
//     std::vector<double> x, y, yaw;
//     std::vector<int> directions;

//     // cout the path ctype
//     std::cout << "path ctypes: ";
//     for (const auto& ctype : paths[0].ctypes) {
//         std::cout << ctype << " ";
//     }
//     std::cout << std::endl;

//     // print the q1 and q0 values
//     std::cout << "q0: " << q0[0] << " " << q0[1] << " " << q0[2] << std::endl;
//     std::cout << "q1: " << q1[0] << " " << q1[1] << " " << q1[2] << std::endl;

//     for (auto& path : paths) {



//         // Convert ctypes (strings) to mode (characters)
//         std::vector<char> mode;
//         for (const auto& ctype : path.ctypes) {
//             if (!ctype.empty()) {
//                 mode.push_back(ctype[0]);  // Take the first character of each string
//             }
//         }

        

//         std::tie(x, y, yaw, directions) = generate_local_course(path.L, path.lengths, 
//                                                                 mode, maxc, step_size * maxc);

//         // print x lenght 
//         // std::cout << "x length: " << x.size() << std::endl;

//         // print path length
//         std::cout << "path length: " << path.L << std::endl;



//         // Convert global coordinates
//         path.x = {cos(-q0[2]) * ix + sin(-q0[2]) * iy + q0[0] for (auto ix, iy) : zip(x, y)};
//         path.y = {-sin(-q0[2]) * ix + cos(-q0[2]) * iy + q0[1] for (auto ix, iy) : zip(x, y)};
//         path.yaw = {pi_2_pi(iyaw + q0[2]) for (auto iyaw) : yaw};
//         path.directions = directions;
//         path.lengths = {l / maxc for (auto l) : path.lengths};
//         path.L = path.L / maxc;
//     }

//     // cout final path lenght 
//     std::cout << "final path length: " << paths[0].L << std::endl;

//     // cout ceparetor
//     std::cout << "--------------------------------" << std::endl;

//     return paths;
// }