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

class Path {
public:
    std::vector<double> lengths;
    std::vector<std::string> ctypes;
    double L;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> yaw;
    std::vector<int> directions;

    Path(std::vector<double> lengths, std::vector<std::string> ctypes, double L, 
         std::vector<double> x, std::vector<double> y, std::vector<double> yaw, 
         std::vector<int> directions)
        : lengths(lengths), ctypes(ctypes), L(L), x(x), y(y), yaw(yaw), directions(directions) {}
};

// Forward declarations without default arguments
std::vector<Path> calc_all_paths(double sx, double sy, double syaw, double gx, double gy, double gyaw, double maxc, double step_size);
std::vector<Path> generate_path(const std::vector<double>& q0, const std::vector<double>& q1, double maxc);
std::pair<std::vector<double>, std::vector<double>> generate_local_course(double L, const std::vector<double>& lengths, 
                                                                          const std::vector<std::string>& mode, double maxc, double step_size);
void interpolate(int ind, double l, const std::string& m, double maxc, double ox, double oy, double oyaw, 
                 std::vector<double>& px, std::vector<double>& py, std::vector<double>& pyaw, std::vector<int>& directions);
std::tuple<bool, double, double, double> SLS(double x, double y, double phi);


double pi_2_pi(double theta) {
    while (theta > PI) {
        theta -= 2.0 * PI;
    }
    while (theta < -PI) {
        theta += 2.0 * PI;
    }
    return theta;
}

std::pair<double, double> R(double x, double y) {
    double r = std::hypot(x, y);
    double theta = std::atan2(y, x);
    return {r, theta};
}

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

Path calc_optimal_path(double sx, double sy, double syaw, double gx, double gy, double gyaw, double maxc, double step_size = STEP_SIZE) {
    std::vector<Path> paths = calc_all_paths(sx, sy, syaw, gx, gy, gyaw, maxc, step_size);
    
    double minL = paths[0].L;
    int mini = 0;
    
    for (int i = 0; i < paths.size(); ++i) {
        if (paths[i].L <= minL) {
            minL = paths[i].L;
            mini = i;
        }
    }
    
    return paths[mini];
}

// Updated function definition without default arguments
std::vector<Path> calc_all_paths(double sx, double sy, double syaw, double gx, double gy, double gyaw, double maxc, double step_size) {
    std::vector<double> q0 = {sx, sy, syaw};
    std::vector<double> q1 = {gx, gy, gyaw};
    
    std::vector<Path> paths = generate_path(q0, q1, maxc);
    
    for (auto& path : paths) {
        auto [px, py] = generate_local_course(path.L, path.lengths, path.ctypes, maxc, step_size * maxc);
        std::vector<double> yaw;
        std::vector<double> directions;  // Assume this is originally a double vector

        for (size_t i = 0; i < px.size(); ++i) {
            double ix = px[i];
            double iy = py[i];
            if (i >= path.yaw.size()) {
                std::cerr << "Error: Out-of-bounds access for path.yaw at index " << i << std::endl;
                break;
            }
            double iyaw = path.yaw[i];
            
            px[i] = cos(-q0[2]) * ix + sin(-q0[2]) * iy + q0[0];
            py[i] = -sin(-q0[2]) * ix + cos(-q0[2]) * iy + q0[1];
            yaw.push_back(pi_2_pi(iyaw + q0[2]));
        }

        std::vector<int> int_directions(directions.begin(), directions.end());
        
        path.x = px;
        path.y = py;
        path.yaw = yaw;
        path.directions = int_directions;  // Now assigning correctly typed vector

        for (auto& l : path.lengths) {
            l /= maxc;
        }
        
        path.L /= maxc;
    }
    
    return paths;
}


std::vector<Path> set_path(std::vector<Path> paths, const std::vector<double>& lengths, const std::vector<std::string>& ctypes) {
    Path path({}, {}, 0.0, {}, {}, {}, {});
    path.ctypes = ctypes;
    path.lengths = lengths;

    for (auto& path_e : paths) {
        if (path_e.ctypes == path.ctypes) {
            double sum_diff = 0;
            for (size_t i = 0; i < path_e.lengths.size(); ++i) {
                sum_diff += std::abs(path_e.lengths[i] - path.lengths[i]);
            }
            if (sum_diff <= 0.01) return paths;
        }
    }

    path.L = 0;
    for (const auto& l : lengths) {
        path.L += std::abs(l);
    }

    if (path.L >= MAX_LENGTH) return paths;
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
    auto [u1, t1] = R(x + std::sin(phi), y - 1.0 - std::cos(phi));
    u1 = u1 * u1;

    if (u1 >= 4.0) {
        double u = std::sqrt(u1 - 4.0);
        double theta = std::atan2(2.0, u);
        double t = M(t1 + theta);
        double v = M(t - phi);

        if (t >= 0.0 && v >= 0.0) {
            return {true, t, u, v};
        }
    }

    return {false, 0.0, 0.0, 0.0};
}

std::tuple<bool, double, double, double> LRL(double x, double y, double phi) {
    auto [u1, t1] = R(x - std::sin(phi), y - 1.0 + std::cos(phi));

    if (u1 <= 4.0) {
        double u = -2.0 * std::asin(0.25 * u1);
        double t = M(t1 + 0.5 * u + PI);
        double v = M(phi - t + u);

        if (t >= 0.0 && u <= 0.0) {
            return {true, t, u, v};
        }
    }

    return {false, 0.0, 0.0, 0.0};
}

std::vector<Path> SCS(double x, double y, double phi, std::vector<Path>& paths) {
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
        double v = std::sqrt(std::pow(x - xd, 2) + std::pow(y, 2)) - std::tan(phi / 2.0);
        return {true, t, u, v};
    } else if (y < 0.0 && 0.0 < phi && phi < PI * 0.99) {
        double xd = -y / std::tan(phi) + x;
        double t = xd - std::tan(phi / 2.0);
        double u = phi;
        double v = -std::sqrt(std::pow(x - xd, 2) + std::pow(y, 2)) - std::tan(phi / 2.0);
        return {true, t, u, v};
    }

    return {false, 0.0, 0.0, 0.0};
}

std::vector<Path> CSC(double x, double y, double phi, std::vector<Path>& paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LSL(x, y, phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"L", "S", "L"});

    std::tie(flag, t, u, v) = LSL(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, -u, -v}, {"L", "S", "L"});

    std::tie(flag, t, u, v) = LSL(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"R", "S", "R"});

    std::tie(flag, t, u, v) = LSL(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, -u, -v}, {"R", "S", "R"});

    std::tie(flag, t, u, v) = LSR(x, y, phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"L", "S", "R"});

    std::tie(flag, t, u, v) = LSR(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, -u, -v}, {"L", "S", "R"});

    std::tie(flag, t, u, v) = LSR(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"R", "S", "L"});

    std::tie(flag, t, u, v) = LSR(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, -u, -v}, {"R", "S", "L"});

    return paths;
}

std::vector<Path> CCC(double x, double y, double phi, std::vector<Path>& paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRL(x, y, phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"L", "R", "L"});

    std::tie(flag, t, u, v) = LRL(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, -u, -v}, {"L", "R", "L"});

    std::tie(flag, t, u, v) = LRL(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, u, v}, {"R", "L", "R"});

    std::tie(flag, t, u, v) = LRL(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, -u, -v}, {"R", "L", "R"});

    double xb = x * std::cos(phi) + y * std::sin(phi);
    double yb = x * std::sin(phi) - y * std::cos(phi);

    std::tie(flag, t, u, v) = LRL(xb, yb, phi);
    if (flag) paths = set_path(paths, {v, u, t}, {"L", "R", "L"});

    std::tie(flag, t, u, v) = LRL(-xb, yb, -phi);
    if (flag) paths = set_path(paths, {-v, -u, -t}, {"L", "R", "L"});

    std::tie(flag, t, u, v) = LRL(xb, -yb, -phi);
    if (flag) paths = set_path(paths, {v, u, t}, {"R", "L", "R"});

    std::tie(flag, t, u, v) = LRL(-xb, -yb, phi);
    if (flag) paths = set_path(paths, {-v, -u, -t}, {"R", "L", "R"});

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

    return {tau, omega};
}

std::tuple<bool, double, double, double> LRLRn(double x, double y, double phi) {
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    double rho = 0.25 * (2.0 + std::sqrt(xi * xi + eta * eta));

    if (rho <= 1.0) {
        double u = std::acos(rho);
        auto [t, v] = calc_tauOmega(u, -u, xi, eta, phi);
        if (t >= 0.0 && v <= 0.0) {
            return {true, t, u, v};
        }
    }

    return {false, 0.0, 0.0, 0.0};
}

std::tuple<bool, double, double, double> LRLRp(double x, double y, double phi) {
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    double rho = (20.0 - xi * xi - eta * eta) / 16.0;

    if (0.0 <= rho && rho <= 1.0) {
        double u = -std::acos(rho);
        if (u >= -0.5 * PI) {
            auto [t, v] = calc_tauOmega(u, u, xi, eta, phi);
            if (t >= 0.0 && v >= 0.0) {
                return {true, t, u, v};
            }
        }
    }

    return {false, 0.0, 0.0, 0.0};
}

std::vector<Path> CCCC(double x, double y, double phi, std::vector<Path>& paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRLRn(x, y, phi);
    if (flag) paths = set_path(paths, {t, u, -u, v}, {"L", "R", "L", "R"});

    std::tie(flag, t, u, v) = LRLRn(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, -u, u, -v}, {"L", "R", "L", "R"});

    std::tie(flag, t, u, v) = LRLRn(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, u, -u, v}, {"R", "L", "R", "L"});

    std::tie(flag, t, u, v) = LRLRn(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, -u, u, -v}, {"R", "L", "R", "L"});

    std::tie(flag, t, u, v) = LRLRp(x, y, phi);
    if (flag) paths = set_path(paths, {t, u, u, v}, {"L", "R", "L", "R"});

    std::tie(flag, t, u, v) = LRLRp(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, -u, -u, -v}, {"L", "R", "L", "R"});

    std::tie(flag, t, u, v) = LRLRp(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, u, u, v}, {"R", "L", "R", "L"});

    std::tie(flag, t, u, v) = LRLRp(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, -u, -u, -v}, {"R", "L", "R", "L"});

    return paths;
}

std::tuple<bool, double, double, double> LRSR(double x, double y, double phi) {
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    auto [rho, theta] = R(-eta, xi);

    if (rho >= 2.0) {
        double t = theta;
        double u = 2.0 - rho;
        double v = M(t + 0.5 * PI - phi);
        if (t >= 0.0 && u <= 0.0 && v <= 0.0) {
            return {true, t, u, v};
        }
    }

    return {false, 0.0, 0.0, 0.0};
}

std::tuple<bool, double, double, double> LRSL(double x, double y, double phi) {
    double xi = x - std::sin(phi);
    double eta = y - 1.0 + std::cos(phi);
    auto [rho, theta] = R(xi, eta);

    if (rho >= 2.0) {
        double r = std::sqrt(rho * rho - 4.0);
        double u = 2.0 - r;
        double t = M(theta + std::atan2(r, -2.0));
        double v = M(phi - 0.5 * PI - t);
        if (t >= 0.0 && u <= 0.0 && v <= 0.0) {
            return {true, t, u, v};
        }
    }

    return {false, 0.0, 0.0, 0.0};
}

std::vector<Path> CCSC(double x, double y, double phi, std::vector<Path>& paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRSL(x, y, phi);
    if (flag) paths = set_path(paths, {t, -0.5 * PI, u, v}, {"L", "R", "S", "L"});

    std::tie(flag, t, u, v) = LRSL(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"L", "R", "S", "L"});

    std::tie(flag, t, u, v) = LRSL(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, -0.5 * PI, u, v}, {"R", "L", "S", "R"});

    std::tie(flag, t, u, v) = LRSL(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"R", "L", "S", "R"});

    std::tie(flag, t, u, v) = LRSR(x, y, phi);
    if (flag) paths = set_path(paths, {t, -0.5 * PI, u, v}, {"L", "R", "S", "R"});

    std::tie(flag, t, u, v) = LRSR(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"L", "R", "S", "R"});

    std::tie(flag, t, u, v) = LRSR(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, -0.5 * PI, u, v}, {"R", "L", "S", "L"});

    std::tie(flag, t, u, v) = LRSR(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, 0.5 * PI, -u, -v}, {"R", "L", "S", "L"});

    double xb = x * std::cos(phi) + y * std::sin(phi);
    double yb = x * std::sin(phi) - y * std::cos(phi);

    std::tie(flag, t, u, v) = LRSL(xb, yb, phi);
    if (flag) paths = set_path(paths, {v, u, -0.5 * PI, t}, {"L", "S", "R", "L"});

    std::tie(flag, t, u, v) = LRSL(-xb, yb, -phi);
    if (flag) paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"L", "S", "R", "L"});

    std::tie(flag, t, u, v) = LRSL(xb, -yb, -phi);
    if (flag) paths = set_path(paths, {v, u, -0.5 * PI, t}, {"R", "S", "L", "R"});

    std::tie(flag, t, u, v) = LRSL(-xb, -yb, phi);
    if (flag) paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"R", "S", "L", "R"});

    std::tie(flag, t, u, v) = LRSR(xb, yb, phi);
    if (flag) paths = set_path(paths, {v, u, -0.5 * PI, t}, {"R", "S", "R", "L"});

    std::tie(flag, t, u, v) = LRSR(-xb, yb, -phi);
    if (flag) paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"R", "S", "R", "L"});

    std::tie(flag, t, u, v) = LRSR(xb, -yb, -phi);
    if (flag) paths = set_path(paths, {v, u, -0.5 * PI, t}, {"L", "S", "L", "R"});

    std::tie(flag, t, u, v) = LRSR(-xb, -yb, phi);
    if (flag) paths = set_path(paths, {-v, -u, 0.5 * PI, -t}, {"L", "S", "L", "R"});

    return paths;
}

std::tuple<bool, double, double, double> LRSLR(double x, double y, double phi) {
    double xi = x + std::sin(phi);
    double eta = y - 1.0 - std::cos(phi);
    auto [rho, theta] = R(xi, eta);

    if (rho >= 2.0) {
        double u = 4.0 - std::sqrt(rho * rho - 4.0);
        if (u <= 0.0) {
            double t = M(std::atan2((4.0 - u) * xi - 2.0 * eta, -2.0 * xi + (u - 4.0) * eta));
            double v = M(t - phi);

            if (t >= 0.0 && v >= 0.0) {
                return {true, t, u, v};
            }
        }
    }

    return {false, 0.0, 0.0, 0.0};
}

std::vector<Path> CCSCC(double x, double y, double phi, std::vector<Path>& paths) {
    bool flag;
    double t, u, v;

    std::tie(flag, t, u, v) = LRSLR(x, y, phi);
    if (flag) paths = set_path(paths, {t, -0.5 * PI, u, -0.5 * PI, v}, {"L", "R", "S", "L", "R"});

    std::tie(flag, t, u, v) = LRSLR(-x, y, -phi);
    if (flag) paths = set_path(paths, {-t, 0.5 * PI, -u, 0.5 * PI, -v}, {"L", "R", "S", "L", "R"});

    std::tie(flag, t, u, v) = LRSLR(x, -y, -phi);
    if (flag) paths = set_path(paths, {t, -0.5 * PI, u, -0.5 * PI, v}, {"R", "L", "S", "R", "L"});

    std::tie(flag, t, u, v) = LRSLR(-x, -y, phi);
    if (flag) paths = set_path(paths, {-t, 0.5 * PI, -u, 0.5 * PI, -v}, {"R", "L", "S", "R", "L"});

    return paths;
}

std::pair<std::vector<double>, std::vector<double>> generate_local_course(double L, const std::vector<double>& lengths, 
                                                                          const std::vector<std::string>& mode, double maxc, double step_size) {
    int point_num = static_cast<int>(L / step_size) + static_cast<int>(lengths.size()) + 3;

    std::vector<double> px(point_num, 0.0);
    std::vector<double> py(point_num, 0.0);
    std::vector<double> pyaw(point_num, 0.0);
    std::vector<int> directions(point_num, 0);
    int ind = 1;

   // Add debug output to verify the size and content of yaw vector
    std::cout << "Generating local course. Initial yaw size: " << pyaw.size() << std::endl;


    directions[0] = lengths[0] > 0.0 ? 1 : -1;
    double d = directions[0] > 0.0 ? step_size : -step_size;
    double pd = d;
    double ll = 0.0;

    for (size_t i = 0; i < mode.size(); ++i) {
        const std::string& m = mode[i];
        double l = lengths[i];

        if (l > 0.0) {
            d = step_size;
        } else {
            d = -step_size;
        }

        double ox = px[ind];
        double oy = py[ind];
        double oyaw = pyaw[ind];

        ind -= 1;
        if (i >= 1 && lengths[i - 1] * lengths[i] > 0) {
            pd = -d - ll;
        } else {
            pd = d - ll;
        }

        while (fabs(pd) <= fabs(l)) {
            ++ind;
            interpolate(ind, pd, m, maxc, ox, oy, oyaw, px, py, pyaw, directions);
            pd += d;
        }

        ll = l - pd - d;

        ++ind;
        interpolate(ind, l, m, maxc, ox, oy, oyaw, px, py, pyaw, directions);
    }

    while (px.back() == 0.0 && py.back() == 0.0) {
        px.pop_back();
        py.pop_back();
        pyaw.pop_back();
        directions.pop_back();
    }

    if (pyaw.empty()) {
        std::cerr << "Error: yaw vector is empty after generation." << std::endl;
    } else {
        std::cout << "Generated yaw size: " << pyaw.size() << std::endl;
    }

    return {px, py};

}

void interpolate(int ind, double l, const std::string& m, double maxc, double ox, double oy, double oyaw, 
                 std::vector<double>& px, std::vector<double>& py, std::vector<double>& pyaw, std::vector<int>& directions) {
    if (m == "S") {
        px[ind] = ox + l / maxc * std::cos(oyaw);
        py[ind] = oy + l / maxc * std::sin(oyaw);
        pyaw[ind] = oyaw;
    } else {
        double ldx = std::sin(l) / maxc;
        double ldy;
        if (m == "L") {
            ldy = (1.0 - std::cos(l)) / maxc;
        } else if (m == "R") {
            ldy = (1.0 - std::cos(l)) / (-maxc);
        }

        double gdx = std::cos(-oyaw) * ldx + std::sin(-oyaw) * ldy;
        double gdy = -std::sin(-oyaw) * ldx + std::cos(-oyaw) * ldy;
        px[ind] = ox + gdx;
        py[ind] = oy + gdy;

        if (m == "L") {
            pyaw[ind] = oyaw + l;
        } else if (m == "R") {
            pyaw[ind] = oyaw - l;
        }
    }

    if (l > 0.0) {
        directions[ind] = 1;
    } else {
        directions[ind] = -1;
    }
}

std::vector<Path> generate_path(const std::vector<double>& q0, const std::vector<double>& q1, double maxc) {
    double dx = q1[0] - q0[0];
    double dy = q1[1] - q0[1];
    double dth = q1[2] - q0[2];
    double c = std::cos(q0[2]);
    double s = std::sin(q0[2]);
    double x = (c * dx + s * dy) * maxc;
    double y = (-s * dx + c * dy) * maxc;

    std::vector<Path> paths;
    paths = SCS(x, y, dth, paths);
    paths = CSC(x, y, dth, paths);
    paths = CCC(x, y, dth, paths);
    paths = CCCC(x, y, dth, paths);
    paths = CCSC(x, y, dth, paths);
    paths = CCSCC(x, y, dth, paths);

    return paths;
}

// The main simulation function
void main_simulation() {
    std::vector<std::tuple<double, double, double>> states = {
        {-3, 3, 120}, {10, -7, 30}, {10, 13, 30}, {20, 5, -25},
        {35, 10, 180}, {32, -10, 180}, {5, -12, 90}
    };

    double max_c = 0.1;
    std::vector<double> path_x, path_y, yaw;

    for (size_t i = 0; i < states.size() - 1; ++i) {
        auto [s_x, s_y, s_yaw_deg] = states[i];
        auto [g_x, g_y, g_yaw_deg] = states[i + 1];

        double s_yaw = M_PI / 180.0 * s_yaw_deg;
        double g_yaw = M_PI / 180.0 * g_yaw_deg;

        std::cout << "Calculating path from state " << i << " to " << i + 1 << std::endl;

        Path path_i = calc_optimal_path(s_x, s_y, s_yaw, g_x, g_y, g_yaw, max_c);

        std::cout << "Generated path sizes: x=" << path_i.x.size() << ", y=" << path_i.y.size() << ", yaw=" << path_i.yaw.size() << std::endl;

        if (path_i.x.empty() || path_i.y.empty() || path_i.yaw.empty()) {
            std::cerr << "Error: Path segment is empty. Check path calculation." << std::endl;
            return;
        }

        path_x.insert(path_x.end(), path_i.x.begin(), path_i.x.end());
        path_y.insert(path_y.end(), path_i.y.begin(), path_i.y.end());
        yaw.insert(yaw.end(), path_i.yaw.begin(), path_i.yaw.end());
    }

    plt::ion();
    plt::figure(1);

    for (size_t i = 0; i < path_x.size(); ++i) {
        std::cout << "Plotting point " << i << ": (" << path_x[i] << ", " << path_y[i] << ", " << yaw[i] << ")" << std::endl;
        
        if (i >= yaw.size()) {
            std::cerr << "Error: yaw index out-of-bounds at index " << i << " with yaw size " << yaw.size() << std::endl;
            break;  // Prevent further out-of-bounds access
        }

        plt::clf();
        plt::plot(path_x, path_y, {{"linewidth", "1"}, {"color", "gray"}});

        for (const auto& [x, y, theta_deg] : states) {
            Arrow(x, y, M_PI / 180.0 * theta_deg, 2, "blueviolet");
        }

        // Uncomment to visualize the car if needed:
        // Car(path_x[i], path_y[i], yaw[i], 1.5, 3);

        plt::axis("equal");
        plt::title("Simulation of Reeds-Shepp Curves");
        plt::axis({-10, 42, -20, 20});
        plt::draw();
        plt::pause(0.001);
    }

    plt::pause(1);
}

int main() {
    main_simulation();
    return 0;
}
