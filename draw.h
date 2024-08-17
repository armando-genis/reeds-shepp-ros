#ifndef DRAW_H
#define DRAW_H

#include <string>

namespace draw {

class Arrow {
public:
    Arrow(double x, double y, double theta, double L, const std::string& c);
private:
    double deg2rad(double deg);
};

class Car {
public:
    Car(double x, double y, double yaw, double w, double L);
};

}  // namespace draw

#endif // DRAW_H
