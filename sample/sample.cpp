#include "trigonometry.h"

#include <cmath>
#include <iostream>

int main()
{
    constexpr size_t polynomial_power = 5;
    using floating_point_type = float;

    using math = trigonometry<float, polynomial_power>;

    std::cout<<"sin(pi): "<<math::sin(M_PI)<<std::endl;
    std::cout<<"cos(pi): "<<math::cos(M_PI)<<std::endl;
    std::cout<<"tan(pi): "<<math::tan(M_PI)<<std::endl;

    std::cout<<"asin(0): "<<math::asin(0)<<std::endl;
    std::cout<<"acos(0): "<<math::acos(0)<<std::endl;
    std::cout<<"atan2(1, 0): "<<math::atan2(1,0)<<std::endl;
}