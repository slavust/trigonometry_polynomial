#pragma once

#include "helper_math.h"
#include "polynomial.h"

#include <cstddef>
#include <cmath>

namespace _math_detail
{
    template<typename computation_type, std::size_t n>
    constexpr polynomial<computation_type, n> create_cos_polynomial()
    {
        constexpr computation_type start = 0;
        constexpr computation_type end = static_cast<computation_type>(M_PI*2);
        constexpr computation_type step = (end - start) / (n-1);

        typename polynomial<computation_type, n>::table_type table;
        for(std::size_t i = 0; i < n; ++i)
        {
            const computation_type cur_param = start + step * i;
            const computation_type cur_value = _math_helper::cos_tailor<computation_type, 10>(cur_param);
            table[i] = std::make_pair(cur_param, cur_value);
        }

        return polynomial<computation_type, n>::interpolate(table);
    }

    template<typename computation_type, std::size_t n>
    constexpr polynomial<computation_type, n> create_asin_polynomial()
    {
        constexpr computation_type start = -1;
        constexpr computation_type end = 1;
        constexpr computation_type step = (end - start) / (n-1);

        typename polynomial<computation_type, n>::table_type table;
        for(std::size_t i = 0; i < n; ++i)
        {
            const computation_type cur_param = start + step * i;
            const computation_type cur_value = _math_helper::asin_tailor<computation_type, 10>(cur_param);
            table[i] = std::make_pair(cur_param, cur_value);
        }

        return polynomial<computation_type, n>::interpolate(table);
    }
}



template<typename computation_type, std::size_t polynomial_size>
struct trigonometry
{
    constexpr static computation_type cos(computation_type x)
    {
        constexpr auto polynomial_interpolation = _math_detail::create_cos_polynomial<computation_type, polynomial_size>();
        constexpr computation_type min = 0;
        constexpr computation_type max = static_cast<computation_type>(M_PI*2);
        constexpr computation_type period = static_cast<computation_type>(M_PI*2);
        if (x < min)
            x -= (static_cast<int>(x / period) -1) * period;
        else if (x > max)
            x -= static_cast<int>(x / period) * period;
        return polynomial_interpolation(x);
    }

    constexpr static computation_type sin(computation_type x)
    {
        constexpr computation_type theta = M_PI_2;
        return cos(theta - x);
    }

    constexpr static computation_type tan(computation_type x)
    {
        const auto _sin = sin(x);
        const auto _cos = cos(x);
        return _sin / _cos;
    }

    constexpr static computation_type asin(computation_type x)
    {
        constexpr auto polynomial_interpolation = _math_detail::create_asin_polynomial<computation_type, polynomial_size>();
        return polynomial_interpolation(x);
    }

    constexpr static computation_type acos(computation_type x)
    {
        constexpr computation_type pi_two = M_PI_2;
        return pi_two - asin(x);
    }
    
    //constexpr static computation_type atan(computation_type x)
    //{
    //    return asin(x / _math_helper::sqrt<computation_type>(1 + x*x));
    //}

    constexpr static computation_type atan2(computation_type rsin, computation_type rcos)
    {
        const computation_type r = _math_helper::sqrt<computation_type>(rcos*rcos + rsin*rsin);
        const computation_type _asin = asin(rsin / r);

        if(rcos >= 0)
            return _asin;
        
        if(rsin > 0)
            return static_cast<computation_type>(M_PI) - _asin;
        return -static_cast<computation_type>(M_PI) - _asin;
    }
};