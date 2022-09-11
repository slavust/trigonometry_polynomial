#pragma once

#include <cstddef>
#include <limits>

namespace _math_helper
{
    // compile-time pow
    template <typename computation_type, std::size_t power>
    constexpr computation_type pow(computation_type base)
    {
        return base * pow<computation_type, power-1>(base);
    }
    template<>
    constexpr double pow<double, 0>(double base)
    {
        return 1.0;
    }
    template<>
    constexpr float pow<float, 0>(float base)
    {
        return 1.0f;
    }
    template<>
    constexpr size_t pow<size_t, 0>(size_t base)
    {
        return 1;
    }

    template <typename computation_type>
    computation_type sqrtNewtonRaphson(computation_type x, computation_type curr, computation_type prev)
	{
		return curr == prev
			? curr
			: sqrtNewtonRaphson(x, static_cast<computation_type>(0.5) * (curr + x / curr), curr);
	}

    //compile-time square root
    template <typename computation_type>
    computation_type sqrt(computation_type x)
    {
    	return x >= 0 && x < std::numeric_limits<computation_type>::infinity()
    		? sqrtNewtonRaphson(x, x, (computation_type)0)
    		: std::numeric_limits<computation_type>::quiet_NaN();
    }

    // compile-time factorial
    template<std::size_t number>
    constexpr std::size_t factorial()
    {
        return number * factorial<number-1>();
    }
    template<>
    constexpr std::size_t factorial<0>()
    {
        return 1;
    }

    //compile-time Tailor series for cosine
    template<typename computation_type, std::size_t n> 
    constexpr computation_type cos_tailor(computation_type x)
    {
        return pow<computation_type, n>(-1) 
            * pow<computation_type, n*2>(x) 
            / static_cast<computation_type>(factorial<n*2>()) 
            + cos_tailor<computation_type, n-1>(x);
    }
    template<>
    constexpr double cos_tailor<double, 0>(double x)
    {
        return 1.0;
    }
    template<>
    constexpr float cos_tailor<float, 0>(float x)
    {
        return 1.0f;
    }

    //compile-time Tailor series for arcsin
    template<typename computation_type, std::size_t n>
    constexpr computation_type asin_tailor(computation_type x)
    {
        return factorial<2*n>() * pow<computation_type, 2*n + 1>(x)
            / (pow<computation_type, n>(4) * pow<std::size_t, 2>(factorial<n>()) * (2*n + 1))
            + asin_tailor<computation_type, n-1>(x);
    }
    template<>
    constexpr double asin_tailor<double, 0>(double x)
    {
        return x;
    }
    template<>
    constexpr float asin_tailor<float, 0>(float x)
    {
        return x;
    }
}