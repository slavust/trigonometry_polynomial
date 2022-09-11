#include "polynomial.h"
#include "trigonometry.h"
#include <cmath>
#include <functional>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <stdlib.h>

template<typename real>
using unary_function = std::function<real(real)>;

template<typename real>
void compare_functions_output(unary_function<real> i_func1, 
                            unary_function<real> i_func2, 
                            const real param_min, 
                            const real param_max,
                            const std::string& i_name,
                            const std::string& i_suffix1,
                            const std::string& i_suffix2)
{
    const std::string points_filename = i_name + "_points.dat";
    const std::string func1_name = i_name + " " + i_suffix1;
    const std::string func2_name = i_name + " " + i_suffix2;

    { // generate points
        constexpr std::size_t num_points_for_plot = 200;
        std::vector<real> param(num_points_for_plot, 0), res1(num_points_for_plot, 0), res2(num_points_for_plot, 0);

        const real step = (param_max - param_min) / (num_points_for_plot - 1);
        for(std::size_t i = 0; i < num_points_for_plot; ++i)
        {
            const real cur_param = param_min + step*i;
            param[i] = cur_param;
            res1[i] = i_func1(cur_param);
            res2[i] = i_func2(cur_param);
        }

        std::ofstream points(points_filename.c_str());
        points<<"#\t"<<"param"<<"\t\t"<<func1_name<<"\t\t"<<func2_name<<std::endl;
        for(std::size_t i = 0; i < num_points_for_plot; ++i)
        {
            points<<"\t"<<param[i]<<"\t\t"<<res1[i]<<"\t\t"<<res2[i]<<std::endl;
        }
    }
    { // generate gnuplot script
        std::ofstream plot(i_name + ".p");
        plot<<"set title \""<<i_name<<"\""<<std::endl;
        plot<<"plot \""<<points_filename<<"\" using 1:2 title \""<<func1_name<<"\" with lines lw 3 lc rgb \"#800000FF\", \"\" using 1:3 title \""<<func2_name<<"\" with lines lw 3 lt 3 lc rgb \"#8000FF00\""<<std::endl; 
    }
}


template <typename real>
void compare_functions_performance(std::vector<std::pair<unary_function<real>, unary_function<real>>> i_functions, 
    std::vector<std::pair<real, real>> i_domains, 
    std::vector<std::string> i_names,
    std::pair<std::string, std::string> i_suffixes)
{
    constexpr std::size_t num_calls = 10000000;

    std::vector<std::pair<std::size_t, std::size_t>> performance(i_functions.size(), std::make_pair(0, 0));
    for(std::size_t i = 0; i < i_functions.size(); ++i)
    {
        const auto& [func1, func2] = i_functions[i];
        const auto [domain_min, domain_max] = i_domains[i];
        const auto step = (domain_max - domain_min) / (num_calls-1);

        using std::chrono::steady_clock;

        const auto func1_start = steady_clock::now();
        for(std::size_t c = 0; c < num_calls; ++c)
        {
            const real param = domain_min + step*c;
            const real result = func1(param);
            static_cast<void>(result);
        }
        const auto func2_start = steady_clock::now();
        for(std::size_t c = 0; c < num_calls; ++c)
        {
            const real param = domain_min + step*c;
            const real result = func2(param);
            static_cast<void>(result);
        }
        const auto func2_end = steady_clock::now();

        const auto func1_ms = std::chrono::duration_cast<std::chrono::milliseconds>(func2_start - func1_start);
        const auto func2_ms = std::chrono::duration_cast<std::chrono::milliseconds>(func2_end - func2_start);
        performance[i] = std::make_pair(func1_ms.count(), func2_ms.count());
    }
    { // performance data
        std::ofstream performance_dat("performance.dat");
        for(std::size_t i = 0; i < performance.size(); ++i)
        {
            const std::string& name = i_names[i];
            const auto [func1_ms, func2_ms] = performance[i];
            performance_dat<<name<<"\t"<<static_cast<double>(func1_ms) / num_calls<<"\t"<<static_cast<double>(func2_ms) / num_calls<<std::endl;
        }
    }
    { // performance plot
        std::ofstream performance_plot("performance.p");
        performance_plot<<"set style data histogram"<<std::endl;
        performance_plot<<"set style histogram cluster gap 1"<<std::endl;
        performance_plot<<"set style fill solid"<<std::endl;
        performance_plot<<"set boxwidth 0.9"<<std::endl;
        performance_plot<<"set xtics format \"\""<<std::endl;
        performance_plot<<"set grid ytics"<<std::endl;
        performance_plot<<"set title \"performance comparison\""<<std::endl;
        performance_plot<<"set xlabel \"function\""<<std::endl;
        performance_plot<<"set ylabel \"avg time (ms)\""<<std::endl;
        performance_plot<<"plot \"performance.dat\" using 2:xtic(1) title \""<<i_suffixes.first<<"\", \\"<<std::endl;
        performance_plot<<"\t\"performance.dat\" using 3:xtic(1) title \""<<i_suffixes.second<<"\""<<std::endl;
        performance_plot<<std::endl;
    }
}

int main()
{
    using real = float;
    constexpr std::size_t lookup_table_size = 5;

    using interpolation = trigonometry<real, lookup_table_size>;

    std::vector<std::pair<unary_function<real>, unary_function<real>>> unary_functions;
    std::vector<std::pair<real, real>> domains;
    std::vector<std::string> names;
    
    unary_functions.emplace_back(std::make_pair([](real x){return std::cos(x);}, [](real x){return interpolation::cos(x);}));
    domains.emplace_back(std::make_pair(-M_PI*4, M_PI*4));
    names.emplace_back("cos");

    unary_functions.emplace_back(std::make_pair([](real x){return std::sin(x);}, [](real x){return interpolation::sin(x);}));
    domains.emplace_back(std::make_pair(-M_PI*4, M_PI*4));
    names.emplace_back("sin");

    unary_functions.emplace_back(std::make_pair([](real x){return std::tan(x);}, [](real x){return interpolation::tan(x);}));
    domains.emplace_back(std::make_pair(-M_PI*4, M_PI*4));
    names.emplace_back("tan");

    unary_functions.emplace_back(std::make_pair([](real x){return std::acos(x);}, [](real x){return interpolation::acos(x);}));
    domains.emplace_back(std::make_pair(-1, 1));
    names.emplace_back("acos");

    unary_functions.emplace_back(std::make_pair([](real x){return std::asin(x);}, [](real x){return interpolation::asin(x);}));
    domains.emplace_back(std::make_pair(-1, 1));
    names.emplace_back("asin");
    
    //unary_functions.emplace_back(std::make_pair([](real x){return std::atan(x);}, [](real x){return interpolation::atan(x);}));
    //domains.emplace_back(std::make_pair(-5, 5));
    //names.emplace_back("atan");
    
    unary_functions.emplace_back(std::make_pair([](real x){return std::atan2(x, std::sqrt(1 - x*x));}, [](real x){return interpolation::atan2(x, std::sqrt(1 - x*x));}));
    domains.emplace_back(std::make_pair(-1, 1));
    names.emplace_back("atan2_1");

    unary_functions.emplace_back(std::make_pair([](real x){return std::atan2(x, -std::sqrt(1 - x*x));}, [](real x){return interpolation::atan2(x, -std::sqrt(1 - x*x));}));
    domains.emplace_back(std::make_pair(-1, 1));
    names.emplace_back("atan2_2");
    
    for(std::size_t i = 0; i < unary_functions.size(); ++i)
    {
        compare_functions_output(unary_functions[i].first, 
            unary_functions[i].second, 
            domains[i].first, 
            domains[i].second, 
            names[i], 
            "std",
            "poly");
    }
    compare_functions_performance(unary_functions, 
                                 domains,
                                 names,
                                 std::make_pair("std", "poly"));

    for(const auto& name : names)
    {
        const std::string command = "gnuplot -p \"" + name + "\".p";
        system(command.c_str());
    }

    const std::string performance_command = "gnuplot -p \"performance.p\"";
    system(performance_command.c_str());

    return 0;
}