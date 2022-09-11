#pragma once

#include "helper_math.h"

#include <cstddef>
#include <array>


template<typename computation_type, std::size_t size>
class polynomial_helper {
public:
    typedef std::array<computation_type, size> vector;
    typedef std::array<vector, size> matrix;


    constexpr static matrix create_Vandermonde_U_inv(const vector &samples_x) {
        matrix U_inv;

        // diagonal
        for (std::size_t i = 0; i < size; i++)
            U_inv[i][i] = static_cast<computation_type>(1);

        // first column
        for (std::size_t i = 1; i < size; i++)
            U_inv[i][0] = static_cast<computation_type>(0);

        // under diagonal
        for (std::size_t i = 2; i < size; i++) {
            for (std::size_t j = 1; j < i; j++) {
                U_inv[i][j] = static_cast<computation_type>(0);
            }
        }

        // first row
        for (std::size_t j = 1; j < size; j++) {
            const computation_type left = U_inv[0][j - 1];
            U_inv[0][j] = -left * samples_x[j - 1];
        }

        // above diagonal, from second row
        for (std::size_t i = 1; i < size - 1; i++) {
            for (std::size_t j = i + 1; j < size; j++) {
                const computation_type up_left = U_inv[i - 1][j - 1];
                const computation_type left = U_inv[i][j - 1];
                U_inv[i][j] = up_left - left * samples_x[j - 1];
            }
        }

        return U_inv;
    }

    constexpr static matrix create_Vandermonde_L_inv(const vector &samples_x) {
        matrix L_inv;

        L_inv[0][0] = static_cast<computation_type>(1);

        // above diagonal
        for (std::size_t i = 0; i < size - 1; i++)
            for (std::size_t j = i + 1; j < size; j++)
                L_inv[i][j] = static_cast<computation_type>(0);

        // diagonal
        for (std::size_t ij = 1; ij < size; ij++) {
            computation_type product = 1.0f;
            for (std::size_t k = 0; k < ij; k++) {
                product *= static_cast<computation_type>(1) /
                            (samples_x[ij] - samples_x[k]);
            }
            L_inv[ij][ij] = product;
        }

        // below diagonal
        for (std::size_t i = 1; i < size; i++) {
            for (std::size_t j = 0; j < i; j++) {
                L_inv[i][j] = L_inv[i - 1][j] * (static_cast<computation_type>(1) / 
                            (samples_x[j] - samples_x[i]));
            }
        }
        return L_inv;
    }

    constexpr static vector mul(const matrix &m, const vector &v) {
        vector r;

        for (std::size_t i = 0; i < r.size(); i++) {
            r[i] = static_cast<computation_type>(0);

            for (int k = 0; k < r.size(); k++) {
                r[i] += m[i][k] * v[k];
            }
        }

        return r;
    }
};

namespace polynomial_detail
{
    template<typename computation_type, std::size_t size, std::size_t idx>
    struct evaluate_polynomial
    {
        constexpr static computation_type value(const typename polynomial_helper<computation_type, size>::vector& coefficients, computation_type param)
        {
            return evaluate_polynomial<computation_type, size, idx-1>::value(coefficients, param) 
                + _math_helper::pow<computation_type, idx>(param) * coefficients[idx];
        }
    };

    template<typename computation_type, std::size_t size>
    struct evaluate_polynomial<computation_type, size, 0>
    {
        constexpr static computation_type value(const typename polynomial_helper<computation_type, size>::vector& coefficients, computation_type param)
        {
            return coefficients[0];
        }
    };
}

/// \brief Polynomial template class
/// coefficients and evaluation precision are set by template parameter 'type'
template<typename type, std::size_t size>
class polynomial {
public:
    template <typename table_type>
    using table_record = std::pair<table_type, table_type>;

    template <typename table_type, std::size_t table_size> 
    using table = std::array<table_record<table_type>, table_size>;

    using table_record_type = table_record<type>;
    using table_type = table<type, size>;
    using coefficient_list_type = std::array<type, size>;

    /// \brief power of polynomial
    /// power of polynomial is max power of variable in polynomial
    constexpr static std::size_t power = size-1;

protected:
    coefficient_list_type mCoefficients;

public:
    polynomial() = delete;

    /// \brief initializes polynomial with given coefficients, min and max domains
    /// power is set explicitly
    constexpr polynomial(const coefficient_list_type &coefficients) : 
        mCoefficients(coefficients) {
    }

    /// \brief copy constructor
    constexpr polynomial(const polynomial& p) :
        mCoefficients(p.mCoefficients)
    {
    }
    
    /// \brief returns coefficient values
    /// 0 is free coefficient:
    /// c0 + c1*x1 + c2*x2^2 + ... + cn*xn^n
    constexpr const coefficient_list_type& getCoefficients() const { return mCoefficients; }

    /// \brief calculates polynomial value in given x
    ///
    constexpr type evaluate(type x) const {
        return polynomial_detail::evaluate_polynomial<type, size, size-1>::value(mCoefficients, x);
    }

    /// \brief calculates polynomial value in given x
    ///
    constexpr type operator() (type x) const
    {
        return evaluate(x);
    }

    /// \brief Interpolates table function
    /// \param points: map of function parameters and values
    constexpr static polynomial<type, size> interpolate(const table_type& sample_table) {
        static_assert(size > 0, "cannot interpolate using empty table");
        constexpr std::size_t sample_count = size;
        constexpr std::size_t coefficient_count = size;
        constexpr std::size_t polynomial_power = size - 1;

        typename polynomial_helper<type, sample_count>::vector X;
        typename polynomial_helper<type, sample_count>::vector Y;

        std::size_t indx = 0;
        for (const auto [param, value] : sample_table) {
            X[indx] = param;
            Y[indx] = value;
            indx++;
        }

        const auto L_inv = polynomial_helper<type, size>::create_Vandermonde_L_inv(X);

        const auto U_inv = polynomial_helper<type, size>::create_Vandermonde_U_inv(X);

        const auto C = polynomial_helper<type, size>::mul(
            U_inv, 
            polynomial_helper<type, size>::mul(L_inv, Y));

        return polynomial<type, size>(C);
    }
};
