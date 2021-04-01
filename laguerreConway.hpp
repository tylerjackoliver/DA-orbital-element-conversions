#ifndef __LAGUERRE_CONWAY_H__
#define __LAGUERRE_CONWAY_H__

#include <cmath>
#include <dace/dace_s.h>

/* @brief Implements the Laguerre-Conway method for solving Kepler's equation. Convergence is generally good within 5 iterations.
   @param[in] E0: Initial guess for the eccentric anomaly
   @param[out] Ef: Final iterate; solution of Kepler's equation.
   @param[in] ecc: Eccentricity for Kepler's equation.
   @param[in] M: Mean anomaly for which the solution is desired.
   @param[in]: eps: Tolerance for the solver
   @param[in]: f: Function that returns the error in Kepler's equation
   @param[in]: fp: Function that returns the error in first derivative of Kepler's equation
   @param[in]: fpp: Function that returns the error in second derivative of Kepler's equation
*/
template <typename Function, typename derivFunction>
void laguerreConway(DACE::DA &E0, DACE::DA &Ef, DACE::DA &ecc, DACE::DA &M, double &eps, Function f, derivFunction fp, derivFunction fpp)
{
    const int n = 5;        // Tuning parameter - 5 for now, as per paper
    DACE::DA tolerance = 100.;  // Stopping tolerance
    DACE::DA xi = E0;
    unsigned num_iters = 0;

    while ( abs( DACE::cons(tolerance) ) >= eps && num_iters < 10) // Iteration should be within 4 iterations for majority of E, ecc
    {   
        // Pre-compute function evaluations and derivatives
        DACE::DA fval = f(xi, ecc, M);
        DACE::DA deriv = fp(xi, ecc);
        DACE::DA dDeriv = fpp(xi, ecc);
        DACE::DA numerator = - n * fval;
        DACE::DA root = sqrt( abs((n-1) * (n-1) * (deriv * deriv) - n * (n-1) * fval * dDeriv) );

        // Denominator is such that absolute value is maximised
        DACE::DA denominator = std::max( abs( deriv + root) , abs( deriv - root) );

        // Compute update to iterate
        DACE::DA delta_n1 = numerator / denominator;
        xi += delta_n1; 
        tolerance = delta_n1; // Change in latest iterate
        num_iters++;
    };
    Ef = xi;
}

/* @brief Computes the error in Kepler's equation.
   @param[in] E The value of eccentric anomaly at the current step
   @param[in] ecc The value of eccentricity at the current step
   @param[in] M The value of Mean anomaly at the current step
   @returns The error in Kepler's equation.
*/
DACE::DA keplersEquation(const DACE::DA& E, const DACE::DA& ecc, const DACE::DA& M)
{
    return E - ecc * sin(E) - M;
}

/* @brief Computes the Hyperbolic Kepler's equation.
   @param[in] F The value of eccentric anomaly at the current step
   @param[in] ecc The value of eccentricity at the current step
   @param[in] M The value of Mean anomaly at the current step
   @returns The error in Kepler's equation.    
*/
DACE::DA hypKeplersEquation(const DACE::DA& F, const DACE::DA& ecc, const DACE::DA& M)
{
    return ecc * sinh(F) - ecc - M;
}

/* @brief Computes the error in the first derivative of Kepler's equation.
   @param[in] E The value of eccentric anomaly at the current step
   @param[in] ecc The value of eccentricity at the current step
   @returns The error in the first derivative of Kepler's equation.
*/
DACE::DA dMdE(const DACE::DA& E, const DACE::DA& ecc)
{
    return 1 - ecc * cos(E);
}

/* @brief Computes the error in the first derivative of the hyperbolic kepler's equation.
   @param[in] F The value of hyperbolic eccentric anomaly at the current step
   @param[in] ecc The value of eccentricity at the current step
   @returns The error in the first derivative of Kepler's equation.
*/
DACE::DA hypdMdE(const DACE::DA& F, const DACE::DA& ecc)
{
    return ecc * cosh(F);
}

/* @brief Computes the error in the second derivative of Kepler's equation.
   @param[in] E The value of eccentric anomaly at the current time-step
   @param[in] ecc The value of eccentricity at the current step.
   @returns The error in the second derivative of Kepler's equation.
 */
DACE::DA dMMdEE(const DACE::DA& E, const DACE::DA& ecc)
{
    return ecc * sin(E);
}

/* @brief Computes the error in the second derivative of hyperbolic Kepler's equation.
   @param[in] F The value of eccentric anomaly at the current time-step
   @param[in] ecc The value of eccentricity at the current step.
   @returns The error in the second derivative of Kepler's equation.
 */
DACE::DA hypdMMdEE(const DACE::DA& F, const DACE::DA& ecc)
{
    return ecc * sinh(F);
}

/* @brief Wrapper function to convert from Mean anomaly to Eccentric Anomaly
   @param[in] M: Mean anomaly to solver for
   @param[in] ecc: Eccentricity to solve for
   @param[in] eps: Tolerance to be used in the solver.
   @returns The eccentric anomaly at that point.
*/
DACE::DA meanToEccentric(DACE::DA &M, DACE::DA &ecc, double eps)
{
    DACE::DA E = M;
    if ( DACE::cons(ecc) > 1 )
    {
        laguerreConway(M, E, ecc, M, eps, hypKeplersEquation, hypdMdE, hypdMMdEE);
    } else {
        laguerreConway(M, E, ecc, M, eps, keplersEquation, dMdE, dMMdEE);
    }
    return E;
}

/* @brief Wrapper function to convert from eccentric to true anomaly
   @param[in] E: Eccentric anomaly to convert
   @param[in] ecc: Eccentricity of orbit to convert
   @returns The true anomaly corresponding to the eccentric anomaly.
*/
template <typename Type>
Type eccentricToTrue(Type &E, Type &ecc)
{
    Type pi = 4.0 * atan(1.0);
    Type root = std::sqrt( (1-ecc) / (1+ecc) );
    Type inv = root * std::tan(E / 2.0);
    Type returnValue = std::fmod( (2.0 * std::atan(inv)), (2. * pi) );
    return returnValue;
}

/* Wrapper: Converts from Mean to True anomaly */
/* @brief Wrapper function to convert from mean anomaly to true anomaly.
   @param[in] M Mean anomaly to convert
   @param[in] ecc Eccentricity of orbit for which to convert
   @returns The true anomaly corresponding to the mean anomaly and eccentricity.
*/
template <typename Type>
void meanToTrue(Type &M, Type &ecc, Type &f)
{
    Type eps = 1.e-012;
    Type E =  meanToEccentric(M, ecc, eps);
    f = eccentricToTrue(E, ecc);
}

#endif
