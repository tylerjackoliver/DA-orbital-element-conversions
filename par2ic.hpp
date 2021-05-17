/****************************************************************************
 *   Copyright (C) 2020 J. Tyler                                             *
 *   University of Southampton, United Kingdom                               *
 *                                                                           *
 *   jack.tyler@soton.ac.uk                                                  *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 **************************************************************************/
#ifndef __ORBITAL_ELEMENT_CONVERSIONS_H__
#define __ORBITAL_ELEMENT_CONVERSIONS_H__

#include <cmath>
#include <vector>
#include <stdexcept>
#include "laguerreConway.hpp"
#include <dace/dace_s.h>
#include "linearAlgebra.h"
#include <algorithm>

extern "C"
{
#include <SpiceUsr.h>
}

DACE::DA max(const DACE::DA &in1, const DACE::DA &in2)
{
    return DACE::cons(in1) >= DACE::cons(in2) ? in1 : in2;
}

DACE::DA min(const DACE::DA &in1, const DACE::DA &in2)
{
    return DACE::cons(in1) <= DACE::cons(in2) ? in1 : in2;
}

/* @brief Find the separation angle in radians between two double 
    precision, 3-dimensional vectors.  This angle is defined as zero
    if either vector is zero.
   @param[in] Vector to compute the angle of.
   @param[in] Vector to compute the angle of.
   @returns Angle between the two vectors.
*/
template <typename T>
DACE::DA separationAngle(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    using DACE::DA;
    const DACE::DA pi = 4.0 * atan(1.0);
    DA v1Mag = vectorNorm(vec1);
    DA v2Mag = vectorNorm(vec2);
    DA retVal;
    if ( fabs(DACE::cons(v1Mag)) == 0. || fabs(DACE::cons(v2Mag)) == 0)
    {
        retVal = 0.;
        return retVal;
    }
    std::vector<T> vTemp(3);
    std::vector<T> v1Normal(3), v2Normal(3);
    normalizeVector(vec1, v1Normal);
    normalizeVector(vec2, v2Normal);
    DA v1DotV2 = dotProduct(v1Normal, v2Normal);
    if (DACE::cons(v1DotV2) > 0.)
    {
        std::transform(v1Normal.begin(), v1Normal.end(), v2Normal.begin(), vTemp.begin(), std::minus<T>());
        retVal = asin(vectorNorm(vTemp) * 0.5) * 2.;
    }
    else if (DACE::cons(v1DotV2) < 0.)
    {
        std::transform(v1Normal.begin(), v1Normal.end(), v2Normal.begin(), vTemp.begin(), std::plus<T>());
        retVal = pi - asin(vectorNorm(vTemp) * 0.5) * 2.;
    }
    else
    {
        retVal = pi / 2.0;
    }
    return retVal;
}

/*      From SPICELIB: The new method uses the criterion: for inclination zero or pi
        the argument of periapse is between zero and pi radians when

            e  *  ( h x n )  >  0
            -       -   -    -

        where

            e  is the eccentricity vector,
            -

            h  is the specific angular momentum vector,
            -

            n  is the node vector.
            -

        The computation of M0 was re-coded for improved accuracy.
        The new computation uses ATAN2 rather than ACOS to find
        the eccentric anomaly for the ellipse case.  The quadrant
        of M0 is now found by converting the position to the
        perifocal frame and finding the corresponding longitude.

        The old method, using the sign of <r,v>, did not work
        for circular orbits and was unreliable for near-circular
        orbits.

        Inclination is now computed using VSEP.

        Also, added error checks for non-positive MU, zero
        position, velocity, and specific angular momentum vectors.
*/
void state_to_oes(const std::vector<DACE::DA> &state, const double &mu, std::vector<DACE::DA> &OEs)
{
    using DACE::DA;
    const double pi = 4.0 * atan(1.0);
    /* Check for bad mu */
    if (mu <= 0.)
    {
        throw std::runtime_error("Negative value of mu passed to state_to_oes.");
    }
    std::vector<DA> r = {state[0], state[1], state[2]};
    std::vector<DA> v = {state[3], state[4], state[5]};
    DA rMag = vectorNorm(r);
    DA vMag = vectorNorm(v);
    std::vector<DA> h(3); // Angular momentum vector
    cross3(r, v, h);
    if ( fabs( DACE::cons(h[0]) ) < 1e-14 && fabs( DACE::cons(h[1]) ) < 1e-14 && fabs( DACE::cons(h[2]) ) < 1e-14)
    {
        throw std::runtime_error("r, v too close to parallel in state_to_oes.");
    }
    DA d_1 = -h[1];
    std::vector<DA> n = {d_1, h[0], 0.0};
    DA d_2 = vMag;
    d_1 = d_2 * d_2 - mu / rMag;
    DA d_3 = -dotProduct(r, v);
    std::vector<DA> eccVector(3);
    for (unsigned i = 0; i < 3; ++i)
        eccVector[i] = d_1 * r[i] + d_3 * v[i];
    d_1 = 1. / mu;
    for (unsigned i = 0; i < 3; ++i)
        eccVector[i] *= d_1;
    /* Determine size and shape of the orbit.
       The eccentricity of the orbit is the magnitude of the
       eccentricity vector. If the eccentricity is "close" to one,
       go ahead and make this a parabola.

       The perifocal distance depends on the eccentricity and the
       semi-latus rectum, which in turn orbit depends only on the
       specific angular momentum of the orbiting object. */
    d_1 = vectorNorm(eccVector);
    DA ecc;
    /* Modified 17/05/2021: want eccentricity to not be rounded - we're losing derivatives here! */
//    if ( fabs(cons(d_1) - 1) < 1e-10)
//        ecc = 1.0;
//    else
//        ecc = d_1;
    ecc = d_1;
    /* </end modify */
    DA p = dotProduct(h, h) / mu;
    DA r_p = p / (ecc + 1);
    /*     Next, the orientation of the orbit. 
                                                           ^ 
         The inclination of the orbit is the angle between k (which is 
         perpendicular to the equator) and h (which is perpendicular to 
         the orbit.                        -

         If close to zero or pi, make it exact. In either case, the node 
         vector becomes undefined. 
    */
    std::vector<DA> zVector = {0., 0., 1.};
    DA inc = separationAngle(h, zVector);
    /* Modification 17/05/2021: Don't change the constant part as this removes derivatives.*/
//    if ( fabs(DACE::cons(inc)) < 1e-10) {
//        inc = 0.;
//        n = {1., 0, 0};
//    }
//    else if ( fabs( fabs(DACE::cons(inc)) - pi ) < 1e-10) {
//        inc = pi;
//        n = {1., 0, 0};
//    }
    if ( fabs(DACE::cons(inc)) < 1e-010 || fabs( fabs(DACE::cons(inc)) - pi ) < 1e-010)
    {
        n = {1., 0., 0.};
    }

    /*                                                              ^ */
    /*     The longitude of the ascending node is the angle between i */
    /*     (the x-axis) and the node vector, n. */
    /*                                       - */
    DA lNode = atan2(n[1], n[0]);
    if (DACE::cons(lNode) < 0.0)
        lNode += 2.0 * pi; // + 2pi
    DA argp;
    if (DACE::cons(ecc) == 0.0) {
        argp = 0.0;
    }
    else
    {
        /* Magnitude first, sign next */
        argp = separationAngle(n, eccVector);
        if (DACE::cons(argp) != 0)
        {
            if (DACE::cons(inc) == 0.0 || DACE::cons(inc) == pi)
            {
                /* Quadrant of argp determined by component of e in direction of  h x n */
                std::vector<DA> hCrossN(3), hCrossNNormalized(3);
                cross3(h, n, hCrossN);
                normalizeVector(hCrossN, hCrossNNormalized);
                if (DACE::cons(dotProduct(eccVector, hCrossNNormalized)) < 0.)
                {
                    argp *= -1;
                    argp += 2. * pi;
                }
            }
            else if (DACE::cons(eccVector[2]) < 0.)
            {
                /* Periapsis below reference plane */
                argp *= -1;
                argp += 2.0 * pi;
            }
        }
    }
    std::vector<DA> perix(3);
    /* Get the true, then mean anomaly. */
    if (DACE::cons(ecc) == 0.0)
    {
        normalizeVector(n, perix);
    }
    else
    {
        normalizeVector(eccVector, perix);
    }
    std::vector<DA> periy(3);
    std::vector<DA> periyNormalized(3);
    cross3(h, perix, periy);
    normalizeVector(periy, periyNormalized);
    DA nu = atan2(dotProduct(r, periy), dotProduct(r, perix));
    DA m0;
    /* Convert true anomaly to mean anomaly */
    if (DACE::cons(ecc) < 1.0)
    {
        DA cosea = (ecc + cos(nu)) / (ecc * cos(nu) + 1);
        DA sinea = rMag / r_p * sqrt((1 - ecc) / (ecc + 1)) * sin(nu);
        DA ea = atan2(sinea, cosea);
        d_1 = ea - ecc * sin(ea);
        if (DACE::cons(d_1 * nu) < 0.)
        {
            m0 = -d_1;
        }
        else
        {
            m0 = d_1;
        }
        if (DACE::cons(m0) < 0.)
            m0 += 2. * pi;
    }
    else if (DACE::cons(ecc) > 1.0)
    {
        DA coshf = (ecc + cos(nu)) / (ecc * cos(nu) + 1);
        d_1 = max(1.0, coshf);
        DA ea = acosh(d_1);
        d_1 = ecc * sinh(ea) - ea;
        if (DACE::cons(d_1 * nu) < 0.)
        {
            m0 = -d_1;
        }
        else
        {
            m0 = d_1;
        }
    }
    else
    {
        DA ea = tan(nu / 2.0);
        d_2 = ea;
        d_1 = ea + d_2 * (d_2 * d_2) / 3.;
        if (DACE::cons(d_1 * nu) < 0)
            m0 = -d_1;
        else
            m0 = d_1;
    }
    OEs.resize(6);
    OEs[0] = r_p;
    OEs[1] = ecc;
    OEs[2] = inc;
    OEs[3] = lNode;
    OEs[4] = argp;
    OEs[5] = m0;
}

void stmp03_(const DACE::DA& x, DACE::DA& c0, DACE::DA& c1, DACE::DA &c2, DACE::DA &c3)
{
    using DACE::DA;
    /* If |magnitude| > 1, compute directly from cosh etc */
    if ( fabs( DACE::cons(x) ) > 1 )
    {
        DA z = sqrt( -x );
        c0 = cosh( z );
        c1 = sinh(z) / z;
        c2 = (1-c0) / x;
        c3 = (1-c1) / x;
    }
    else
    {
        c3 = 1.;
        for (int i = 19; i > 3; i -= 2)
        {
            c3 = 1 - x/( i * (i-1) ) * c3;
        }
        c3 *= (1./6);
        c2 = 1.;
        for (int i = 18; i > 2; i -= 2)
        {
            c2 = 1. - c2 * x / ( i * (i-1) );
        }
        c2 /= 2;
        c1 = 1 - x * c3;
        c0 = 1 - x * c2;
    }
}

void stmp03(const DACE::DA& x, DACE::DA &c0, DACE::DA &c1, DACE::DA &c2, DACE::DA&c3)
{
    stmp03_(x, c0, c1, c2, c3);
}

DACE::DA brcktd(DACE::DA &number, DACE::DA &end1, DACE::DA &end2)
{
    /*      Bracket a number. That is, given a number and an acceptable */
    /*      interval, make sure that the number is contained in the */
    /*      interval. (If the number is already in the interval, leave it */
    /*      alone. If not, set it to the nearest endpoint of the interval.) */
    if (DACE::cons(end1) < DACE::cons(end2))
    {
        DACE::DA b = DACE::cons(end2) <= DACE::cons(number) ? end2 : number;
        return max(end1, b);
    }
    else
    {
        DACE::DA b = min(end1, number);
        return max(end2, b);
    }
}

void prop2b(const double &mu, const std::vector<DACE::DA> &pState, const double dt, std::vector<DACE::DA> &state)
{
    using DACE::DA;
    /* First, perform error checking */
    std::vector<DACE::DA> pos(3), vel(3);
    for (int i = 0; i < 3; ++i)
    {
        pos[i] = pState[i];
        vel[i] = pState[i + 3];
    }
    if (mu <= 0)
        throw std::runtime_error("A non-positive value of mu was passed into prop2b.");
    if (DACE::cons(vectorNorm(pos)) == 0)
        throw std::runtime_error("A zero-magnitude position vector was passed into prop2b.");
    if (DACE::cons(vectorNorm(vel)) == 0)
        throw std::runtime_error("A zero-magnitude velocity vector was passed into prop2b.");
    /*        Obvious problems have been checked. Here are the relevant */
    /*        equations. Let ... */

    /*           GM        be the gravitational attraction of the central */
    /*                     mass. */

    /*           POS and   be the initial position and velocity respectively */
    /*           VEL       of the orbiting object. */

    /*           R0       be the magnitude of the position vector POS */

    /*           RV       be the value of the dot product  POS * VEL */
    DA r0 = vectorNorm(pos);
    DA rv = dotProduct(pos, vel);
    /*        Let HVEC be the specific angular momentum vector and let Q be */
    /*        the distance at periapse. */

    /*                   1)    HVEC  =   POS  x  VEL */

    /*                                       2 */
    /*                   2)    H2    = |HVEC|  =  GM*(1+E)*Q */
    std::vector<DA> hVec(3);
    cross3(pos, vel, hVec);
    DA h2 = dotProduct(hVec, hVec);
    if (DACE::cons(h2) == 0.)
    {
        throw std::runtime_error("Magnitude of the angular momentum vector is zero in prop2b.");
    }
    std::vector<DA> tmpVector(3);
    cross3(vel, hVec, tmpVector);
    DA d_1 = 1. / mu;
    DA d_2 = -1. / r0;
    std::vector<DA> eqVec(3);
    for (int i = 0; i < 3; ++i)
    {
        eqVec[i] = d_1 * tmpVector[i] + d_2 * pos[i];
    }
    DA e = vectorNorm(eqVec);
    DA q = h2 / (mu * (e + 1));
    /*        From the discussion of the universal variables formulation in */
    /*        Danby's book on pages 174 and 175 (see the reference listed */
    /*        above) you can show that by making the substitutions */

    /*              F  =  1 - E */

    /*        and */

    /*                       _____ */
    /*                      /  Q */
    /*              S  =   / -----    X   = B * X */
    /*                   \/   GM */

    /*        that DT satisfies the universal variables Kepler's equation: */

    /*                                   2     2     2        2 */
    /*              DT =  B*R0*X*C_1( F*X ) + B *RV*X C_2( F*X ) */

    /*                                               3        2 */
    /*                                      +   B*Q*X C_3( F*X ) */

    /*                 =  KFUN( X ) */

    /*        (where C_k is used to denote the Stumpff functions. This is */
    /*        the universal variables formulation of Kepler's equation. */
    /*        KFUN is our abbreviation for "Kepler function.") */

    /*        (One might wonder, "Why make such a change of variables?" */
    /*        By making this substitution early in the derivation supplied */
    /*        in Danby's book, you can always deal with physically */
    /*        meaningful quantities --- the pure numeric value of F and the */
    /*        distance of periapse.  Thus one does not need to be concerned */
    /*        about infinite or negative semi-major axes or with discussing */
    /*        how to interpret these somewhat artificial artifacts of the */
    /*        classical derivations for two body motion.) */

    /*        Given the unique X for which this Kepler's equation is */
    /*        satisfied, we can compute the state of the orbiting object */
    /*        at a time DT past the epoch of the state POS and VEL. */
    /*        Evidently we will need the constants: */
    DA f = 1. - e;
    DA b = sqrt(q / mu);
    DA br0 = b * r0;
    DA b2rv = b * b * rv;
    DA bq = b * q;
    /*        The state corresponding to the value of X that solves this */
    /*        equation is given by */

    /*              PC * POS + VC * VEL              ( position ) */

    /*        and */

    /*              PCDOT * POS + VCDOT * VEL        ( velocity ) */

    /*        where */
    /*                                            2        2 */
    /*           ( 1 )    PC    =  1  -  ( Q/R0 )X C_2( F*X ) */

    /*                                            3        2 */
    /*           ( 2 )    VC    =  DT -  ( B*Q  )X C_3( F*X ) */

    /*                                       Q               2 */
    /*           ( 3 )    PCDOT =     -  ( ------ ) X C_1( F*X ) */
    /*                                     B*R*R0 */

    /*                                      B*Q     2        2 */
    /*           ( 4 )    VCDOT =  1  -  (  ---  ) X C_2( F*X ) */
    /*                                      B*R */

    /*        Here R denotes the distance from the center of CP*POS + CV*VEL */
    /*        It turns out that R can be computed as: */

    /*                                        2     2             2 */
    /*           ( 5 )   B*R    = B*R0 C_0(F*X ) + B *RV X C_1(F*X ) */

    /*                                                 2       2 */
    /*                                        +   B*Q X C_2(F*X ) */

    /*        Therefore we will also need the constant */
    DA qovr0 = q / r0;
    /*        We will have to find the unique value of X such that */

    /*             DT = KFUN ( X ) */

    /*        where KFUN stands for the "Kepler function" defined by the */
    /*        equation below: */

    /*                                   2 */
    /*        KFUN(X) =   B*R0*X * C_1(FX ) */

    /*                   2     2        2 */
    /*                + B *RV*X * C_2(FX ) */

    /*                         3        2 */
    /*                +   B*Q*X * C_3(FX ) */

    /*        (There is a unique solution to this equation. KFUN(X) is */
    /*        unbounded above and below and is an increasing function */
    /*        over all real X for all non-rectilinear orbits. To see this */
    /*        we note that the variable X is a function of DT and is given */
    /*        by the integral from 0 to DT of the differential: */

    /*                   dt */
    /*                 ------ */
    /*                 B*R(t) */

    /*        where R(t) is the range of the body as a function of time. */
    /*        Therefore X is an increasing function of DT, and DT must */
    /*        also be an increasing function of X. */

    /*        Thus, there is a unique value of X  that solves this */
    /*        equation). */

    /*        If F is less than zero, we can have the computation of C0,... */
    /*        overflow.  This is because for X < 0 */

    /*               C_0(X) = COSH( DSQRT(-X) ) */

    /*               C_1(X) = SINH( DSQRT(-X) ) */
    /*                        ----------------- */
    /*                              DSQRT(-X) */

    /*        and from the recursion relationship we know that */

    /*               C_2(X) =  ( 1/0! - C_0(X) ) / X */

    /*               C_3(X) =  ( 1/1! - C_1(X) ) / X */

    /*                         1 - COSH( DSQRT(-X) ) */
    /*               C_2(X) = ------------------------ */
    /*                                  X */

    /*                         1  - SINH( DSQRT(-X) ) / DSQRT(-X) */
    /*               C_3(X) = ----------------------------------- */
    /*                                    X */

    /*        Clearly for negative values of F*X*X having large magnitude, */
    /*        it is easy to get an overflow. */

    /*        In the case when F is less than 0 we choose X so that we can */
    /*        compute all of the following: */

    /*               | COEF_0 * X**0 * C_0(FX**2) | */

    /*               | COEF_1 * X**1 * C_1(FX**2) | */

    /*               | COEF_2 * X**2 * C_2(FX**2) | */

    /*               | COEF_3 * X**3 * C_3(FX**2) | */

    /*         where COEF_n are coefficients that will be used in forming */
    /*         linear combinations of X**n C_n(FX**2) terms. */

    /*         The variable portion of the last 3 terms above can be */
    /*         rewritten as: */

    /*                                   SINH ( DSQRT(-F)*|X| ) */
    /*        | X**1 * C_1(FX**2) |  =   ---------------------- */
    /*                                          DSQRT(-F) */

    /*                                   1 - COSH( DSQRT(-F)*|X| ) */
    /*        | X**2 * C_2(FX**2) |  =  ---------------------------- */
    /*                                             -F */

    /*                                  DSQRT(-F)*|X|   - SINH(DSQRT(-F)*|X|) */
    /*        | X**3 * C_3(FX**2) |  =  ------------------------------------- */
    /*                                              F*DSQRT(-F) */

    /*        For large |X| the absolute values of these expressions are well */
    /*        approximated by */

    /*                                         0.0 */
    /*               COSH( DSQRT(-F)|X| ) * |F| */

    /*                                         -0.5 */
    /*               SINH( DSQRT(-F)|X| ) * |F| */

    /*                                         -1.0 */
    /*               COSH( DSQRT(-F)|X| ) * |F| */

    /*                                         -1.5 */
    /*               SINH( DSQRT(-F)|X| ) * |F| */

    /*        For large |X| the logarithms of these expressions are well */
    /*        approximated by: */

    /*               DSQRT(-F)|X| - LOG(2) - 0.0*LOG(-F) */

    /*               DSQRT(-F)|X| - LOG(2) - 0.5*LOG(-F) */

    /*               DSQRT(-F)|X| - LOG(2) - 1.0*LOG(-F) */

    /*               DSQRT(-F)|X| - LOG(2) - 1.5*LOG(-F) */

    /*        respectively. */

    /*        To ensure that we can form a linear combination of these terms */
    /*        we will require that: */

    /*           |COEF_N*X**N * C_N(FX**2)| < DPMAX / 4 */

    /*        for N=0,1,2,3.  This is equivalent to */

    /*              LOG ( X**N * C_N(FX**2) )   <      LOG ( DPMAX ) */
    /*            + LOG (|COEF_N|)                   - 2 LOG ( 2     ) */

    /*        or */

    /*              LOG ( X**N * C_N(FX**2) )   <      LOG ( DPMAX    ) */
    /*                                             -   LOG ( |COEF_N| ) */
    /*                                             - 2*LOG ( 2        ). */

    /*        Replacing the left hand side with the magnitude expressions */
    /*        computed above we have: */

    /*            DSQRT(-F)|X| - LOG(2) - N*0.5*LOG( -F )  <   LOG ( DPMAX  ) */
    /*                                                      -  LOG (|COEF_N|) */
    /*                                                      -2*LOG ( 2      ) */

    /*         So that: */

    /*            |X|  <    {   LOG ( DPMAX  ) */
    /*                        - LOG (|COEF_N|) */
    /*                        - LOG (  2     ) */
    /*                        + LOG ( -F     )*N*0.5 } / DSQRT(-F) */

    /*         Let MAXC be the maximum of 1.0D0 and the various coefficients */
    /*         of the Stumpff functions.  We can then set our absolute value */
    /*         bound on X to be: */

    /*             MIN        LOG(DPMAX/2) - LOG(MAXC) + (n/2)LOG(-F) */
    /*            n = 0,3  {  -----------------------------------------  } */
    /*                               DSQRT(-F) */

    /*        (Actually we know that the minimum must occur for n = 0 or */
    /*        for n = 3). */
    DA d_3;
    d_2 = 1., d_3 = abs(br0), d_2 = max(d_2, d_3), d_3 = abs(b2rv),
    d_2 = max(d_2, d_3), d_3 = abs(bq), d_2 = max(d_2, d_3),
    d_3 = (d_1 = qovr0 / bq, abs(d_1));
    DA maxc = max(d_2, d_3);
    DA logmxc, logdpm, fixed, rootf, logf, bound, logbnd;
    if (DACE::cons(f) < 0.)
    {
        logmxc = log(maxc);
        logdpm = log(dpmax_() / 2.);
        fixed = logdpm - logmxc;
        rootf = sqrt(-f);
        logf = log(-f);
        /* Computing MIN */
        d_1 = fixed / rootf, d_2 = (fixed + logf * 1.5) / rootf;
        bound = min(d_1, d_2);

        /*           Note that in the above, we can always perform the division */
        /*           by ROOTF.  To see this we note that -F is at least the */
        /*           machine precision (we got it by subtracting E from 1.) */
        /*           Thus its square root is a reasonably large number (if F is */
        /*           10**-N then ROOTF is 10**(-N/2) )  The value of FIXED is */
        /*           about 3*M where M is the largest exponent such that 2**M */
        /*           is representable on the host machine.  Thus BOUND is at */
        /*           worst M*10**(N/2)  This will always be computable. */
    }
    else
    {

        /*           In the case when F is non-negative we must be sure we */
        /*           can compute all of the following. */

        /*               | COEF_0 * X**0 * C_0(FX**2) | < | COEF_0          | */

        /*               | COEF_1 * X**1 * C_1(FX**2) | < | COEF_1*|X|      | */

        /*               | COEF_2 * X**2 * C_2(FX**2) | < | COEF_2*X**2 / 2 | */

        /*               | COEF_3 * X**3 * C_3(FX**2) | < | COEF_3*X**3 / 6 | */

        /*           If we assume that COEF_0 is computable, all of these are */
        /*           bounded above by: */

        /*                       | MAX(COEF_1,...COEF_3) * X**3 / 6 | */

        /*           We want to make sure we can add these terms so we need to */
        /*           make sure that */

        /*              | MAX(COEF_1,...,COEF_3) * X**3 / 6 | < DPMAX() / 4. */

        /*           Thus we need: */

        /*              |X**3| <          1.5*DPMAX / MAX(COEF_1,...,COEF_3) */
        /*              |X|    <  DCBRT ( 1.5*DPMAX / MAX(COEF_1,...,COEF_3) ) */

        /*           (We'll use logarithms to compute the upper bound for |X|.) */

        logbnd = (log(1.5) + log(dpmax_()) - log(maxc)) / 3.;
        bound = exp(logbnd);
    }
    /*     We are now ready to find the unique value of X such that */

    /*             DT = KFUN ( X ) */

    /*     First we must bracket the root. The basic idea is this: */

    /*     1) KFUN(0) = 0 so we will let one endpoint of our initial */
    /*        guess of a bracketing interval be 0. */

    /*     2) We get our initial guess at the other endpoint of the */
    /*        bracketing interval by recalling that */

    /*                   dt */
    /*         dX  =   ------ */
    /*                 B*R(t) */

    /*        From this observation it follows that */

    /*                   DT */
    /*          X  <  ------- */
    /*                   B*Q */

    /*        Thus the solution to */

    /*             DT = KFUN ( X ) */

    /*        Satisifies */

    /*                     DT */
    /*         0 < X  <  ------- */
    /*                    B*Q */

    /*        We now have a guess at a bracketing interval. In the case */
    /*        DT is positive it looks like */

    /*                0        X */
    /*         -------[--------]----------------------------- */

    /*        This is ok mathematically, but due to rounding etc it is */
    /*        conceivable that we might not have bracketed the root. */
    /*        We check and if not we will double the */
    /*        endpoint farthest from zero and call this X, and make */
    /*        the other endpoint the old value of X. */

    /*                0 */
    /*         -------+--------[--------]-------------------- */

    /*        We continue this process ... */

    /*                0 */
    /*         -------+-----------------[-----------------]-- */

    /*        ...until the root is bracketed. (One shift is certain */
    /*        to do the job). */

    /*        If we perform this interval shift, we will have to take */
    /*        care that X does not run out of the domain for which */
    /*        we can safely compute KFUN.  Thus we will make sure that */
    /*        the endpoints of these shifted intervals always stay safely */
    /*        inside the domain for which KFUN can be computed. */
    DA x = dt / bq;
    d_1 = -bound;
    x = brcktd(x, d_1, bound);
    DA fx2 = f * x * x;
    DA c0, c1, c2, c3;
    stmp03( fx2, c0, c1, c2, c3 );
    DA kFun = x * (br0 * c1 + x * (b2rv * c2 + x * (bq * c3)));
    DA upper = 0.;
    DA lower = x;
    if (dt < 0.)
    {

        while (DACE::cons(kFun) > dt)
        {
            upper = lower;
            lower *= 2;
            DA oldx = x;
            d_1 = -bound;
            x = brcktd(lower, d_1, bound);
            /*           Make sure we are making progress. (In other words make sure */
            /*           we don't run into the boundary of values that X can assume. */
            /*           If we do run into the boundary, X will be unchanged and */
            /*           there's nothing further we can do.  We'll have to call it */
            /*           quits and tell the user what happened.) */
            if (DACE::cons(x) == DACE::cons(oldx))
            {
                throw std::runtime_error("Propagation not possible.");
            }
            fx2 = f * x * x;
            stmp03( fx2, c0, c1, c2, c3);
            kFun = x * (br0 * c1 + x * (b2rv * c2 + x * (bq * c3)));
        }
    }
    else if (dt > 0)
    {
        lower = 0.;
        upper = x;
        DA oldx = x;
        d_1 = -bound;
        x = brcktd(upper, d_1, bound);
        if (DACE::cons(x) == DACE::cons(oldx))
        {
            throw std::runtime_error("Propagation not possible.");
        }
        fx2 = f * x * x;
        stmp03( fx2, c0, c1, c2, c3);
    }
    else
    {
        state = pState;
        return;
    }
    /*     Ok. We've bracketed the root.  Now for lack of anything more */
    /*     clever, we just bisect to find the solution. */

    /*     We add a loop counter so that we can ensure termination of the */
    /*     loop below. */

    /*     On some systems the computed midpoint is stored in an extended */
    /*     precision register.  Thus the midpoint is always different from */
    /*     UPPER and LOWER.  Yet when the new value of LOWER and UPPER */
    /*     are assigned UPPER and LOWER do not change and hence the */
    /*     loop fails to terminate.  With the loop counter we force */
    /*     termination of the loop. */
    DA d_4;
    d_3 = lower, d_4 = (lower + upper) / 2.;
    d_1 = upper, d_2 = max(d_3, d_4);
    x = min(d_1, d_2);
    fx2 = f * x * x;
    stmp03( fx2, c0, c1, c2, c3);
    int lCount = 0;
    int mostC = 1000;
    while (DACE::cons(x) > DACE::cons(lower) && DACE::cons(x) < DACE::cons(upper) && lCount < mostC)
    {
        kFun = x * (br0 * c1 + x * (b2rv * c2 + x * bq * c3));
        if (DACE::cons(kFun) > dt)
        {
            upper = x;
        }
        else if (DACE::cons(kFun) < dt)
        {
            lower = x;
        }
        else
        {
            upper = x;
            lower = x;
        }
        /*        As soon as the bracketting values move away from */
        /*        zero we can modify the count limit. */
        if (mostC > 64)
        {
            if (DACE::cons(upper) != 0. && DACE::cons(lower) != 0.)
            {
                mostC = 64;
                lCount = 0;
            }
        }
        /* Computing MIN */
        /* Computing MAX */
        d_3 = lower, d_4 = (lower + upper) / 2.;
        d_1 = upper, d_2 = max(d_3, d_4);
        x = min(d_1, d_2);
        fx2 = f * x * x;
        stmp03(fx2, c0, c1, c2, c3);
        ++lCount;
    }
    DA x2 = x * x;
    DA x3 = x2 * x;
    DA br = br0 * c0 + x * (b2rv * c1 + x * (bq * c2));
    DA pc = 1. - qovr0 * x2 * c2;
    DA vc = dt - bq * x3 * c3;
    DA pcdot = -(qovr0 / br) * x * c1;
    DA vcdot = 1. - bq / br * x2 * c2;

    for (int i = 0; i < 3; ++i)
    {
        state[i] = pc * pos[i] + vc * vel[i];
        state[i + 3] = pcdot * pos[i] + vcdot * vel[i];
    }
}

void conics(std::vector<DACE::DA> &OEs, double &mu, std::vector<DACE::DA> &r, std::vector<DACE::DA> &v)
{
    using DACE::DA;
    const double pi = 4.0 * atan(1.0);
    DA r_p = OEs[0];
    DA ecc = OEs[1];
    DA inc = OEs[2];
    DA lNode = OEs[3];
    DA argP = OEs[4];
    DA M0 = OEs[5];

    /* Handle exceptions */
    if (DACE::cons(ecc) < 0)
    {
        throw std::runtime_error("The eccentricity supplied was negative.");
    }
    if (DACE::cons(r_p) <= 0)
    {
        throw std::runtime_error("The radius of periapsis was not positive.");
    }
    if (mu <= 0)
    {
        throw std::runtime_error("Value of mu provided was non-positive.");
    }
    /* Construct orthonormal basis vector */
    DA cosi = DACE::cos(inc);
    DA sini = DACE::sin(inc);
    DA cosn = DACE::cos(lNode);
    DA sinn = DACE::sin(lNode);
    DA cosw = DACE::cos(argP);
    DA sinw = DACE::sin(argP);
    DA snci = sinn * cosi;
    DA cnci = cosn * cosi;
    std::vector<DA> basisP(3);
    std::vector<DA> basisQ(3);
    basisP[0] = cosn * cosw - snci * sinw;
    basisP[1] = sinn * cosw + cnci * sinw;
    basisP[2] = sini * sinw;
    basisQ[0] = -cosn * sinw - snci * cosw;
    basisQ[1] = -sinn * sinw + cnci * cosw;
    basisQ[2] = sini * cosw;

    /*     Next construct the state at periapse. */

    /*     The position at periapse is just BASISP scaled by the distance */
    /*     at periapse. */

    /*     The velocity must be constructed so that we can get an orbit */
    /*     of this shape.  Recall that the magnitude of the specific angular */
    /*     momentum vector is given by DSQRT ( MU*RP*(1+ECC) ) */
    /*     The velocity will be given by V * BASISQ.  But we must have the */
    /*     magnitude of the cross product of position and velocity be */
    /*     equal to DSQRT ( MU*RP*(1+ECC) ). So we must have */

    /*        RP*V = DSQRT( MU*RP*(1+ECC) ) */

    /*     so that: */
    std::vector<DA> pState(6);
    DA velocity = sqrt(mu * (ecc + 1) / r_p);
    for (int i = 0; i < 3; ++i)
    {
        pState[i] = r_p * basisP[i];
        pState[i + 3] = velocity * basisQ[i];
    }

    /*     Finally compute DT the elapsed time since the epoch of periapse. */
    /*     Ellipses first, since they are the most common. */
    DA d_1; // Holding variable
    DA dt;
    DA ainvrs;
    DA n;
    DA period;
    if (DACE::cons(ecc) < 1.)
    {
        /*        Recall that: */

        /*        N ( mean motion ) is given by DSQRT( MU / A**3 ). */
        /*        But since, A = RP / ( 1 - ECC ) ... */
        ainvrs = (1. - ecc) / r_p;
        n = sqrt( mu * ainvrs) * ainvrs;
        period = 2. * pi / n;
        /*        In general the mean anomaly is given by */

        /*           M  = (T - TP) * N */

        /*        Where TP is the time of periapse passage.  M0 is the mean */
        /*        anomaly at time T0 so that */
        /*        Thus */

        /*           M0 = ( T0 - TP ) * N */

        /*        So TP = T0-M0/N hence the time since periapse at time ET */
        /*        is given by ET - T0 + M0/N.  Finally, since elliptic orbits are */
        /*        periodic, we can mod this value by the period of the orbit. */
        d_1 = M0 / n;
        dt = mod(d_1, DACE::cons(period));
    }
    else if (DACE::cons(ecc) > 1.)
    {
        ainvrs = (ecc - 1.) / r_p;
        n = sqrt(mu * ainvrs) * ainvrs;
        dt = M0 / n;
    }
    else
    {
        n = sqrt(mu / (r_p * 2.)) / r_p;
        dt = M0 / n;
    }
    std::vector<DA> state(6);
    prop2b(mu, pState, DACE::cons(dt), state);
    for (int i = 0; i < 3; ++i)
    {
        r[i] = state[i];
        v[i] = state[i + 3];
    }
}

void conics(std::vector<DACE::DA>& OEs, double& mu, std::vector<DACE::DA>& x)
{
    std::vector<DACE::DA> r(3), v(3);
    conics(OEs, mu, r, v);
    if ( x.size() != (r.size() + v.size())) x.resize( r.size() + v.size() );
    for (uint8_t i = 0; i < 3; ++i)
    {
        x[i] = r[i];
        x[i+3] = v[i];
    }
}

void conics( DACE::AlgebraicVector<DACE::DA>& OEs, double& mu, DACE::AlgebraicVector<DACE::DA>& x)
{
    std::vector<DACE::DA> OEsVector( OEs.size() );
    std::vector<DACE::DA> xVector( x.size() );
    for (unsigned long i = 0; i < OEs.size(); ++i) OEsVector[i] = OEs[i];
    conics( OEsVector, mu, xVector );
    if ( x.size() != xVector.size() ) x.resize( xVector.size() );
    for (unsigned long i = 0; i < xVector.size(); ++i) x[i] = xVector[i];
}

#endif
