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

/* @brief Compute the cross product of two vectors of arbitrary type and dimension.

   @param[in] vec1 The first component vector.
   @param[in] vec2 The second component vector.
   @param[out] result The cross product.
*/
template <typename T>
void cross3( const std::vector<T>& vec1, const std::vector<T>& vec2, std::vector<T>& result )
{
    /* Make sure vec1 and vec2 are of the same (correct) size.*/
#ifdef DEBUG
    if ( vec1.size() != 3 ) throw std::runtime_error("One of vec1 or vec2 are not of the same size.")
    if ( vec1.size() != vec2.size() ) throw std::runtime_error("vec1 and vec2 are not the same size in cross.\n");
#endif
    /* Resize result to ensure correct size & memory allocated. */
    result.resize( vec1.size() );
    /* Perform cross product */
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

/* @brief Compute the dot product of two arbitrarily-sized vectors.
   @param[in] vec1 First component vector.
   @param[in] vec2 Second component vector.
*/
template <typename T>
T dotProduct( const std::vector<T>& vec1, const std::vector<T>& vec2 )
{
#ifdef DEBUG
    if ( vec1.size() != vec2.size() ) throw std::runtime_error("One of vec1 or vec2 are not the same size in dotProduct.");
#endif
    T result = 0.0;
    for (size_t idx = 0; idx < vec1.size(); ++idx)
    {
        result += vec1[idx] * vec2[idx];
    }
    return result;
}

/* @brief Comput the Euclidean vector norm (length) for a given vector.
   @param[in] Desired vector to compute the norm for.
   @return Euclidean norm.
*/
template <typename T>
T vectorNorm( const std::vector<T>& vector )
{
    T result = 0.0;
    for (size_t idx = 0; idx < vector.size(); ++idx)
    {
        result += vector[idx] * vector[idx];
    }
    return sqrt( DACE::cons(result) );
}

/* @brief Compute the normalized version of an input vector.
   @param[in] Vector to compute the norm in.
   @param[out] Normalised vector
*/
template <typename T>
void normalizeVector( const std::vector<T>& vec1, std::vector<T>& unitVector )
{
    /* Ensure unit vector is the correct size. */
    unitVector.resize( vec1.size() );
    /* Get vector norm */
    T norm = vectorNorm( vec1 );
    for (size_t idx = 0; idx < vec1.size(); ++idx)
    {
        unitVector[idx] = vec1[idx] / norm;
    }
}

/* @brief Find the separation angle in radians between two double 
    precision, 3-dimensional vectors.  This angle is defined as zero
    if either vector is zero.
   @param[in] Vector to compute the angle of.
   @param[in] Vector to compute the angle of.
   @returns Angle between the two vectors.
*/
template <typename T>
DACE::DA separationAngle( const std::vector<T>& vec1, const std::vector<T>& vec2 )
{
    using DACE::DA;
    const DACE::DA pi = 4.0 * atan(1.0);
    DA v1Mag = vectorNorm( vec1 );
    DA v2Mag = vectorNorm( vec2 );
    DA retVal;
    if ( abs(v1Mag) == 0.)
    {
        retVal = 0.;
        return retVal;
    } else if (  abs(v2Mag) == 0.)
    {
        retVal = 0.;
        return retVal;
    }
    std::vector<T> vTemp(3);
    std::vector<T> v1Normal(3), v2Normal(3);
    normalizeVector(vec1, v1Normal);
    normalizeVector(vec2, v2Normal);
    DA v1DotV2 = dotProduct(v1Normal, v2Normal);
    if ( DACE::cons(v1DotV2) > 0. )
    {
        std::transform( v1Normal.begin(), v1Normal.end(), v2Normal.begin(), vTemp.begin(), std::minus<T>() );
        retVal = asin( vectorNorm(vTemp) * 0.5 ) * 2.;
    } else if ( DACE::cons(v1DotV2) < 0. )
    {
        std::transform( v1Normal.begin(), v1Normal.end(), v2Normal.begin(), vTemp.begin(), std::plus<T>() );
        retVal = pi - asin( vectorNorm(vTemp) * 0.5 ) * 2.;
    } else
    {
        retVal = pi / 2.0;
    }
    return retVal;
}

template <typename T, typename V>
void par2ic( const std::vector<T>& OEs, const V& mu, std::vector<T>& r0, std::vector<T>& v0)
{
    T r_p = OEs[0];
    T e = OEs[1];
    T i = OEs[2];
    T omg = OEs[3];
    T omp = OEs[4];
    T M = OEs[5];
    T b, n, xper, yper, xdotper, ydotper;
    T R[3][3];
    T cosomg, cosomp, sinomg, sinomp, cosi, sini;
    T dNdZeta;
    T a = r_p / (1 - e);
    T EA = meanToEccentric( M, e, 1e-012 ); // 1e-012 => tolerance

    /* SMA asumed positive */
    if ( DACE::cons(e) > 1)
    {
        a = -a;
    }
    if (DACE::cons(e) < 1.0)
    {
        b = a * sqrt( 1- e*e );
        n = mu / (a);
        n /= (a*a);
        n = sqrt(n);
        xper = a * (cos(EA) - e);
        yper = b * sin(EA);
        xdotper = -( a * n * sin(EA) ) / ( 1 - e * cos(EA) );
        ydotper = ( b * n * cos(EA) ) / ( 1 - e * cos(EA) );
    } else 
    {
        b = -a * sqrt( e * e - 1 );
        n = sqrt(-mu / ( a * a * a ));
        dNdZeta = e * ( 1 + tan(EA) * tan(EA) ) - (0.5 + 0.5 * pow(tan(0.5 * EA + M_PI_4), 2)) / tan(0.5 * EA + M_PI_4);
        xper = a / cos(EA) - a * e;
        yper = b * tan(EA);
        xdotper = a * tan(EA) / cos(EA) * n / dNdZeta;
        ydotper = b / pow(cos(EA), 2) * n / dNdZeta;
    }
    /* Construct rotation matrix */
    cosomg = cos(omg);
    cosomp = cos(omp);
    sinomg = sin(omg);
    sinomp = sin(omp);
    cosi = cos(i);
    sini = sin(i);
    R[0][0] = cosomg * cosomp - sinomg * sinomp * cosi;
    R[0][1] = -cosomg * sinomp - sinomg * cosomp * cosi;
    R[0][2] = sinomg * sini;
    R[1][0] = sinomg * cosomp + cosomg * sinomp * cosi;
    R[1][1] = -sinomg * sinomp + cosomg * cosomp * cosi;
    R[1][2] = -cosomg * sini;
    R[2][0] = sinomp * sini;
    R[2][1] = cosomp * sini;
    R[2][2] = cosi;
    /* Transform using naive matmul */
    r0.resize( 3 );
    v0.resize( 3 );
    std::fill( r0.begin(), r0.end(), 0 );
    std::fill( v0.begin(), v0.end(), 0 );
    T temp[3] = {xper, yper, 0.0};
    T temp2[3] = {xdotper, ydotper, 0.0};
    for (size_t j = 0; j < 3; ++j)
    {
        for (size_t k = 0; k < 3; ++k)
        {
            r0[j] += R[j][k] * temp[k];
            v0[j] += R[j][k] * temp2[k];
        }
    }
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
void ic2par(const std::vector<DACE::DA>& state, const double& mu, std::vector<DACE::DA>&OEs)
{
    using DACE::DA;
    /* Check for bad mu */
    if (mu <= 0.)
    {
        throw std::runtime_error("Negative value of mu passed to state_to_oes.");
    }
    const double pi = 4.0 * atan(1.0);
    std::vector<DA> r = {state[0], state[1], state[2]};
    std::vector<DA> v = {state[3], state[4], state[5]};
    DA rMag = vectorNorm( r );
    DA vMag = vectorNorm( v );
    std::vector< DA > h(3); // Angular momentum vector
    cross3( r, v, h );
    if ( abs( h[0] ) < 1e-14 && abs( h[1] ) < 1e-14 && abs( h[2] ) < 1e-14 )
    {
        throw std::runtime_error("r, v too close to parallel in state_to_oes.");
    }
    DA d_1 = -h[1];
    std::vector<DA> n = {d_1, h[0], 0.0};
    DA d_2 = vMag;
    d_1 = d_2 * d_2 - mu / rMag;
    DA d_3 = -dotProduct(r,  v);
    std::vector<DA> eccVector(3);
    for (unsigned i = 0; i < 3; ++i) eccVector[i] = d_1 * r[i] + d_3 * v[i];
    d_1 = 1./mu;
    for (unsigned i = 0; i < 3; ++i) eccVector[i] *= d_1;
    /* Determine size and shape of the orbit.
       The eccentricity of the orbit is the magnitude of the
       eccentricity vector. If the eccentricity is "close" to one,
       go ahead and make this a parabola.

       The perifocal distance depends on the eccentricity and the
       semi-latus rectum, which in turn orbit depends only on the
       specific angular momentum of the orbiting object. */
    d_1 = vectorNorm(eccVector);
    DA ecc;
    if ( abs( cons(d_1)  - 1 ) < 1e-10 ) ecc = 1.0;
    else ecc = d_1;
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
    if (abs (DACE::cons(inc)) < 1e-10) inc = 0.;
    else if ( abs(DACE::cons(inc) - pi < 1e-10 ) ) inc = pi;
    n = {1., 0, 0};

/*                                                              ^ */
/*     The longitude of the ascending node is the angle between i */
/*     (the x-axis) and the node vector, n. */
/*                                       - */
    DA lNode = atan2( n[1], n[0] );
    if ( DACE::cons(lNode) < 0.0 ) lNode += 2.0 * pi; // + 2pi
    DA argp;
    if ( DACE::cons(ecc) == 0.0 ) argp = 0.0;
    else
    {
        /* Magnitude first, sign next */
        argp = separationAngle(n, eccVector);
        if ( DACE::cons(argp) != 0)
        {
            if ( DACE::cons(inc) == 0.0 || DACE::cons(inc) == pi )
            {
                /* Quadrant of argp determined by component of e in direction of  h x n */
                std::vector<DA> hCrossN(3), hCrossNNormalized(3);
                cross3(h, n, hCrossN);
                normalizeVector(hCrossN, hCrossNNormalized);
                if ( DACE::cons( dotProduct(eccVector, hCrossNNormalized) ) < 0.)
                {
                    argp *= -1;
                    argp += 2. * pi;
                }
            } else if ( DACE::cons(eccVector[2]) < 0.){
                /* Periapsis below reference plane */
                argp *= -1;
                argp += 2.0 * pi;
            }
        }
    }
    std::vector<DA> perix(3);
    /* Get the true, then mean anomaly. */
    if ( DACE::cons(ecc) == 0.0 )
    {
        normalizeVector(n, perix);
    } else
    {
        normalizeVector(eccVector, perix);
    }
    std::vector<DA> periy(3);
    std::vector<DA> periyNormalized(3);
    cross3(h, perix, periy);
    normalizeVector(periy, periyNormalized);
    DA nu = atan2( dotProduct(r, periy), dotProduct(r, perix) );
    DA m0;
    /* Convert true anomaly to mean anomaly */
    if ( DACE::cons(ecc) < 1.0 )
    {
        DA cosea = (ecc + cos(nu)) / (ecc * cos(nu) + 1);
        DA sinea = rMag / r_p * sqrt( (1-ecc) / (ecc+1) ) * sin(nu);
        DA ea = atan2(sinea, cosea);
        d_1 = ea - ecc * sin(ea);
        if ( DACE::cons(d_1 * nu) < 0.)
        {
            m0 = -d_1;
        } else
        {
            m0 = d_1;
        }
        if ( DACE::cons(m0) < 0.) m0 += 2. * pi;
    } else if ( DACE::cons (ecc) > 1.0)
    {
        DA coshf = (ecc + cos(nu)) / (ecc * cos(nu) + 1);
        d_1 = std::max( 1.0, DACE::cons(coshf) );
        DA ea = acosh( d_1 );
        d_1 = ecc * sinh(ea) - ea;
        if ( DACE::cons(d_1 * nu) < 0.)
        {
            m0 = -d_1;
        } else 
        {
            m0 = d_1;
        }
    } else
    {
        DA ea = tan( nu / 2.0 );
        DA d_2 = ea;
        d_1 = ea + d_2 * (d_2 * d_2) / 3.;
        if ( DACE::cons(d_1 * nu) < 0 ) m0 = -d_1;
        else m0 = d_1;
    }
    OEs.resize(6);
    OEs[0] = r_p;
    OEs[1] = ecc;
    OEs[2] = inc;
    OEs[3] = lNode;
    OEs[4] = argp;
    OEs[5] = m0;
}
#endif
