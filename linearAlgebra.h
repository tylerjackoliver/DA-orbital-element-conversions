/*
 * Created by Jack Tyler on 14/05/2021.
 * Copyright (c) 2021 University of Southampton. All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef ACROBAT_LINEARALGEBRA_H
#define ACROBAT_LINEARALGEBRA_H


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
    return sqrt( result );
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

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T>& vec)
{
    for (size_t i = 0; i < vec.size()-1; ++i)
    {
        out << vec[i] << ", ";
    }
    return out << vec[ vec.size() - 1 ];
}

/* @brief Sometimes, DACE's extra operations can cause numbers to be different to 2 ULP. This function makes values 'close' to a number equal to that number.
   @param[in] DA object to test.
   @param[out] Rounded object
*/
DACE::DA roundDA(const DACE::DA& in, double tol=1e-14)
{
    if ( DACE::cons(in) - 1 < tol && DACE::cons(in) - 1 > -tol )
    {
        return in - tol;
    }
    return in;
}

#endif //ACROBAT_LINEARALGEBRA_H
