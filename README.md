# Orbital element conversions using differential algebra

This repository contains routines for converting to-and-from Keplerian orbital elements and cartesian states using differential algebra (DA). Such an approach allows one to compute the derivatives of the final state with respect to any of the original orbital elements to machine precision, and without the need for analytical expressions.

Two routines are provided, and are based upon the methods provided in \[1\]. The first, `par2ic` converts the orbital elements to its equivalent cartesian position and velocity vectors. The orbital elements must be expressed relative to an inertial frame. The second, `ic2par`, converts a cartesian state to its equivalent orbital elements.

These routines are designed for use with the Differential Algebra Computational Engine (DACE) library \[2\].

## References and Acknowledgements

\[1\]: Battin, R.H.: An Introduction to the Mathematics and Methods of Astrodynamics, p. 125. AIAA, New York (1987).

\[2\]: <a href="https://github.com/dacelib/dace.git">Differential Algebra Computational Engine on GitHub.</a>

