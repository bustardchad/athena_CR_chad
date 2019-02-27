//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ideal.cpp
//  \brief implements ideal EOS in general EOS framework, mostly for debuging
//======================================================================================

// C headers

// C++ headers
#include <algorithm>
#include <cmath>   // sqrt()
#include <fstream>
#include <iostream> // ifstream
#include <limits>
#include <sstream>
#include <stdexcept> // std::invalid_argument

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../parameter_input.hpp"
#include "../eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::RiemannAsq(Real rho, Real hint)
//  \brief Return adiabatic sound speed squared for use in Riemann solver.
Real EquationOfState::RiemannAsq(Real rho, Real hint) {
  return (gamma_ - 1.) * hint;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimplePres(Real rho, Real egas)
//  \brief Return gas pressure
Real EquationOfState::SimplePres(Real rho, Real egas) {
  return (gamma_ - 1.) * egas;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimpleEgas(Real rho, Real pres)
//  \brief Return internal energy density
Real EquationOfState::SimpleEgas(Real rho, Real pres) {
  return pres / (gamma_ - 1.);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::SimpleAsq(Real rho, Real pres)
//  \brief Return adiabatic sound speed squared
Real EquationOfState::SimpleAsq(Real rho, Real pres) {
  return gamma_ * pres / rho;
}
