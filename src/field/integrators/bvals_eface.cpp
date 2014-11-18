//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Primary header
#include "../field.hpp"    // Field
#include "field_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()

// Athena headers
#include "../../athena.hpp"         // enums, macros, Real
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../mesh.hpp"           // MeshBlock
#include "../../fluid/fluid.hpp"    // Fluid

//======================================================================================
//! \file bvals_eface.cpp
//  \brief
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  void FieldIntegrator::BoundaryValuesFaceCenteredE3
//  \brief

void FieldIntegrator::BoundaryValuesFaceCenteredE3(MeshBlock *pmb)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e_x1f = pmb->pfield->e.x1f.ShallowCopy();
  AthenaArray<Real> e_x2f = pmb->pfield->e.x2f.ShallowCopy();
  AthenaArray<Real> w_x1f = pmb->pfield->wght.x1f.ShallowCopy();
  AthenaArray<Real> w_x2f = pmb->pfield->wght.x2f.ShallowCopy();

// boundary conditions for E3 on x2-face at inner/outer x1

  switch(pmb->block_bcs.ix1_bc){
    case 4:
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        e_x2f(X2E3,k,j,is-1) = e_x2f(X2E3,k,j,ie); 
        w_x2f(k,j,is-1) = w_x2f(k,j,ie);
      }}
      break;
    default:
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        e_x2f(X2E3,k,j,is-1) = e_x2f(X2E3,k,j,is); 
        w_x2f(k,j,is-1) = w_x2f(k,j,is);
      }}
  }

  switch(pmb->block_bcs.ox1_bc){
    case 4:
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        e_x2f(X2E3,k,j,ie+1) = e_x2f(X2E3,k,j,is); 
        w_x2f(k,j,ie+1) = w_x2f(k,j,is);
      }}
      break;
    default:
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        e_x2f(X2E3,k,j,ie+1) = e_x2f(X2E3,k,j,ie); 
        w_x2f(k,j,ie+1) = w_x2f(k,j,ie);
      }}
  }

// boundary conditions for E3 on x1-face at inner/outer x2

  switch(pmb->block_bcs.ix2_bc){
    case 4:
      for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie+1; ++i) {
        e_x1f(X1E3,k,js-1,i) = e_x1f(X1E3,k,je,i); 
        w_x1f(k,js-1,i) = w_x2f(k,je,i);
      }}
      break;
    default:
      for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie+1; ++i) {
        e_x1f(X1E3,k,js-1,i) = e_x1f(X1E3,k,js,i); 
        w_x1f(k,js-1,i) = w_x2f(k,js,i);
      }}
  }

  switch(pmb->block_bcs.ox2_bc){
    case 4:
      for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie+1; ++i) {
        e_x1f(X1E3,k,je+1,i) = e_x1f(X1E3,k,js,i); 
        w_x1f(k,je+1,i) = w_x2f(k,js,i);
      }}
      break;
    default:
      for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie+1; ++i) {
        e_x1f(X1E3,k,je+1,i) = e_x1f(X1E3,k,je,i); 
        w_x1f(k,je+1,i) = w_x2f(k,je,i);
      }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn  void FieldIntegrator::BoundaryValuesFaceCenteredE1
//  \brief

void FieldIntegrator::BoundaryValuesFaceCenteredE1(MeshBlock *pmb)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e_x2f = pmb->pfield->e.x2f.ShallowCopy();
  AthenaArray<Real> e_x3f = pmb->pfield->e.x3f.ShallowCopy();

  return;
}

//--------------------------------------------------------------------------------------
//! \fn  void FieldIntegrator::BoundaryValuesFaceCenteredE2
//  \brief

void FieldIntegrator::BoundaryValuesFaceCenteredE2(MeshBlock *pmb)
{
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> e_x1f = pmb->pfield->e.x1f.ShallowCopy();
  AthenaArray<Real> e_x3f = pmb->pfield->e.x3f.ShallowCopy();

  return;
}
