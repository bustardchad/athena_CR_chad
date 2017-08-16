#ifndef MGGRAVITY_HPP
#define MGGRAVITY_HPP

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mggravity.hpp
//  \brief defines MGGravity class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../multigrid/multigrid.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class Multigrid;

//! \class MGGravity
//  \brief Multigrid gravity solver for each block

class MGGravity : public Multigrid {
public:
  MGGravity(MultigridDriver *pmd, RegionSize isize, MGBoundaryFunc_t *MGBoundary,
            enum BoundaryFlag *input_bcs, bool root)
   : Multigrid(pmd,1,1,isize,MGBoundary,input_bcs,root), omega_(1.15)
  { btype=BNDRY_MGGRAV; btypef=BNDRY_MGGRAVF; };
  ~MGGravity() {};
  void Smooth(int color);
  void CalculateDefect(void);

private:
  const Real omega_;
};


//! \class GravityDriver
//  \brief Multigrid gravity solver

class GravityDriver : public MultigridDriver{
public:
  GravityDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary, ParameterInput *pin);
  void LoadSourceAndData(void);
  Multigrid* AllocateNewMultigrid(RegionSize isize, MGBoundaryFunc_t *MGBoundary,
                                  enum BoundaryFlag *input_bcs, bool root);

private:
  Real four_pi_G_;
};

#endif // MGGRAVITY_HPP