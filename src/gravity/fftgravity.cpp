//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fftgravity.cpp
//  \brief implementation of functions in class FFTGravity

// Athena++ headers
#include "fftgravity.hpp"
#include "gravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../fft/athena_fft.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../task_list/grav_task_list.hpp"

#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>


//----------------------------------------------------------------------------------------
//! \fn FFTGravityDriver::FFTGravityDriver(Mesh *pm, ParameterInput *pin)
//  \brief FFTGravityDriver constructor

FFTGravityDriver::FFTGravityDriver(Mesh *pm, ParameterInput *pin)
 : FFTDriver(pm, pin)
{
  four_pi_G_=pmy_mesh_->four_pi_G_;
  if(four_pi_G_==0.0) {
   std::stringstream msg;
   msg << "### FATAL ERROR in GravityDriver::GravityDriver" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  int igid=Globals::my_rank;
  pmy_fb=new FFTGravity(this, fft_loclist_[igid], igid, fft_mesh_size_, fft_block_size_);

  pmy_fb->SetNormFactor(four_pi_G_/gcnt_);

  QuickCreatePlan();

  gtlist_ = new GravitySolverTaskList(pin, pm);
}


//----------------------------------------------------------------------------------------
//! \fn void GravityDriver::Solve(int step)
//  \brief load the data and solve

void FFTGravityDriver::Solve(int step)
{
  FFTBlock *pfb=pmy_fb;
  AthenaArray<Real> in;
  // Load the source 
  int nbs=nslist_[Globals::my_rank];
  int nbe=nbs+nblist_[Globals::my_rank]-1;
  for(int igid=nbs;igid<=nbe;igid++){
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if(pmb!=NULL) {
      if(step == 1) in.InitWithShallowSlice(pmb->phydro->u,4,IDN,1);
      else if(step == 2) in.InitWithShallowSlice(pmb->phydro->u1,4,IDN,1);
      pfb->LoadSource(in, 1, NGHOST, pmb->loc, pmb->block_size);
    }
//    else { // on another process
//    }
  }

  pfb->ExecuteForward();
  pfb->ApplyKernel(0);
  pfb->ExecuteBackward();

  // Return the result
  for(int igid=nbs;igid<=nbe;igid++){
    MeshBlock *pmb=pmy_mesh_->FindMeshBlock(igid);
    if(pmb!=NULL) {
      pfb->RetrieveResult(pmb->pgrav->phi, 1, NGHOST, 
                          pmb->loc, pmb->block_size);
    }
//    else { // on another process
//    }
  }

  gtlist_->DoTaskListOneSubstep(pmy_mesh_,step);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FFTGravity::ApplyKernel(const AthenaArray<Real> &src, int ns)
//  \brief Apply kernel
void FFTGravity::ApplyKernel(int mode)
{
  Real pcoeff;
  Real dx1sq=rdx_*rdx_;
  Real dx2sq=rdy_*rdy_;
  Real dx3sq=rdz_*rdz_;
  for(int k=0; k<knx_[2]; k++) {
    for(int j=0; j<knx_[1]; j++) {
      for(int i=0; i<knx_[0]; i++) {
        long int gidx = GetGlobalIndex(i,j,k);
        if(gidx == 0){ pcoeff = 0.0;}
        else {
          if(mode == 0){ // Discrete FT
            pcoeff = ((2.0*std::cos((i+kdisp_[0])*dkx_[0])-2.0)/dx1sq);
            if(dim_ > 1)
              pcoeff += ((2.0*std::cos((j+kdisp_[1])*dkx_[1])-2.0)/dx2sq);
            if(dim_ > 2)
              pcoeff += ((2.0*std::cos((k+kdisp_[2])*dkx_[2])-2.0)/dx3sq);
          } else if(mode == 1) { // Continous FT
            Real kx=(i+kdisp_[0])*dkx_[0];
            Real ky=(j+kdisp_[1])*dkx_[1];
            Real kz=(k+kdisp_[2])*dkx_[2];
            if(kx > PI) kx -= 2*PI;
            if(ky > PI) ky -= 2*PI;
            if(kz > PI) kz -= 2*PI;
            pcoeff = -(kx*kx/dx1sq+ky*ky/dx2sq+kz*kz/dx3sq);
          }
          pcoeff = 1.0/pcoeff;
        }

        long int idx_in=GetIndex(i,j,k,b_in_);
        long int idx_out=GetIndex(i,j,k,f_out_);
        in_[idx_in][0] = out_[idx_out][0]*pcoeff;
        in_[idx_in][1] = out_[idx_out][1]*pcoeff;
      }
    }
  }
  return;
}