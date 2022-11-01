//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"


//======================================================================================
/*! \file coolpulse.cpp
 *  \brief Inject CRs in a pulse at the left boundary, toward a cool cloud
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief CR injection with radiative cooling
//======================================================================================

static Real sigma=1.e30; // some large number to prevent diffusion (large enough?)

// prototypes for user-defined diffusion, source terms, and boundary conditions
//void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
//        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt);

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

/*void RadCool(MeshBlock *pmb, const Real time, const Real dt, 
     const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
     AthenaArray<Real> &cons);
*/

//void CRFluxSrcInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
//     int is, int ie, int js, int je, int ks, int ke);

//void CRReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
//     int is, int ie, int js, int je, int ks, int ke);

//void CROutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
//     int is, int ie, int js, int je, int ks, int ke);

void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
     int ks, int ke, int ngh);

void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);


void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 
    
}

// enroll user-defined boundary conditions and source terms
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserBoundaryFunction(inner_x1, FixMHDLeft);
  if(CR_ENABLED)
    EnrollUserCRBoundaryFunction(inner_x1, FixCRsourceLeft);
 // if(CR_ENABLED){
  //  EnrollUserBoundaryFunction(INNER_X1, CRFluxSrcInnerX1);
  //  EnrollUserBoundaryFunction(OUTER_X1, CROutflowOuterX1);
 //   EnrollUserBoundaryFunction(inner_x1, CRFluxSrcInnerX1);
 //   EnrollUserBoundaryFunction(outer_x1, CROutflowOuterX1);
 // }
 // EnrollUserExplicitSourceFunction(RadCool);
}

// enroll user-defined CR diffusion function
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
    if(CR_ENABLED){
    pcr->EnrollOpacityFunction(Diffusion);
  }
 // if(CR_ENABLED){
 //   pcr->EnrollDiffFunction(Diffusion);
 // }

}


// problem setup
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real gamma = peos->GetGamma();

  // get parameters from input file
 // Real pgas  = pin->GetOrAddReal("problem","pgas", 1.);         // gas pressure
 // Real rho_c = pin->GetOrAddReal("problem","rhocold",  1.);     // gas density inside cloud
 // Real rho_h = pin->GetOrAddReal("problem","rhohot", .01);      // gas density outside cloud
 // Real delta_r = pin->GetOrAddReal("problem","deltar", .01);    // thickness of cloud interface
 // Real r_cloud = pin->GetOrAddReal("problem","rcloud", .1);     // radius of cloud
 // Real x_cloud = pin->GetOrAddReal("problem","xcloud", 1.);     // x-coordinate of cloud center
 // Real y_cloud = pin->GetOrAddReal("problem","ycloud", 0.);     // y-coordinate of cloud center

  Real pgas = .00494;         // gas pressure
  Real rho_c = 34.31;     // gas density inside cloud
  Real rho_h = 0.3066;      // gas density outside cloud
  Real delta_r = .00025;    // thickness of cloud interface
  Real r_cloud = .05;     // radius of cloud
  Real x_cloud = 1.0;     // x-coordinate of cloud center
  Real y_cloud = 0.0;     // y-coordinate of cloud center
  Real z_cloud = 0.0;     // y-coordinate of cloud center

  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        Real x1 = pcoord->x1v(i) - x_cloud;
        Real x2 = pcoord->x2v(j) - y_cloud;
        Real x3 = pcoord->x3v(k) - z_cloud;
        Real r2 = x1*x1 + x2*x2 + x3*x3;
        Real r  = sqrt(r2);
        Real density = rho_h + (rho_c - rho_h) * 0.5 *
                               (1.0 - 1.0*tanh((r-r_cloud)/delta_r));

        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = 0.0;                             // no initial momentum
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){
          phydro->u(IEN,k,j,i) = pgas/(gamma-1.0);
        }
        
      }// end i
    }
  }

  if(CR_ENABLED){
    // set basic CR quantities
   // Real ecr0 = pin->GetOrAddReal("problem","ecr0", 1.e-10);  // initial CR energy density (a small number)
    Real ecr0 = 1.e-10;  // initial CR energy density (a small number)
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
           pcr->u_cr(CRE,k,j,i) = ecr0;
           pcr->u_cr(CRF1,k,j,i) = 0.0;                       // no initial CR flux
           pcr->u_cr(CRF2,k,j,i) = 0.0;
           pcr->u_cr(CRF3,k,j,i) = 0.0;
 //          pcr->q_cr(k,j,i) = 0.0;                            // no explicit CR source
         }
       }
     }

    // set opactiy sigma in the ghost zones
    int nz1 = block_size.nx1 + 2*(NGHOST);
    int nz2 = block_size.nx2;
    if(nz2 > 1) nz2 += 2*(NGHOST);
    int nz3 = block_size.nx3;
    if(nz3 > 1) nz3 += 2*(NGHOST);
    for(int k=0; k<nz3; ++k){
      for(int j=0; j<nz2; ++j){
        for(int i=0; i<nz1; ++i){
          pcr->sigma_diff(0,k,j,i) = sigma;
          pcr->sigma_diff(1,k,j,i) = sigma;
          pcr->sigma_diff(2,k,j,i) = sigma;
        }
      }
    }// end k,j,i

  }// End CR

  // set initial magnetic field (default: uniform horizontal)
  if(MAGNETIC_FIELDS_ENABLED){

   // Real bx = pin->GetOrAddReal("problem","bx", 1.);
   // Real by = pin->GetOrAddReal("problem","by", 0.);
   // Real bz = pin->GetOrAddReal("problem","bz", 0.);
    Real bx = .0349;
    Real by = 0.;
    Real bz = 0.;

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = bx;
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = by;
          }
        }
      }

    }

    if(block_size.nx3 > 1){

      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = bz;
          }
        }
      }
    }// end nx3

    // set cell centerd magnetic field
    // add magnetic energy density to the total energy
    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);

    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));
      
        }
      }
    }

  }// end MHD
  
  
  return;
}


// user-defined CR diffusion function for streaming
// Currently this function does the same thing as the default CR diffusion function,
// which is defined by DefaultDiff in /src/cr/cr.cpp, so there's not much point in enrolling it,
// but this is a useful place to experiment if you need to.
//void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
//        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt)
//{ 

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc)
{

  // set the default opacity to be a large value in the default hydro case
  CosmicRay *pcr=pmb->pcr;
  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is-1, iu=pmb->ie+1;
  if(pmb->block_size.nx2 > 1){
    jl -= 1;
    ju += 1;
  }
  if(pmb->block_size.nx3 > 1){
    kl -= 1;
    ku += 1;
  }

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){

        pcr->sigma_diff(0,k,j,i) = sigma;
        pcr->sigma_diff(1,k,j,i) = sigma;
        pcr->sigma_diff(2,k,j,i) = sigma;  

      }
    }
  }

  Real invlim=1.0/pcr->vmax;

  // The information stored in the array
  // b_angle is
  // b_angle[0]=sin_theta_b
  // b_angle[1]=cos_theta_b
  // b_angle[2]=sin_phi_b
  // b_angle[3]=cos_phi_b


  if(MAGNETIC_FIELDS_ENABLED){
    // First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
        // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                         + pcr->cwidth(i);
         // Real dprdx=(pcr->prtensor_cr(PC11,k,j,i+1) * u_cr(CRE,k,j,i+1)
         //              - pcr->prtensor_cr(PC11,k,j,i-1) * u_cr(CRE,k,j,i-1));
          Real dprdx=(u_cr(CRE,k,j,i+1) - u_cr(CRE,k,j,i-1))/3.0;
          dprdx /= distance;
          pcr->b_grad_pc(k,j,i) = bcc(IB1,k,j,i) * dprdx;
        }
        //y component
        pmb->pcoord->CenterWidth2(k,j-1,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth2(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth2(k,j+1,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                         + pcr->cwidth(i);
         // Real dprdy=(pcr->prtensor_cr(PC22,k,j+1,i) * u_cr(CRE,k,j+1,i)
         //                  - pcr->prtensor_cr(PC22,k,j-1,i) * u_cr(CRE,k,j-1,i));
          Real dprdy=(u_cr(CRE,k,j+1,i) - u_cr(CRE,k,j-1,i))/3.0;
	  dprdy /= distance;
          pcr->b_grad_pc(k,j,i) += bcc(IB2,k,j,i) * dprdy;

        }
        // z component
        pmb->pcoord->CenterWidth3(k-1,j,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth3(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth3(k+1,j,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                          + pcr->cwidth(i);
         // Real dprdz=(pcr->prtensor_cr(PC33,k+1,j,i) * u_cr(CRE,k+1,j,i)
         //                  - pcr->prtensor_cr(PC33,k-1,j,i) * u_cr(CRE,k-1,j,i));
          Real dprdz=(u_cr(CRE,k+1,j,i) - u_cr(CRE,k-1,j,i))/3.0;
	  dprdz /= distance;
          pcr->b_grad_pc(k,j,i) += bcc(IB3,k,j,i) * dprdz;

        }

      // now calculate the streaming velocity
      // streaming velocity is calculated with respect to the current coordinate system
      // diffusion coefficient is calculated with respect to B direction
        for(int i=il; i<=iu; ++i){
          Real pb= bcc(IB1,k,j,i)*bcc(IB1,k,j,i)
                  +bcc(IB2,k,j,i)*bcc(IB2,k,j,i)
                  +bcc(IB3,k,j,i)*bcc(IB3,k,j,i);
          Real inv_sqrt_rho = 1.0/sqrt(prim(IDN,k,j,i));
          Real va1 = bcc(IB1,k,j,i)*inv_sqrt_rho;
          Real va2 = bcc(IB2,k,j,i)*inv_sqrt_rho;
          Real va3 = bcc(IB3,k,j,i)*inv_sqrt_rho;

          Real va = sqrt(pb/prim(IDN,k,j,i));

          Real dpc_sign = 0.0;
          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = 1.0;
          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = -1.0;
          
          pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
          pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
          pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

          // now the diffusion coefficient

          if(va < TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
          }else{
            pcr->sigma_adv(0,k,j,i) = fabs(pcr->b_grad_pc(k,j,i))
                                   /(va * (1.0 + 1.0/3.0) 
                                    * invlim * u_cr(CRE,k,j,i) * sqrt(pb)); 
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;  

          // now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(pb);
          if(btot > TINY_NUMBER){
            pcr->b_angle(0,k,j,i) = bxby/btot;
            pcr->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            pcr->b_angle(0,k,j,i) = 1.0;
            pcr->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            pcr->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            pcr->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            pcr->b_angle(2,k,j,i) = 0.0;
            pcr->b_angle(3,k,j,i) = 1.0;            
          }

        }
      }
    }// end i,j,k

  }// end MHD  
  else{

    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
        // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
           Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                          + pcr->cwidth(i);
          // Real grad_pr=(pcr->prtensor_cr(PC11,k,j,i+1) * u_cr(CRE,k,j,i+1)
          //             - pcr->prtensor_cr(PC11,k,j,i-1) * u_cr(CRE,k,j,i-1));
           Real grad_pr=(u_cr(CRE,k,j,i+1) - u_cr(CRE,k,j,i-1))/3.0;
           grad_pr /= distance;

           Real va = 1.0;

           if(va < TINY_NUMBER){
             pcr->sigma_adv(0,k,j,i) = sigma;
             pcr->v_adv(0,k,j,i) = 0.0;
           }else{
             Real sigma2 = fabs(grad_pr)/(va * (1.0 + 1.0/3.0) 
                               * invlim * u_cr(CRE,k,j,i)); 
             if(fabs(grad_pr) < TINY_NUMBER){
               pcr->sigma_adv(0,k,j,i) = 0.0;
               pcr->v_adv(0,k,j,i) = 0.0;
             }else{
               pcr->sigma_adv(0,k,j,i) = sigma2;
               pcr->v_adv(0,k,j,i) = -va * grad_pr/fabs(grad_pr);     
             }
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
       
          pcr->v_adv(1,k,j,i) = 0.0;
          pcr->v_adv(2,k,j,i) = 0.0;

        }
      }
    }// end i,j,k

  }
}

// user-defined boundary condition
// outflow boundary condition for hydro variables
// reflecting boundary condition for CRs
// do this by reflecting CR flux in the desired direction
// I'm not using this BC anymore, but it may be useful

// user-defined boundary condition
// outflow boundary condition for hydro variables
// constant inward CR flux
// This controls the amount of CRs flowing into the domain through the boundary
// according to input parameters. It replaces the old method of injection which used
// an explicit CR source term
//void CRFluxSrcInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
//     int is, int ie, int js, int je, int ks, int ke)
//{

 void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{


  if(CR_ENABLED){
   // CosmicRay *pcr=pmb->pcr;
  //  Real ramptime = 3.0;
  //  Real time_factor = (1.-exp(-time/ramptime))*(1.-exp((time-30.0)/ramptime)); // time-dependence of CR source
  //  time_factor = std::max(time_factor, 0.);
  //  Real xflux = 2.4e-5*time_factor;
    Real ramptime = pcr->fluxdt;
    Real time_factor = (1.-exp(-time/ramptime))*(1.-exp((time-pcr->flux_t)/ramptime)); // time-dependence of CR source
    time_factor = std::max(time_factor, 0.);
    Real xflux = pcr->influx*time_factor;

    for (int n=0; n<(NCR); ++n) {
      if (n==(CRF1)) {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRF1,k,j,is-i) = xflux;  // set x-crflux to set value
          }
        }}
      }
      else {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(n,k,j,is-i) = u_cr(n,k,j,is);  // copy everything else as for hydro variables
          }
        }}
      }
    }
  }


}

/*
void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{

  Real fix_u = 3.0;


  if(CR_ENABLED){

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          u_cr(CRE,k,j,is-i) = fix_u;
          u_cr(CRF1,k,j,is-i) = u_cr(CRF1,k,j,is);
          u_cr(CRF2,k,j,is-i) = u_cr(CRF2,k,j,is);
          u_cr(CRF3,k,j,is-i) = u_cr(CRF3,k,j,is);

        }
      }
    }


  }
}
*/

void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
     int ks, int ke, int ngh)
{

  // copy hydro variables into ghost zones, reflecting v1

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        prim(IDN,k,j,is-i) = prim(IDN,k,j,is);
        prim(IVX,k,j,is-i) = -prim(IVX,k,j,is); // reflect 1-velocity
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        if(NON_BAROTROPIC_EOS)
          prim(IEN,k,j,is-i) = prim(IEN,k,j,is);
//          prim(IEN,k,j,is-i) = 10.0;
//         prim(IEN,k,j,is-i) = prim(IEN,k,j,is+i-1);
      }
    }
  }



  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
//        b.x1f(k,j,(is-i)) = sqrt(2.0*const_pb);  // reflect 1-field
          b.x1f(k,j,(is-i)) =  b.x1f(k,j,is);
      }
    }}
    if(je > js){
     for (int k=ks; k<=ke; ++k) {
     for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,is);
      }
     }}
    }
    if(ke > ks){
     for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
       for (int i=1; i<=(NGHOST); ++i) {
         b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
       }
      }}
    }
  }


}



// user-defined boundary condition
// outflow boundary condition for hydro variables
// outflow boundary condition for CRs
// do this by keeping P_CR in outer ghost zones slightly less than in domain
// This is a trick carried over from my ZEUS simulations - CRs won't stream out of the domain
// unless there is a slight CR gradient at the boundary.
// This is currently NOT used on the y-boundaries.
// (I'm not positive this BC is necessary - CRs might be able to leave normally with the default
// outflow BC)

/*
// user-defined source function
// analytic approximation to radiative cooling at solar metallicity
void RadCool(MeshBlock *pmb, const Real time, const Real dt, 
     const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
     AthenaArray<Real> &cons)
{
  // some constants
  Real mp  = 1.6726e-24;    // proton mass
  Real mu  = 0.6;           // mean molecular weight (fully ionized)
  Real kb  = 1.38065e-16;   // boltzmann constant
  Real Tfloor = 1.e4;       // temperature floor (no radiative cooling below this temperature)
  Real gHeat  = 1.e-25;     // heat constant (heating rate = n_H*gHeat)

  // some unit conversion factors - (add functionality later:  get these from param file)
  // current code units are:
  // length: 1 kpc
  // time: 1 Myr
  // mass: 10^5 M_sun
  Real dens2cgs = 6.85e-27;
  Real pres2cgs = 6.54e-11;
  Real heat2code = 4.81e23;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real dcgs = prim(IDN,k,j,i)*dens2cgs;  // density in cgs units
        Real pcgs = prim(IEN,k,j,i)*pres2cgs;  // gas pressure in cgs units
                                               // not that the variable IEN in the primitives vector is PRESSURE, not energy density
                                               // so no gamma factor is needed

        Real temp = mu*mp*pcgs/(dcgs*kb);      // temperature in Kelvin
        Real nH = 0.7*dcgs/mp;                 // proton density for fully ionized medium with ~25% Helium by mass

        Real tx = log10(temp)-5.;
        Real theta = 0.4*tx - 3. +5.2/(exp(tx+.08)+exp(-1.5*(tx+.08)));
        Real lambda = 1.1e-21*pow(10.,theta);
        if (temp < Tfloor) { lambda = 0.; }
        cons(IEN,k,j,i) -= dt*(nH*nH*lambda - nH*gHeat)*heat2code;  // apply cooling
      }
    }
  } // end i,j,k
  return;
}
*/
