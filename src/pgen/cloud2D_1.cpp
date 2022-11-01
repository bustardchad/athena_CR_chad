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
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================

static Real sigma=1.e8;

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt);

void CRFluxSrcInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);

void CRReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);

void CROutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 
    
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{

  EnrollUserBoundaryFunction(INNER_X1, CRFluxSrcInnerX1);
  EnrollUserBoundaryFunction(OUTER_X1, CROutflowOuterX1);

}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  if(CR_ENABLED){
//    pcr->EnrollDiffFunction(Diffusion);
  }

}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real gamma = peos->GetGamma();

  Real vx=0.0;
  Real pgas  = pin->GetReal("problem","pgas");
  Real rho_c = pin->GetReal("problem","rhocold");
  Real rho_h = pin->GetReal("problem","rhohot");
  Real delta_r = pin->GetReal("problem","deltar");
  Real r_cloud = pin->GetReal("problem","rcloud");
  Real x_cloud = pin->GetReal("problem","xcloud");
  Real y_cloud = pin->GetReal("problem","ycloud");

  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        Real x1 = pcoord->x1v(i) - x_cloud;
        Real x2 = pcoord->x2v(j) - y_cloud;
        Real r2 = x1*x1 + x2*x2;
        Real r  = sqrt(r2);
        Real density = rho_h + (rho_c - rho_h) * 0.5 *
                               (1.0 - 1.0*tanh((r-r_cloud)/delta_r));

        phydro->u(IDN,k,j,i) = density;

        phydro->u(IM1,k,j,i) = vx;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){
          phydro->u(IEN,k,j,i) = 0.5*vx*vx+pgas/(gamma-1.0);
        }
        
        if(CR_ENABLED){
            Real ecr0 = pin->GetReal("problem","ecr0");
            Real q0 = pin->GetReal("problem","q0");
            Real xleft = pin->GetReal("problem","xleft");
            Real xright = pin->GetReal("problem","xright");
            Real yleft = pin->GetReal("problem","yleft");
            Real yright = pin->GetReal("problem","yright");

            Real x1left = pcoord->x1f(i);
            Real x1right= pcoord->x1f(i+1);
            Real x2left = pcoord->x1f(j);
            Real x2right= pcoord->x1f(j+1);

            pcr->u_cr(CRE,k,j,i) = ecr0;
            pcr->u_cr(CRF1,k,j,i) = 0.0;
            pcr->u_cr(CRF2,k,j,i) = 0.0;
            pcr->u_cr(CRF3,k,j,i) = 0.0;

            // CR source function Q
            // set a constant CR source in a rectangular region defined by input parameters xleft/right, yleft/right

            // simplified version - ignore edge effects
            //pcr->q_cr(k,j,i) = 0.0;
            //if((xleft < x1) && (x1 < xright)){
            //  if((yleft < x2) && (x2 < yright)){
            //    pcr->q_cr(k,j,i) = q0;
            //  }
            //}

            // better version - account for volume factors when a cell is only partially inside the injection region
            pcr->q_cr(k,j,i) = q0;
            if(x1left < xleft){
              if(x1right < xleft){
                //these cells lie to the left of the injection region
                pcr->q_cr(k,j,i) = 0.0;
              }
              else if(x1right < xright){
                //these cells contain the left edge of the injection region
                pcr->q_cr(k,j,i) *= (x1right - xleft)/(x1right - x1left);
              }
              else {
                //these cells contain the entire x-extent of the injection region
                pcr->q_cr(k,j,i) *= (xright - xleft)/(x1right - x1left);
              }
            }
            else if(x1left < xright){
              if(x1right > xright){
                //these cells contain the right edge of the injection region
                pcr->q_cr(k,j,i) *= (xright - x1left)/(x1right - x1left);
              }
            }
            else {
              //these cells lie to the right of the injection region
              pcr->q_cr(k,j,i) = 0.0;
            }

            if(x2left < yleft){
              if(x2right < yleft){
                //these cells lie below the injection region
                pcr->q_cr(k,j,i) = 0.0;
              }
              else if(x2right < yright){
                //these cells contain the bottom edge of the injection region
                pcr->q_cr(k,j,i) *= (x2right - yleft)/(x2right - x2left);
              }
              else {
                //these cells contain the entire y-extent of the injection region
                pcr->q_cr(k,j,i) *= (yright - yleft)/(x2right - x2left);
              }
            }
            else if(x2left < yright){
              if(x2right > yright){
                //these cells contain the top edge of the injection region
                pcr->q_cr(k,j,i) *= (yright - x2left)/(x2right - x2left);
              }
            }
            else {
              //these cells lie above the injection region
              pcr->q_cr(k,j,i) = 0.0;
            }
        }
      }// end i
    }
  }
  //Need to set opactiy sigma in the ghost zones
  if(CR_ENABLED){

  // Default values are 1/3
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

    // Add horizontal magnetic field lines, to show streaming and diffusion 
  // along magnetic field ines
  if(MAGNETIC_FIELDS_ENABLED){

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = 1.0;
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }

    }

    if(block_size.nx3 > 1){

      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }// end nx3

    // set cell centerd magnetic field
    // Add magnetic energy density to the total energy
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



void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt)
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
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
    // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                         + pcr->cwidth(i);
          Real dprdx=(pcr->prtensor_cr(PC11,k,j,i+1) * u_cr(CRE,k,j,i+1)
                       - pcr->prtensor_cr(PC11,k,j,i-1) * u_cr(CRE,k,j,i-1));
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
          Real dprdy=(pcr->prtensor_cr(PC22,k,j+1,i) * u_cr(CRE,k,j+1,i)
                           - pcr->prtensor_cr(PC22,k,j-1,i) * u_cr(CRE,k,j-1,i));
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
          Real dprdz=(pcr->prtensor_cr(PC33,k+1,j,i) * u_cr(CRE,k+1,j,i)
                           - pcr->prtensor_cr(PC33,k-1,j,i) * u_cr(CRE,k-1,j,i));
          dprdz /= distance;
          pcr->b_grad_pc(k,j,i) += bcc(IB3,k,j,i) * dprdz;

          // now only get the sign
//          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) = 1.0;
//          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) 
//            = -1.0;
//          else pcr->b_grad_pc(k,j,i) = 0.0;
        }

      // now calculate the streaming velocity
      // streaming velocity is calculated with respect to the current coordinate 
      //  system
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
                                   /(va * (1.0 + pcr->prtensor_cr(PC11,k,j,i)) 
                                    * invlim * u_cr(CRE,k,j,i) * sqrt(pb)); 
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;  

          // Now calculate the angles of B
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

        }//        

      }// end j
    }// end k

  }// End MHD  
  else{



  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
  // x component
      pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
      for(int i=il; i<=iu; ++i){
         Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                        + pcr->cwidth(i);
         Real grad_pr=(pcr->prtensor_cr(PC11,k,j,i+1) * u_cr(CRE,k,j,i+1)
                     - pcr->prtensor_cr(PC11,k,j,i-1) * u_cr(CRE,k,j,i-1));
         grad_pr /= distance;

         Real va = 1.0;

         if(va < TINY_NUMBER){
           pcr->sigma_adv(0,k,j,i) = sigma;
           pcr->v_adv(0,k,j,i) = 0.0;
         }else{
           Real sigma2 = fabs(grad_pr)/(va * (1.0 + pcr->prtensor_cr(PC11,k,j,i)) 
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
  }

  }
}

// outflow boundary condition for hydro variables
// reflecting boundary condition for CRs
// do this by reflecting CR flux in the desired direction
void CRReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke)
{

  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        prim(n,k,j,is-i) = prim(n,k,j,is);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) { 
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
      } 
    }}
  
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,is);
      }
    }}  
        
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) =  b.x3f(k,j,is);
      }
    }}
  }


  if(CR_ENABLED){
    for (int n=0; n<(NCR); ++n) {
      if (n==(CRF1)) {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRF1,k,j,is-i) = -u_cr(CRF1,k,j,(is+i-1));  // mirror 1-crflux values across boundary, flipping sign
          }
        }}
      }
      else {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(n,k,j,is-i) = u_cr(n,k,j,(is+i-1));  // mirror everything else, not flipping sign
          }
        }}
      }
    }
  }


}

// outflow boundary condition for hydro variables
// constant inward CR flux
void CRFluxSrcInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke)
{

  CosmicRay *pcr=pmb->pcr;

  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        prim(n,k,j,is-i) = prim(n,k,j,is);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) { 
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) { 
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
      } 
    }}
  
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) =  b.x2f(k,j,is);
      }
    }}  
        
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) =  b.x3f(k,j,is);
      }
    }}
  }


  if(CR_ENABLED){
    for (int n=0; n<(NCR); ++n) {
      if (n==(CRF1)) {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(CRF1,k,j,is-i) = pcr->influx;  // set cr1-flux to set value
          }
        }}
      }
      else {
        for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
#pragma simd
          for (int i=1; i<=(NGHOST); ++i) {
            u_cr(n,k,j,is-i) = u_cr(n,k,j,(is));  // copy everything else as for hydro variables
          }
        }}
      }
    }
  }


}

// outflow boundary condition for hydro variables
// outflow boundary condition for CRs
// do this by keeping P_CR in outer ghost zones slightly less than in domain
void CROutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        prim(n,k,j,ie+i) = prim(n,k,j,ie);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
      }
    }}
  }

  if(CR_ENABLED){
    for (int n=0; n<(NCR); ++n) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          // enforce slight outward gradient to allow CRs to leave
          if(n==CRE) {
            u_cr(n,k,j,ie+i) = 0.999*u_cr(n,k,j,ie);
          } else {
            u_cr(n,k,j,ie+i) = u_cr(n,k,j,ie);
          }
        }
      }}
    }
  }  

}

