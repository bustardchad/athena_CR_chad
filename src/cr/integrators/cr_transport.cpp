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
//! \file rad_transport.cpp
//  \brief implementation of radiation integrators
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../reconstruct/reconstruction.hpp"
#include "../../utils/utils.hpp"
#include <algorithm>   // min,max

#include <math.h>

// class header
#include "cr_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void CRIntegrator::CalculateFluxes(AthenaArray<Real> &w,
            AthenaArray<Real> &bcc, AthenaArray<Real> &cr, const int order)
{
  CosmicRay *pcr=pmy_cr;
  MeshBlock *pmb=pcr->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Real invlim = 1.0/pcr->vmax;


  double neTable[25] = {0.0023,0.0142,0.0704,0.2534,0.5678,0.8177,0.9324, 
			0.9780,1.004,1.035,1.067,1.086,1.094,1.098,1.100,1.102,
		        1.110,1.128,1.154,1.176,1.188,1.195,1.198,1.199,1.200};	
  double tempTable[25] = {4.0,4.050,4.1,4.150,4.200,4.250,4.300,4.35,
			4.40,4.45,4.50,4.55,4.60,4.65,4.70,4.75,4.80,4.85,4.90,4.95,
			5.00,5.05,5.10,5.15,5.20}; 

  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, 
  ncells3 = pmb->ncells3; 

  AthenaArray<Real> &x1flux=pcr->flux[X1DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  jl = js, ju=je, kl=ks, ku=ke;

  if(ncells2 > 1)
  {
    if(ncells3 == 1){
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    }else{
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
 
  }

//--------------------------------------------------------------------------------------
  for (int k=0; k<ncells3; ++k){
    for(int j=0; j<ncells2; ++j){

      // diffusion velocity along the direction of sigma vector
      // We first assume B is along x coordinate
      // Then rotate according to B direction to the actual acooridnate

      for(int i=0; i<ncells1; ++i){
        Real eddxx=1.0/3.0;
        Real totsigma = pcr->sigma_diff(0,k,j,i);
        if(pcr->stream_flag > 0)
          totsigma = 1.0/(1.0/pcr->sigma_diff(0,k,j,i) 
                        + 1.0/pcr->sigma_adv(0,k,j,i));
        Real taux = taufact_ * totsigma * pco->dx1f(i);
        taux = taux * taux/(2.0 * eddxx);
        Real diffv = sqrt((1.0 - exp(-taux)) / taux);

        if(taux < 1.e-3)
          diffv = sqrt((1.0 - 0.5* taux));

        pcr->v_diff(0,k,j,i) = pcr->vmax * sqrt(eddxx) * diffv;
      }// end i direction

       // y direction
      if(ncells2 >1){
        pco->CenterWidth2(k,j,0,ncells1-1,cwidth2_);
        // get the optical depth across the cell
        for(int i=0; i<ncells1; ++i){
          Real eddyy=1.0/3.0;
          Real totsigma = pcr->sigma_diff(1,k,j,i);
          if(pcr->stream_flag > 0)
            totsigma = 1.0/(1.0/pcr->sigma_diff(1,k,j,i) 
                        + 1.0/pcr->sigma_adv(1,k,j,i));
          Real tauy = taufact_ * totsigma * cwidth2_(i);
          tauy = tauy * tauy/(2.0 * eddyy);
          Real diffv = sqrt((1.0 - exp(-tauy)) / tauy);

          if(tauy < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauy));

          pcr->v_diff(1,k,j,i) = pcr->vmax * sqrt(eddyy) * diffv;            
        }// end i
      }else{
        for(int i=0; i<ncells1; ++i)
          pcr->v_diff(1,k,j,i) = 0.0;
      }

     // z direction
      if(ncells3 > 1){
        pco->CenterWidth3(k,j,0,ncells1-1,cwidth3_);
        // get the optical depth across the cell
        for(int i=0; i<ncells1; ++i){
          Real eddzz=1.0/3.0;
          Real totsigma = pcr->sigma_diff(2,k,j,i);
          if(pcr->stream_flag > 0)
            totsigma = 1.0/(1.0/pcr->sigma_diff(2,k,j,i) 
                        + 1.0/pcr->sigma_adv(2,k,j,i));
          Real tauz = taufact_ * totsigma * cwidth3_(i);
          tauz = tauz * tauz/(2.0 * eddzz);
          Real diffv = sqrt((1.0 - exp(-tauz)) / tauz);

          if(tauz < 1.e-3)
            diffv = sqrt((1.0 - 0.5* tauz));

          pcr->v_diff(2,k,j,i) = pcr->vmax * sqrt(eddzz) * diffv;            
        }
      }else{
        for(int i=0; i<ncells1; ++i)
          pcr->v_diff(2,k,j,i) = 0.0;
      }

        //rotate the v_diff vector to the local coordinate
      if(MAGNETIC_FIELDS_ENABLED){
        for(int i=0; i<ncells1; ++i){

          InvRotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                        pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i), 
                          pcr->v_diff(0,k,j,i),pcr->v_diff(1,k,j,i),
                                              pcr->v_diff(2,k,j,i));
            // take the absolute value
            // Also add the Alfven velocity for the streaming flux
          pcr->v_diff(0,k,j,i) = fabs(pcr->v_diff(0,k,j,i));

          pcr->v_diff(1,k,j,i) = fabs(pcr->v_diff(1,k,j,i));
                                   
          pcr->v_diff(2,k,j,i) = fabs(pcr->v_diff(2,k,j,i));

        }

      }// end MHD

       // need to add additional sound speed for stability
      for(int i=0; i<ncells1; ++i){
         Real cr_sound_x = vel_flx_flag_ * sqrt((4.0/9.0) * cr(CRE,k,j,i)/w(IDN,k,j,i)); 

         pcr->v_diff(0,k,j,i) += cr_sound_x;

         if(ncells2 > 1){
           Real cr_sound_y = vel_flx_flag_ * sqrt((4.0/9.0) * cr(CRE,k,j,i)/w(IDN,k,j,i));

           pcr->v_diff(1,k,j,i) += cr_sound_y;

         }

         if(ncells3 > 1){
           Real cr_sound_z = vel_flx_flag_ * sqrt((4.0/9.0) * cr(CRE,k,j,i)/w(IDN,k,j,i)); 
         
           pcr->v_diff(2,k,j,i) += cr_sound_z;
         }
      }

    }// end j
  }// end k

  // prepare Array for reconstruction
  for(int n=0; n<NCR; ++n){
    for (int k=0; k<ncells3; ++k){
      for(int j=0; j<ncells2; ++j){
        for(int i=0; i<ncells1; ++i){
           ucr_vel_(n,k,j,i) = cr(n,k,j,i);
        }// end i
      }// end j
    }// end k
  }// end n

//--------------------------------------------------------------------------------------
// i-direction

  // add vx velocity
  for (int k=0; k<ncells3; ++k){
    for(int j=0; j<ncells2; ++j){
      for(int i=0; i<ncells1; ++i){
         ucr_vel_(NCR,k,j,i) = w(IVX,k,j,i);
      }// end i
    }// end j
  }// end k    

  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      // First, need to do reconstruction
      // to reconstruct Ec, Fc, vel, v_a and 
      // return Ec,Fc and signal speed at left and right state
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, ucr_vel_, ucr_l_, ucr_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, ucr_vel_, ucr_l_, ucr_r_);
      } else {
        pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, ucr_vel_, ucr_l_, ucr_r_);
      }

      // get the optical depth across the cell
#pragma omp simd
      for(int i=is; i<=ie+1; ++i){
        vdiff_l_(i) = pcr->v_diff(0,k,j,i-1);
        vdiff_r_(i) = pcr->v_diff(0,k,j,i);
      }

      // calculate the flux
      CRFlux(CRF1, is, ie+1, ucr_l_, ucr_r_, vdiff_l_, vdiff_r_, dflx_);
      // store the flux
      for(int n=0; n<NCR; ++n){
#pragma omp simd
        for(int i=is; i<=ie+1; ++i){
          x1flux(n,k,j,i) = dflx_(n,i);
        }
      }

    }
  }

//--------------------------------------------------------------------------------------
// j-direction
  if(pmb->pmy_mesh->f2){

    AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
  // add vy velocity
    for (int k=0; k<ncells3; ++k){
      for(int j=0; j<ncells2; ++j){
        for(int i=0; i<ncells1; ++i){
           ucr_vel_(NCR,k,j,i) = w(IVY,k,j,i);
        }// end i
      }// end j
    }// end k    


    il=is-1; iu=ie+1; kl=ks; ku=ke;
    if (ncells3 ==  1) // 2D
      kl = ks, ku = ke;
    else // 3D
      kl = ks-1, ku = ke+1;    
    for (int k=kl; k<=ku; ++k){
      //reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else {
        pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      }


      for (int j=js; j<=je+1; ++j){
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        }

        // get the optical depth across the cell
#pragma omp simd
        for(int i=il; i<=iu; ++i){
          vdiff_l_(i) = pcr->v_diff(1,k,j-1,i);
          vdiff_r_(i) = pcr->v_diff(1,k,j,i);
        }
      // calculate the flux
        CRFlux(CRF2, il, iu, ucr_l_, ucr_r_, vdiff_l_, vdiff_r_, dflx_);           
      // store the flux
        for(int n=0; n<NCR; ++n){
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            x2flux(n,k,j,i) = dflx_(n,i);
          }
        }
       // swap the array for next cycle
        ucr_l_.SwapAthenaArray(ucr_lb_);

      }// end j from js to je+1
    }
  }// finish j direction



//  k-direction
  if(pmb->pmy_mesh->f3){
    AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
    il =is-1, iu=ie+1, jl=js-1, ju=je+1;
  // add vz velocity
    for (int k=0; k<ncells3; ++k){
      for(int j=0; j<ncells2; ++j){
        for(int i=0; i<ncells1; ++i){
           ucr_vel_(NCR,k,j,i) = w(IVZ,k,j,i);
        }// end i
      }// end j
    }// end k   

    for(int j=jl; j<=ju; ++j){

      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      } else {
        pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, ucr_vel_, ucr_l_, ucr_r_);
      }

      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        } else {
          pmb->precon->PiecewiseParabolicX3(k, j, il, iu, ucr_vel_, ucr_lb_, ucr_r_);
        }

#pragma omp simd
        for(int i=il; i<=iu; ++i){
          vdiff_l_(i) = pcr->v_diff(2,k-1,j,i);
          vdiff_r_(i) = pcr->v_diff(2,k,j,i);
        }
      // calculate the flux
        CRFlux(CRF3, il, iu, ucr_l_, ucr_r_, vdiff_l_, vdiff_r_, dflx_);   
        for(int n=0; n<NCR; ++n){
#pragma omp simd
          for(int i=il; i<=iu; ++i){
            x3flux(n,k,j,i) = dflx_(n,i);
          }
        }

       // swap the array for next cycle
        ucr_l_.SwapAthenaArray(ucr_lb_);                

      }// end k loop
    }// end j loop 

  }// finish k direction

//-------------------------------------------------------------------------------------------------
// Now calculate Grad Pc and the associated heating term
// the flux divergence term is Grad P_c for \partial F_c/\partial t
  // only do this for the MHD case and along direction perpendicular 
  // to the magnetic field
  if(MAGNETIC_FIELDS_ENABLED){

    for (int k=ks; k<=ke; ++k) { 
      for (int j=js; j<=je; ++j) {
        pmb->pcoord->Face1Area(k,j,is,ie+1,x1face_area_);
        pmb->pcoord->CellVolume(k,j,is,ie,cell_volume_);
      // x1 direction
        for(int n=0; n<3; ++n){
  #pragma omp simd
          for(int i=is; i<=ie; ++i){
            grad_pc_(n,k,j,i) = (x1face_area_(i+1)*x1flux(CRF1+n,k,j,i+1) 
                               - x1face_area_(i)  *x1flux(CRF1+n,k,j,i))/cell_volume_(i);
          }
        } 

        if(pmb->block_size.nx2 > 1){
          AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
          pmb->pcoord->Face2Area(k,j  ,is,ie,x2face_area_   );
          pmb->pcoord->Face2Area(k,j+1,is,ie,x2face_area_p1_);
          for(int n=0; n<3; ++n){
  #pragma omp simd
            for(int i=is; i<=ie; ++i){
              grad_pc_(n,k,j,i) += (x2face_area_p1_(i)*x2flux(CRF1+n,k,j+1,i) 
                                 -  x2face_area_(i)  *x2flux(CRF1+n,k,j,i))/cell_volume_(i);
            }// end i
          }  
        }// end nx2

        if(pmb->block_size.nx3 > 1){
          AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
          pmb->pcoord->Face3Area(k  ,j,is,ie,x3face_area_);
          pmb->pcoord->Face3Area(k+1,j,is,ie,x3face_area_p1_);
          for(int n=0; n<3; ++n){
  #pragma omp simd
            for(int i=is; i<=ie; ++i){
              grad_pc_(n,k,j,i) += (x3face_area_p1_(i) *x3flux(CRF1+n,k+1,j,i) 
                                  - x3face_area_(i)*x3flux(CRF1+n,k,j,i))/cell_volume_(i);
            } 
          } 
        }// end nx3


        for(int n=0; n<3; ++n){
  #pragma omp simd
          for(int i=is; i<=ie; ++i){
            grad_pc_(n,k,j,i) *= invlim;
          }
        } 

        //need to subtract the coordinate source term to get the actual grad Pc for c
        // curlinear coordinate system
        pmb->pcoord->AddCoordTermsDivergence(cr, grad_pc_);

     // calculate streaming velocity with magnetic field
        for(int i=is; i<=ie; ++i){

          Real inv_sqrt_rho = 1.0/sqrt(w(IDN,k,j,i));

          Real pb= bcc(IB1,k,j,i)*bcc(IB1,k,j,i)
                  +bcc(IB2,k,j,i)*bcc(IB2,k,j,i)
                  +bcc(IB3,k,j,i)*bcc(IB3,k,j,i);

          Real b_grad_pc = bcc(IB1,k,j,i) * grad_pc_(0,k,j,i) 
                         + bcc(IB2,k,j,i) * grad_pc_(1,k,j,i) 
                         + bcc(IB3,k,j,i) * grad_pc_(2,k,j,i);

          // new to calculate pressure anisotropy contribution -- Chad, Jan 10, 2021
          // need magnetic length scale lb = |B/nabla B|
          //
          //

          Real lcgs = 3.0856e21; // 1 kpc is a length unit of 1? 
     //     pmb->pcoord->Face1Area(k,j,is,ie+1,x1face_area_);
     //     pmb->pcoord->CellVolume(k,j,is,ie,cell_volume_);
          Real nablaB1 = (x1face_area_(i+1)*bcc(IB1,k,j,i+1)- x1face_area_(i)*bcc(IB1,k,j,i))/cell_volume_(i);
         // std::cout << "nablaB1: " << nablaB1;
          Real nablaB2 = 0.0;
          if(pmb->block_size.nx2 > 1){
            nablaB2 = (x2face_area_p1_(i)*bcc(IB2,k,j+1,i)- x2face_area_(i)*bcc(IB2,k,j,i))/cell_volume_(i);
          } 
         // std::cout << "nablaB2: " << nablaB2;
          Real nablaB3 = 0.0;
	  if(pmb->block_size.nx3 > 1){
            nablaB3 = (x3face_area_p1_(i)*bcc(IB3,k+1,j,i)- x3face_area_(i)*bcc(IB3,k,j,i))/cell_volume_(i);
          }

          Real BnablaB = bcc(IB1,k,j,i)*nablaB1 + bcc(IB2,k,j,i)*nablaB2 + bcc(IB3,k,j,i)*nablaB3;
          Real lb = pb / std::max(TINY_NUMBER,sqrt(pow(BnablaB,2.)));  // magnetic length scale
	  lb *= lcgs;
          pcr->bscale(k,j,i)=lb;

	  Real lcr = 1./3.*cr(CRE,k,j,i) / std::max(TINY_NUMBER,sqrt(pow(b_grad_pc,2.)));
          lcr *= sqrt(pb);  // to get Pcr/|grad_parallel(Pcr)|
	  lcr *= lcgs;
          pcr->crscale(k,j,i)=lcr;
          /////////////////////////////////////////////////////////////////////////////

          Real va1 = bcc(IB1,k,j,i) * inv_sqrt_rho;
          Real va2 = bcc(IB2,k,j,i) * inv_sqrt_rho;
          Real va3 = bcc(IB3,k,j,i) * inv_sqrt_rho;

          Real va = sqrt(pb) * inv_sqrt_rho;
          Real dpc_sign = 0.0;

          if(b_grad_pc > TINY_NUMBER) dpc_sign = 1.0;
          else if(-b_grad_pc > TINY_NUMBER) dpc_sign = -1.0;

          pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
          pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
          pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

	  Real stream_boost_factor = 1.0;
        //  std::cout << "lb: " << lb << "lcr: " << lcr << " use_piso: " << pcr->use_piso;

          // is pressure anisotropy condition satisfied?
          Real piso_diff = 1.E-10;
          if ((float(lcr/lb) > (1.e2/va)) && (pcr->use_piso > 0)) {
            piso_diff = 1.e10*lb / 3.e29;
           // std::cout << "lb: " << lb << " piso_diff: " << piso_diff;
          }
	  pcr->sigma_diff(0,k,j,i) = 1./(pcr->diffInput + std::min(piso_diff, pcr->max_in_diff)); // if assuming fully ionized
	  pcr->uncapped_diffusion(k,j,i) = (piso_diff)*3.0e29; // in cm^2/s
          

	  if (pcr->ion_alfven_speed > 0) { // Chad, March 8, 2020
            Real fneutral = 1.2;
            Real dens2cgs = 6.85e-27;
            Real pres2cgs = 6.54e-11; 
            Real dcgs = w(IDN,k,j,i)*dens2cgs;  // density in cgs units
            Real pcgs = w(IEN,k,j,i)*pres2cgs;  // gas pressure in cgs units
	    Real temp = 1.67e-24*pcgs/(dcgs*1.38e-16);
	    Real nH = 0.6*dcgs/1.67e-24;
            // How do I calculate fneutral?
            // x = ne/nH from eqn 25 of Ferriere+ 1988
            // fneutral = (1.19e-4 - 1.36e-8 * pow(temp,0.845)/nH) + pow(1.42e-8 + 2.72e-8 * pow(temp,0.845)/nH + 1.85e-16 * pow(temp,1.69)/(nH*nH),0.5);
            //
            //
            /*
	    // from Sutherland-Dopita 1993 table:
	    int count = 0;
	    if (log10(temp) < tempTable[0]) {
		fneutral = std::max(1.e-4,neTable[0] + 
		(neTable[1]-neTable[0])/(tempTable[1]-tempTable[0]) * (log10(temp)-tempTable[0]));
	    }
	    else if (log10(temp) > tempTable[24]) {
		fneutral = 1.20;
	    }
	    else {
		while (count < 24) {
		    if ((log10(temp) < tempTable[count + 1]) && (log10(temp) > tempTable[count])) {
			fneutral = neTable[count] + (neTable[count+1]-neTable[count])/(tempTable[count+1]-tempTable[count]) * (log10(temp)-tempTable[count]);
		        break;
	    	    }
		    else {count += 1;}
		}
	    }
	    */
	//    if (temp < 2.e4) {
	 //     std::cout << "ne: " << fneutral;
          //  }
          //
            // New tanh function for fneutral: Chad, March 27, 2020
           // Real tUnit = 0.6*1.67e-24*pres2cgs/(dens2cgs*1.38e-16);
            Real tUnit = 1.67e-24*pres2cgs/(dens2cgs*1.38e-16);
            Real tempCutoff = 1.6e4;
	   // Real coeff = 3.e3;
	    Real coeff = 2.e3;
            fneutral = 0.5*(pcr->maxNeut)*(1.0 + tanh((-temp + tempCutoff)/(coeff)));

            Real fion = 1.0 - fneutral;
            fion = std::max(fion,1.e-8);
          //  fneutral = 1.0 - fneutral/1.20;
	    
            // Real temp = (1.377e14*1.67e-24/1.38e-16)*w(IDN,k,j,i)/w(IEN,k,j,i);
            // assumes the fiducial temp outside the cloud is 1e6 when P = 1/gamma, rho = 1, and 
            // I assume that rho = 1 corresponds to an actual rho of 1e-26
            //Real dampRate = 1.e-9 * fneutral * (temp / 
	    //		    1000.0)**0.5 * w(IDN,k,j,i)/1.e2; 
            Real dampRate = 1.e-9 * fneutral * sqrt(temp / 1000.0) * dcgs/1.e-24; 
           // stream_boost_factor = 1.0 + (4.0*3.0e10*dampRate*0.5*pb*pres2cgs
	   //			/(4.803e-10*pow(pb*pres2cgs,1.5)*2.0*(sqrt(8*3.1415)*3.1415/dcgs)*w(IDN,k,j,i)*cr(CRE,k,j,i)*pres2cgs));
            stream_boost_factor = 1.0/sqrt(fion);

	   // std::cout << "b_grad_pc: " << b_grad_pc << "\n";
	   // std::cout << "b_grad_pc: " << sqrt(pow(b_grad_pc,2.)) << "\n";
	   // Real lcr = 1./3.*w(IDN,k,j,i)*cr(CRE,k,j,i) / std::max(TINY_NUMBER,sqrt(pow(b_grad_pc,2.)));
	   // Real lcr = 1./3.*cr(CRE,k,j,i) / std::max(TINY_NUMBER,sqrt(pow(b_grad_pc,2.)));
           // lcr *= sqrt(pb);  // to get Pcr/|grad_parallel(Pcr)|
	   // lcr *= lcgs;

	    if (pcr->NLLD > 0) {
              dampRate = dampRate + 7.0e-12*sqrt(cr(CRE,k,j,i)*pres2cgs*624.15*6.2415e11*3.0856e21/lcr)*sqrt(sqrt((temp/10000.0)/(fion*dcgs/1.e-24)));
            }
	    Real damping_diff_term = (4.0*3.0e10*dampRate*0.5*pb*pres2cgs*fion*(4./3.)*lcr
				/(4.803e-10*pow(pb*pres2cgs,1.5)*2.0*(sqrt(8*3.1415)*3.1415/dcgs)*cr(CRE,k,j,i)*pres2cgs));
				///(4.803e-10*pow(pb*pres2cgs,1.5)*2.0*(sqrt(8*3.1415)*3.1415/dcgs)*w(IDN,k,j,i)*cr(CRE,k,j,i)*pres2cgs));
	    

	    damping_diff_term *= va/sqrt(fion);
            damping_diff_term *= 3.155e13/(lcgs*lcgs);  // back into code units now?

	    damping_diff_term *= 1./624.15; // because gamma_L is in GeV, and 1/624 = 1 GeV in ergs

			//	/(4.803e-10*pow(pb*pres2cgs,1.5)*2.0*(sqrt(8*3.1415)*3.1415/dcgs)*w(IDN,k,j,i)*cr(CRE,k,j,i)*pres2cgs));
			//	/(4.803e-10*sqrt(4*3.1415*pb*pres2cgs)*3.1415*(pow(va,2.0)*pres2cgs/dens2cgs)*w(IDN,k,j,i)*cr(CRE,k,j,i)*pres2cgs));
  	   // std::cout << "eb/ecr: " << 0.5*pb/cr(CRE,k,j,i);
  	   // std::cout << "stream boost factor: " << (1.0-stream_boost_factor)*2.0*cr(CRE,k,j,i)/pb;
           // stream_boost_factor = std::min(stream_boost_factor, 100.0);

	   /*
	   if (damping_diff_term > 10000.0) {
		std::cout << "||  density: " << dcgs << "  temperature: " << temp << "  eb: " << 0.5*pb << "  ecr: " << w(IDN,k,j,i)*cr(CRE,k,j,i);
		std::cout << "fion: " << fion << " lcr: " << lcr << "  va: " << va << "  damp rate: " << dampRate << "  diffusive term: " << damping_diff_term;
	   }
	   */

           // pcr->v_adv(0,k,j,i) = -(std::min(va1*stream_boost_factor,pcr->vmax/50.0)) * dpc_sign;
           // pcr->v_adv(1,k,j,i) = -(std::min(va2*stream_boost_factor,pcr->vmax/50.0)) * dpc_sign;
           // pcr->v_adv(2,k,j,i) = -(std::min(va3*stream_boost_factor,pcr->vmax/50.0)) * dpc_sign;
	    
           // va *= stream_boost_factor;  // Chad, March 28, 2020
            va1 *= stream_boost_factor;  // Chad, March 28, 2020
            va2 *= stream_boost_factor;  // Chad, March 28, 2020
            va3 *= stream_boost_factor;  // Chad, March 28, 2020

	    va1 = std::min(va1,pcr->max_ion_va);
	    va2 = std::min(va2,pcr->max_ion_va);
	    va3 = std::min(va3,pcr->max_ion_va);
	    va = sqrt(pow(va1,2.) + pow(va2,2.) + pow(va3,2.));

            
            if ((float(lcr/lb) > (1.e2/va)) && (pcr->use_piso > 0)) {piso_diff = 1.e10*lb / 3.e29;}
	    pcr->uncapped_diffusion(k,j,i) = (damping_diff_term + piso_diff)*3.0e29; // in cm^2/s
		
	   // std::cout << "||  density: " << dcgs << "  temperature: " << temp << "  eb: " << 0.5*pb << "  ecr: " << w(IDN,k,j,i)*cr(CRE,k,j,i);
	   // std::cout << "fion: " << fion << " lcr: " << lcr << "  va: " << va << "  damp rate: " << dampRate << "  diffusive term: " << damping_diff_term << "\n";
	    pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
            pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
            pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;
	   // pcr->v_adv(0,k,j,i) *= stream_boost_factor;
	   // pcr->v_adv(1,k,j,i) *= stream_boost_factor;
	   // pcr->v_adv(2,k,j,i) *= stream_boost_factor;
	    if (pcr->ion_neutral_friction > 0.) {
	      if (fneutral > TINY_NUMBER) {
	       // std::cout << "sigma_diff: " << 1./damping_diff_term << " k parallel: " << damping_diff_term << "\n";
	       // pcr->sigma_diff(0,k,j,i) = 1./(std::min(damping_diff_term, pcr->max_in_diff));  // Chad, March 28, 2020
	        pcr->sigma_diff(0,k,j,i) =  1./(pcr->diffInput + std::min(damping_diff_term + piso_diff, pcr->max_in_diff));  // Chad, March 28, 2020
	      }
	    }
	    else {
	      pcr->sigma_diff(0,k,j,i) = 1./(pcr->diffInput + std::min(piso_diff, pcr->max_in_diff));
	    }

	    /*
	    if (pcr->v_adv(0,k,j,i) < 0) {
		pcr->v_adv(0,k,j,i) = -std::min(std::abs(pcr->v_adv(0,k,j,i)),1.0);
	    }
	    else {pcr->v_adv(0,k,j,i) = std::min(std::abs(pcr->v_adv(0,k,j,i)),1.0);}
	    
	    if (pcr->v_adv(1,k,j,i) < 0) {
		pcr->v_adv(1,k,j,i) = -std::min(std::abs(pcr->v_adv(1,k,j,i)),1.0);
	    }
	    else {pcr->v_adv(1,k,j,i) = std::min(std::abs(pcr->v_adv(1,k,j,i)),1.0);}
	    
	    if (pcr->v_adv(2,k,j,i) < 0) {
		pcr->v_adv(2,k,j,i) = -std::min(std::abs(pcr->v_adv(2,k,j,i)),1.0);
	    }
	    else {pcr->v_adv(2,k,j,i) = std::min(std::abs(pcr->v_adv(2,k,j,i)),1.0);}
	    */
	   // pcr->v_adv(1,k,j,i) = std::min(std::abs(pcr->v_adv(1,k,j,i)),1.0);
	   // pcr->v_adv(2,k,j,i) = std::min(std::abs(pcr->v_adv(2,k,j,i)),1.0);

	  }
          
          if(va > TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = fabs(b_grad_pc)/(sqrt(pb) * va * 
                                   (4.0/3.0) * invlim * cr(CRE,k,j,i));
            pcr->dedt_damping(k,j,i) = va1*grad_pc_(0,k,j,i) + va2*grad_pc_(1,k,j,i) + va3*grad_pc_(2,k,j,i);
           /* if (pcr->ion_neutral_friction > 0) {
              pcr->sigma_adv(0,k,j,i) = fabs(b_grad_pc)/(sqrt(pb) * 
				   std::min(va * stream_boost_factor,1.0) * 
                                   (4.0/3.0) * invlim * cr(CRE,k,j,i));
            } */
            pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
            pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
          }


          if ((float(lcr/lb) > (1.e2/va)) && (pcr->use_piso > 0)) {
            pcr->sigma_adv(0,k,j,i) = pcr->max_opacity; 
            pcr->dedt_damping(k,j,i) = 0.0;
            pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
            pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
            pcr->v_adv(0,k,j,i) = TINY_NUMBER;  // these go into vtot, which is in the flux equation
            pcr->v_adv(1,k,j,i) = TINY_NUMBER;
            pcr->v_adv(2,k,j,i) = TINY_NUMBER;
          }

            

          Real v1 = w(IVX,k,j,i);
          Real v2 = w(IVY,k,j,i);
          Real v3 = w(IVZ,k,j,i);

          Real dpcdx = grad_pc_(0,k,j,i);
          Real dpcdy = grad_pc_(1,k,j,i);
          Real dpcdz = grad_pc_(2,k,j,i);

          


          RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                   pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),dpcdx,dpcdy,dpcdz);

          RotateVec(pcr->b_angle(0,k,j,i),pcr->b_angle(1,k,j,i),
                   pcr->b_angle(2,k,j,i),pcr->b_angle(3,k,j,i),v1,v2,v3);

            // only calculate v_dot_gradpc perpendicular to B
            // perpendicular direction only has flow velocity, no streaming velocity
          Real v_dot_gradpc = v2 * dpcdy + v3 * dpcdz;

          ec_source_(k,j,i) = v_dot_gradpc;

	  pcr->vgradpc(k,j,i) = v_dot_gradpc + v1*dpcdx; // Chad, Nov. 2020
	
         // if ((float(lcr/lb) > (1.e2/va)) && (pcr->use_piso > 0)) {
         //   pcr->vgradpc(k,j,i) = v_dot_gradpc; 
         // }

	  // New part to add in source from collisional losses -- Chad
	  //
	  if (pcr->gamma_flag > 0){
            Real dens2cgs = 6.85e-27;
            Real dcgs = w(IDN,k,j,i)*dens2cgs;  // density in cgs units
           /* Real pcgs = w(IEN,k,j,i)*pres2cgs;  // gas pressure in cgs units
            Real temp = 1.67e-24*pcgs/(dcgs*1.38e-16);
            Real nH = 0.6*dcgs/1.67e-24;
           
	    Real tUnit = 1.67e-24*pres2cgs/(dens2cgs*1.38e-16);
            Real lcgs = 3.0856e21; // 1 kpc is a length unit of 1? 
            Real tempCutoff = 1.6e4;
            Real coeff = 3.e3;
            fneutral = 0.5*(pcr->maxNeut)*(1.0 + tanh((-temp + tempCutoff)/(coeff)));

            Real fion = 1.0 - fneutral;
            fion = std::max(fion,1.e-4);
	    */
	    Real ne = (dcgs/1.67e-24)/1.18;
	    
            pcr->dedt_coul(k,j,i) = 2.78e-16 * ne * cr(CRE,k,j,i)*6.48e-11; // erg/s/cm^3
            pcr->dedt_pion(k,j,i) = 7.44e-16 * ne * cr(CRE,k,j,i)*6.48e-11; // erg/s/cm^3
 	  }
	
/*
	  if (pcr->shockInject > 0){
	        bool okShock = false;
      		Real dx = (pmb->pcoord->x1v(i)) - (pmb->pcoord->x1v(i-1));
     		// Real rho = w(IDN,k,j,i);
      		Real rho = w(IDN,k,j,i);
      		Real epsr = (std::abs(w(IDN,k,j,i+1) - 2.0*w(IDN,k,j,i) + w(IDN,k,j,i-1))
                    + std::abs(w(IDN,k,j+1,i) - 2.0*w(IDN,k,j,i) + w(IDN,k,j-1,i)))
                    /w(IDN,k,j,i);
      		Real epsp = (std::abs(w(IPR,k,j,i+1) - 2.0*w(IPR,k,j,i) + w(IPR,k,j,i-1))
                    + std::abs(w(IPR,k,j+1,i) - 2.0*w(IPR,k,j,i) + w(IPR,k,j-1,i)))
                    /w(IPR,k,j,i);
      		Real eps = std::max(epsr, epsp);
     		// maxeps = std::max(maxeps, eps);
		// refine : curvature > 0.01
	        if (eps > epsShock) okShock = true;
      		// derefinement: curvature < 0.005
      		if (eps < epsShock) okShock = false;
      		// otherwise, stay

      		Real pCRAdded = rampup*dt*(etaInject*rhoUp*velUp*velUp*velUp)/abs(dx);
      		// output variable is dP_cr/dt
   		//   std::cout << "added CR pressure: " << pCRAdded;
      		if (okShock == true) {
        		pmb->user_out_var(5,k,j,i) = rampup*(etaInject*rhoUp*velUp*velUp*velUp)/abs(dx);
       			// pcr->dedt_shock(k,j,i) = rampup*(etaInject*rhoUp*velUp*velUp*velUp)/abs(dx);
        		u_cr(CRE,k,j,i) = u_cr(CRE,k,j,i) + pCRAdded; //  /(4./3. - 1.);
        		cons(IEN,k,j,i) = cons(IEN,k,j,i) - pCRAdded/(1./3.); // /(1./3.); //  take same ENERGY away from gas
       			// u_cr(CRE,k,j,i) = u_cr(CRE,k,j,i) + pCRAdded/(rho*(4./3. - 1.));
      		}
      		else {
        		pmb->user_out_var(5,k,j,i) = 0.0;
        		u_cr(CRE,k,j,i) = u_cr(CRE,k,j,i);
      		}
      		std:cout << "eps: " << eps << "epsShock: " << epsShock << "okShock: " << okShock;


 	  }
*/
	  // Note: units are probably incorrect on loss terms
	  // ec_source_(k,j,i) += (dedt_pion + dedt_coul);

        //cons(IEN,k,j,i) += dt*(1.0/6.0)*dedt_pion + dt*dedt_coul;
	/*
        pmb->user_out_var(0,k,j,i) = dedt_coul;
        pmb->user_out_var(1,k,j,i) = dedt_pion;
        pmb->user_out_var(2,k,j,i) = dedt_coul + dedt_pion;
        pmb->user_out_var(3,k,j,i) = dt*(1.0/6.0)*dedt_pion + dt*dedt_coul;

        cout << "coulomb term: " << dedt_coul << endl;
        cout << "user variable 0: " << pmb->user_out_var(0,k,j,i) << endl;
	  w(IDN,k,j,i)

	*/
        }// end i

      }// end j

    }// end k
  }// end MHD

  //-----------------------------------------------------------------------
  // calculate coordinate source terms for Cosmic ray
  pco->AddCoordTermsDivergence(1,cr,coord_source_);

}


void CRIntegrator::FluxDivergence(const Real wght, AthenaArray<Real> &cr_out)
{
  CosmicRay *pcr=pmy_cr;
  MeshBlock *pmb = pcr->pmy_block;
  Coordinates *pco = pmb->pcoord;

  AthenaArray<Real> &x1flux=pcr->flux[X1DIR];
  AthenaArray<Real> &x2flux=pcr->flux[X2DIR];
  AthenaArray<Real> &x3flux=pcr->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

 
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) {

      // calculate x1-flux divergence 
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);

      for(int n=0; n<NCR; ++n){
 #pragma omp simd
        for(int i=is; i<=ie; ++i){
          dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }// end n
      }// End i

     // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        for(int n=0; n<NCR; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }// end nx2

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);

        for(int n=0; n<NCR; ++n){
#pragma omp simd
          for(int i=is; i<=ie; ++i){
            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }// end nx3
      // update variable with flux divergence
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for(int n=0; n<NCR; ++n){
#pragma omp simd
        for(int i=is; i<=ie; ++i){
          cr_out(n,k,j,i) -= wght*dflx(n,i)/vol(i);
        }
      }


    }// end j
  }// End k

  // Add coordinate source term
  for(int n=0; n<NCR; ++n)
    for(int k=ks; k<=ke; ++k)
      for(int j=js; j<=je; ++j)
#pragma omp simd
        for(int i=is; i<=ie; ++i)
          cr_out(n,k,j,i) += wght * coord_source_(n,k,j,i);  

  // check Ec is positive
  for(int k=ks; k<=ke; ++k)
    for(int j=js; j<=je; ++j)
#pragma omp simd
      for(int i=is; i<=ie; ++i){
        if(cr_out(CRE,k,j,i) < TINY_NUMBER)
          cr_out(CRE,k,j,i) = TINY_NUMBER;
      }

}
