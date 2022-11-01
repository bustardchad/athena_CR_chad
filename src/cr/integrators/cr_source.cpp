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
//! \file rad_source.cpp
//  \brief Add radiation source terms to both radiation and gas
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../hydro/hydro.hpp"
#include "../../field/field.hpp"
#include "../../eos/eos.hpp"
#include "../../utils/utils.hpp"
#include <iostream>   // endl
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


using namespace std;
//add the source terms implicitly



void CRIntegrator::AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &bcc,
        AthenaArray<Real> &u_cr)
{
  CosmicRay *pcr=pmb->pcr;

  Real vlim = pcr->vmax;
  Real invlim = 1.0/vlim;

// The information stored in the array
// b_angle is
// b_angle[0]=sin_theta_b
// b_angle[1]=cos_theta_b
// b_angle[2]=sin_phi_b
// b_angle[3]=cos_phi_b
  
    
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
 
  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){

         Real fxx = 1.0/3.0;
         Real fyy = 1.0/3.0;
         Real fzz = 1.0/3.0;
         Real fxy = 0.0;
         Real fxz = 0.0;
         Real fyz = 0.0;

         Real *ec = &(u_cr(CRE,k,j,0));
         Real *fc1 = &(u_cr(CRF1,k,j,0));
         Real *fc2 = &(u_cr(CRF2,k,j,0));
         Real *fc3 = &(u_cr(CRF3,k,j,0));

         // The angle of B
         Real *sint_b = &(pcr->b_angle(0,k,j,0));
         Real *cost_b = &(pcr->b_angle(1,k,j,0));
         Real *sinp_b = &(pcr->b_angle(2,k,j,0));
         Real *cosp_b = &(pcr->b_angle(3,k,j,0));

 // adv1 is dPc/dx, adv2 is dPc/dy, adv3 is dPc/dz

      for(int i=is; i<=ie; ++i){
         Real rho = u(IDN,k,j,i);
         Real v1 = u(IM1,k,j,i)/rho;
         Real v2 = u(IM2,k,j,i)/rho;
         Real v3 = u(IM3,k,j,i)/rho;
         Real vtot1 = v1;
         Real vtot2 = v2;
         Real vtot3 = v3;

         // add the streaming velocity
         if(pcr->stream_flag > 0){
           vtot1 += pcr->v_adv(0,k,j,i);
           vtot2 += pcr->v_adv(1,k,j,i);
           vtot3 += pcr->v_adv(2,k,j,i);
         }

         // Chad -- new v_A vectors, these are defined in x, y, z coordinates, not with respect to B field
         Real vst1 = pcr->v_adv(0,k,j,i);
         Real vst2 = pcr->v_adv(1,k,j,i);
         Real vst3 = pcr->v_adv(2,k,j,i);

         Real fr1 = fc1[i];
         Real fr2 = fc2[i];
         Real fr3 = fc3[i];

         // in the case with magnetic field
        // rotate the vectors to oriante to the B direction
         if(MAGNETIC_FIELDS_ENABLED){
           // Apply rotation of the vectors
           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],v1,v2,v3);

           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],fr1,fr2,fr3);

           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],vtot1,vtot2,vtot3);
           
           RotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],vst1,vst2,vst3); // Chad -- rotate new va vectors too

           // perpendicular energy source term is already added via ec_source_
           // we will still need to include this in general for diffusion case
           vtot2 = 0.0;
           vtot3 = 0.0;

         }

         Real sigma_x = pcr->sigma_diff(0,k,j,i);
         Real sigma_y = pcr->sigma_diff(1,k,j,i);
         Real sigma_z = pcr->sigma_diff(2,k,j,i);

         if(pcr->stream_flag > 0){
           sigma_x = 1.0/(1.0/pcr->sigma_diff(0,k,j,i) + 
                             1.0/pcr->sigma_adv(0,k,j,i));

           sigma_y = 1.0/(1.0/pcr->sigma_diff(1,k,j,i) +
                             1.0/pcr->sigma_adv(1,k,j,i));

           sigma_z = 1.0/(1.0/pcr->sigma_diff(2,k,j,i) + 
                             1.0/pcr->sigma_adv(2,k,j,i));
         }

         // Now update the momentum equation
         //\partial F/\partial t=-V_m\sigma (F-v(E+Pc_)/v_m)) 
         // And the energy equation
         //\partial E_c/\partial t = -vtot sigma (F- v(E_c+P_c)/v_m)

         Real rhs1 = ec[i];
         Real rhs2 = fr1;
         Real rhs3 = fr2;
         Real rhs4 = fr3;

         Real coef_11 = 1.0 - dt * sigma_x * vtot1 * v1 * invlim * 4.0/3.0
                            - dt * sigma_y * vtot2 * v2 * invlim * 4.0/3.0
                            - dt * sigma_z * vtot3 * v3 * invlim * 4.0/3.0;
         Real coef_12 = dt * sigma_x * vtot1;
         Real coef_13 = dt * sigma_y * vtot2;
         Real coef_14 = dt * sigma_z * vtot3;

         Real coef_21 = -dt * v1 * sigma_x * 4.0/3.0;
         Real coef_22 = 1.0 + dt * vlim * sigma_x;

         Real coef_31 = -dt * v2 * sigma_y * 4.0/3.0;
         Real coef_33 = 1.0 + dt * vlim * sigma_y;

         Real coef_41 = -dt * v3 * sigma_z * 4.0/3.0;
         Real coef_44 = 1.0 + dt * vlim * sigma_z;

        //newfr1 = (rhs2 - coef21 * newEc)/coef22
        // newfr2= (rhs3 - coef31 * newEc)/coef33
        // newfr3 = (rhs4 - coef41 * newEc)/coef44
        // coef11 - coef21 * coef12 /coef22 - coef13 * coef31 /coef33 - coef41 * coef14 /coef44)* newec
        //    =rhs1 - coef12 *rhs2/coef22 - coef13 * rhs3/coef33 - coef14 * rhs4/coef44

        Real e_coef = coef_11 - coef_12 * coef_21/coef_22 - coef_13 * coef_31/coef_33 
                       - coef_14 * coef_41/coef_44;
        Real new_ec = rhs1 - coef_12 * rhs2/coef_22 - coef_13 * rhs3/coef_33 
                      - coef_14 * rhs4/coef_44;

        new_ec /= e_coef;

        Real newfr1 = (rhs2 - coef_21 * new_ec)/coef_22;
        Real newfr2 = (rhs3 - coef_31 * new_ec)/coef_33;
        Real newfr3 = (rhs4 - coef_41 * new_ec)/coef_44;
   
        pcr->dedt_damping(k,j,i) = vst1*sigma_x*(newfr1 - new_ec * v1 * invlim * 4.0/3.0)
                                  + vst2*sigma_y*(newfr2 - new_ec * v2 * invlim * 4.0/3.0)
                                  + vst3*sigma_z*(newfr3 - new_ec * v3 * invlim * 4.0/3.0);


        Real heat2code = 4.81e23; // Chad
        // Now apply the invert rotation
        if(MAGNETIC_FIELDS_ENABLED){
         // Apply rotation of the vectors
          InvRotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i], 
                                         newfr1,newfr2,newfr3);

          // updated dedt_damping (heating term) using flux, sigma, etc. instead of naive va \cdot grad Pcr
          Real s1 = sigma_x*(rhs2 - rhs1 * 4./3. * invlim * v1); // Chad -- need to check this
         // pcr->dedt_damping(k,j,i) = abs(vst1*s1);
          if (pcr->heating_off > 0){ 
            new_ec += dt * (ec_source_(k,j,i) + pcr->dedt_damping(k,j,i));
	  } 
          else {
            new_ec += dt * ec_source_(k,j,i);
          }
        }        

         // Add the energy source term
         // Is this the heating term? If so, super-Alfvenic streaming will have to be handled a little carefully
         if (NON_BAROTROPIC_EOS && (pcr->src_flag > 0)){
           Real new_eg = u(IEN,k,j,i) - (new_ec - ec[i]); 
           if(new_eg < 0.0) new_eg = u(IEN,k,j,i);
          // u(IEN,k,j,i) = new_eg;   
          // Chad    
           u(IEN,k,j,i) = new_eg + heat2code*dt*((1.0/6.0)*pcr->dedt_pion(k,j,i) + pcr->dedt_coul(k,j,i));       
         }         

         if(new_ec < 0.0) new_ec = ec[i];
      
         // JW: comment next 3 lines out to decouple gas from CR pressure gradients
         if((pcr->src_flag > 0) && (pmb->pmy_mesh->time > pcr->start_cr_force_time)){ 
        // if(pcr->src_flag > 0){ 
           u(IM1,k,j,i) += (-(newfr1 - fc1[i]) * invlim);
           u(IM2,k,j,i) += (-(newfr2 - fc2[i]) * invlim);
           u(IM3,k,j,i) += (-(newfr3 - fc3[i]) * invlim);
         }


	// cout << pcr->dedt_coul(k,j,i) << endl;
        // u_cr(CRE,k,j,i) = new_ec;
        // Chad
         u_cr(CRE,k,j,i) = new_ec - heat2code*dt*(pcr->dedt_coul(k,j,i) + pcr->dedt_pion(k,j,i));
         u_cr(CRF1,k,j,i) = newfr1;
         u_cr(CRF2,k,j,i) = newfr2;
         u_cr(CRF3,k,j,i) = newfr3;

	 // uov are in cgs units
	// pmb->user_out_var(0,k,j,i) = pcr->dedt_coul(k,j,i);
	 pmb->user_out_var(0,k,j,i) = pcr->uncapped_diffusion(k,j,i); //lcr in cgs units
         pmb->user_out_var(1,k,j,i) = pcr->dedt_pion(k,j,i);
         pmb->user_out_var(2,k,j,i) = -pcr->dedt_coul(k,j,i) - pcr->dedt_pion(k,j,i); // CR loss part
         pmb->user_out_var(3,k,j,i) = (1.0/6.0)*pcr->dedt_pion(k,j,i) + pcr->dedt_coul(k,j,i); // thermalizes the gas
         pmb->user_out_var(5,k,j,i) = sqrt(pow(pcr->v_adv(0,k,j,i),2.0) + pow(pcr->v_adv(1,k,j,i),2.0) 
			+ pow(pcr->v_adv(2,k,j,i),2.0));
         pmb->user_out_var(6,k,j,i) = sqrt(pow(pcr->sigma_diff(0,k,j,i),-2.0) + pow(pcr->sigma_diff(1,k,j,i),-2.0) 
			+ pow(pcr->sigma_diff(2,k,j,i),-2.0));
	// pmb->user_out_var(7,k,j,i) = ec_source_(k,j,i); //NOT the energy loss rate from wave damping/heating of gas
	 pmb->user_out_var(7,k,j,i) = pcr->dedt_damping(k,j,i); //CORRECT energy loss rate from wave damping/heating of gas
	 pmb->user_out_var(8,k,j,i) = pcr->crscale(k,j,i); //lcr in cgs units
	 pmb->user_out_var(9,k,j,i) = pcr->vgradpc(k,j,i); //flow velocity dot grad Pc
	 pmb->user_out_var(10,k,j,i) = pcr->bscale(k,j,i); //magnetic length scale in cgs units
        // std::cout << pmb->user_out_var(0,k,j,i) << endl;
	// if (pmb->user_out_var(5,k,j,i) > 5) {
	//   std::cout << pcr->v_adv(0,k,j,i);
	// }

        // pmb->user_out_var(5,k,j,i) = pcr->dedt_shock(k,j,i);
        // cout << pcr->dedt_coul(k,j,i) << endl;
        // cout << pmb->user_out_var(0,k,j,i) << endl;
         
      }// end 

    }// end j
  }// end k


  // Add user defined source term for cosmic rays
  if(pcr->cr_source_defined)
   // std::cout << "calling user defined source term: ";
   // pcr->UserSourceTerm_(pmb, pmb->pmy_mesh->time, dt, w,bcc, u_cr);
   // std::cout << "called first one: ";
    pcr->UserSourceTerm2_(pmb, pmb->pmy_mesh->time, dt, w,bcc, u_cr, u);
   // std::cout << "called second one: ";
      
}


