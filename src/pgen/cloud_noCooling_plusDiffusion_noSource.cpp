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

using namespace std;  // Chad
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

static int n_user_var=9;


static AthenaArray<Real> cool_t;
static AthenaArray<Real> cool_coef;
static AthenaArray<Real> cool_index;
static AthenaArray<Real> cool_tef;
static int nline=4;

static Real cool_coef1=0.1414;
//static Real heat_coef=2;
static Real heat_coef=0.0481;
static Real cool_alpha1=-1.0;

static Real sigmaInput;
static Real tunit;
static Real rhounit;
static Real time_unit;
static Real gamma_idx;

// prototypes for user-defined diffusion, source terms, and boundary conditions
//void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
//        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt);

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void RadCool(MeshBlock *pmb, const Real time, const Real dt, 
     const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
     AthenaArray<Real> &cons);

void CRLosses(MeshBlock *pmb, const Real time, const Real dt,
     const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u_cr, AthenaArray<Real> &cons);

//void CRFluxSrcInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
//     int is, int ie, int js, int je, int ks, int ke);

//void CRReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
//     int is, int ie, int js, int je, int ks, int ke);

void CROutflowOuterX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

void FixMHDRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
     int ks, int ke, int ngh);

void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
     int ks, int ke, int ngh);

void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

void ExactCooling2(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);


//void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//{ 
//    
//}

// enroll user-defined boundary conditions and source terms
void Mesh::InitUserMeshData(ParameterInput *pin)
{
 // EnrollUserBoundaryFunction(inner_x1, FixMHDLeft);
 // EnrollUserBoundaryFunction(outer_x1, FixMHDRight);
  if(CR_ENABLED) {
  //  EnrollUserCRBoundaryFunction(inner_x1, FixCRsourceLeft);
 // if(CR_ENABLED){
  //  EnrollUserBoundaryFunction(INNER_X1, CRFluxSrcInnerX1);
  //  EnrollUserBoundaryFunction(OUTER_X1, CROutflowOuterX1);
 //   EnrollUserBoundaryFunction(inner_x1, CRFluxSrcInnerX1);
   // EnrollUserCRBoundaryFunction(outer_x1, CROutflowOuterX1);
   // EnrollUserExplicitSourceFunction(ExactCooling2);
  }
 // EnrollUserExplicitSourceFunction(RadCool);
 // EnrollUserCRGasSource(CRLosses); //? will this work?
 //
 //
 gamma_idx = pin->GetOrAddReal("hydro", "gamma", 5.0/3.0);
 sigmaInput = pin->GetOrAddReal("problem","sigmaInput",1.e30);
 
 if(CR_ENABLED){
    cool_t.NewAthenaArray(nline);
    cool_coef.NewAthenaArray(nline);
    cool_index.NewAthenaArray(nline);
    cool_tef.NewAthenaArray(nline);

    /*
    if ((fp = fopen("cool_func.dat","r")) == NULL) {
        std::cout << "### ERROR to read in the cooling function data" << std::endl
              << "Cannot open cool_func.dat" << std::endl;
      return;
    }
    */
    // temperature has nline numbers
    // coef and index only has nline-1 numbers

    /*

    for(int i=0; i<nline; ++i){
      fscanf(fp,"%lf",&(cool_t(i)));
      cool_t(i) /= tunit;
    }
    for(int i=0; i<nline; ++i){
      fscanf(fp,"%lf",&(cool_coef(i)));
    }
    for(int i=0; i<nline; ++i){
      fscanf(fp,"%lf",&(cool_index(i)));
    }
    */

    rhounit = 6.85e-27;
    Real pres2cgs = 6.54e-11;

    // use X=0.7, Z=0.02
    Real miu =0.617284;
    Real miue=1.17647;
    Real miuh=1.42857;
    Real kb=1.3807e-16;
    Real mp=1.6605e-24;
    
    tunit = miu*mp*pres2cgs/(rhounit*kb); // = 7.09E7 K
   // tunit = mp*pres2cgs/(rhounit*kb); // = 1.15E8 K
    time_unit = 3.155e13; // 1 Myr is how Josh set it up?? Note to self


    cool_coef(0) = 1.0012e-30; //* pow(2000.0,1.5);
    cool_coef(1) = 4.6240e-36; //* pow(8000.0,2.867);
    cool_coef(2) = 1.780e-18; //* pow(1.e5,-0.65);
    
    cool_index(0) = 1.5;
    cool_index(1) = 2.867;
    cool_index(2) = -0.65;
      
    cool_t(0) = 2000.0/tunit;
    cool_t(1) = 8000.0/tunit;
    cool_t(2) = 1.e5/tunit;
    cool_t(3) = 4.e7/tunit;

   // cool_t(0) = 1.e-10;
   // cool_t(1) = 0.02;
   // cool_t(2) = 1.0;
   // cool_t(3) = 1.e8;

   // cool_coef(0) = 0.3/pow(0.02,7.7);
   // cool_coef(1) = 0.3;
   // cool_coef(2) = 0.3;

   // cool_index(0) = 6.0;
   // cool_index(1) = -1.7;
   // cool_index(2) = 0.5;

    cool_tef(nline-1) = 0.0;



    // Scale the unit
    for(int i=0; i<nline-1; ++i){
      cool_coef(i) *= (pow(tunit,cool_index(i)) * (gamma_idx-1.0) * miu 
                       * rhounit * time_unit/(tunit * kb * mp * miue * miuh));
     // cool_coef(i) *= (pow(tunit,cool_index(i)) * (gamma_idx-1.0) 
      //                 * rhounit * time_unit/(tunit * kb * mp));
                      // * rhounit * time_unit/(tunit * kb * miu * mp));
                      // * rhounit * time_unit/(tunit * kb * mp * miue * miuh));
    }
    // After scale the unit, the equation for cooling is just:
    // dT/dt = -coef T^alpha in code unit
    // The solution is (T_ref/Lambda_ref)(Y(T^n) - Y(T^n+1))=-Delta t

    // The TEF is just Lambda(T_ref)/T_{ref} \int _T ^ref dT
    // Starting from Npoint, d
    cool_tef(nline-1) = 0.0;
    Real lambda_t_n=cool_coef(nline-2)*pow(cool_t(nline-1),cool_index(nline-2))
                                                            /cool_t(nline-1);

    for(int i=nline-2; i>=0; i--){
      Real slope = cool_index(i);
      Real coef = cool_coef(i);
      if(fabs(slope-1.0) < TINY_NUMBER){
        cool_tef(i) = cool_tef(i+1) + lambda_t_n*log(cool_t(i+1)/cool_t(i))/coef;
      }else{
        cool_tef(i) = cool_tef(i+1) + lambda_t_n*(pow(cool_t(i+1),1.0-slope) - 
                                       pow(cool_t(i),1.0-slope))/(coef*(1.0-slope));        
      }

    }
 
  //  fclose(fp);
    

  }          

 
 //
}


void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{

  cool_t.DeleteAthenaArray();
  cool_coef.DeleteAthenaArray();
  cool_index.DeleteAthenaArray();
  cool_tef.DeleteAthenaArray();

//  input_random_pot.DeleteAthenaArray();

}


// enroll user-defined CR diffusion function
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
    AllocateUserOutputVariables(n_user_var);
   // AllocateUserOutputVariables(2);
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
    
 // std::cout << "sigma in pgen: " << sigmaInput << "\n";
  Real gamma = peos->GetGamma();

  // get parameters from input file
  Real pgas  = pin->GetOrAddReal("problem","pgas", 1.);         // gas pressure
  Real rho_c = pin->GetOrAddReal("problem","rhocold",  1.);     // gas density inside cloud
  Real rho_h = pin->GetOrAddReal("problem","rhohot", .01);      // gas density outside cloud
  Real delta_r = pin->GetOrAddReal("problem","deltar", .01);    // thickness of cloud interface
  Real r_cloud = pin->GetOrAddReal("problem","rcloud", .1);     // radius of cloud
  Real x_cloud = pin->GetOrAddReal("problem","xcloud", 1.);     // x-coordinate of cloud center
  Real y_cloud = pin->GetOrAddReal("problem","ycloud", 0.);     // y-coordinate of cloud center

 // Real pgas = .00494;         // gas pressure
 // Real rho_c = 34.31;     // gas density inside cloud
 // Real rho_h = 0.3066;      // gas density outside cloud
 // Real delta_r = .00025;    // thickness of cloud interface
 // Real r_cloud = .05;     // radius of cloud
 // Real x_cloud = 1.0;     // x-coordinate of cloud center
 // Real y_cloud = 0.0;     // y-coordinate of cloud center

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
    Real ecr0 = pin->GetOrAddReal("problem","ecr0", 1.e-10);  // initial CR energy density (a small number)
   // Real ecr0 = 1.e-6;  // initial CR energy density (a small number)
    for(int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
           pcr->u_cr(CRE,k,j,i) = ecr0;
           pcr->u_cr(CRF1,k,j,i) = 0.0;                       // no initial CR flux
           pcr->u_cr(CRF2,k,j,i) = 0.0;
           pcr->u_cr(CRF3,k,j,i) = 0.0;
           pcr->dedt_coul(k,j,i) = 0.0;                            // Chad
           pcr->dedt_pion(k,j,i) = 0.0;                            // Chad
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
         // pcr->sigma_diff(0,k,j,i) = sigmaInput;
          pcr->sigma_diff(0,k,j,i) = 1./(pcr->diffInput);
         // pcr->sigma_diff(1,k,j,i) = sigma;
         // pcr->sigma_diff(2,k,j,i) = sigma;
          pcr->sigma_diff(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_diff(2,k,j,i) = pcr->max_opacity;
        }
      }
    }// end k,j,i

  }// End CR

  // set initial magnetic field (default: uniform horizontal)
  if(MAGNETIC_FIELDS_ENABLED){

    Real bx = pin->GetOrAddReal("problem","bx", 1.);
    Real by = pin->GetOrAddReal("problem","by", 0.);
    Real bz = pin->GetOrAddReal("problem","bz", 0.);
   // Real bx = .0349;
   // Real by = 0.;
   // Real bz = 0.;

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


void ExactCooling2(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

    // scale the temperature unit and cooling time
    // so the cooling will be just 
    // dT/dt = coef() Lambda(T)


  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is, iu=pmb->ie;


  Real lambda_t_n=cool_coef(nline-2)*pow(cool_t(nline-1),cool_index(nline-2))
                                                            /cool_t(nline-1);


  Real ncells1=pmb->block_size.nx1 + 2*(NGHOST);

//  cell_vol.NewAthenaArray(ncells1);


//  Real tot_heating = 0.0;
//  Real tot_cooling = 0.0;
//  Real tot_vol=0.0;

  Real Tfloor = 1.e4; 
  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
//      pmb->pcoord->CellVolume(k,j,il,iu,cell_vol);
      for(int i=il; i<=iu; ++i){
        // Need to take density and temperature at time step n from 
        // cons, not from prim
        // Because we do not need intermediate step to calculate the cooling 
        // function
        Real rho = cons(IDN,k,j,i);
        Real eint = cons(IEN,k,j,i)
                    - 0.5 *(cons(IM1,k,j,i)*cons(IM1,k,j,i)
                      + cons(IM2,k,j,i)*cons(IM2,k,j,i)
                      + cons(IM3,k,j,i)*cons(IM3,k,j,i))/rho;
        if(MAGNETIC_FIELDS_ENABLED){
             eint -= 0.5 * (bcc(IB1,k,j,i) * bcc(IB1,k,j,i)
                      + bcc(IB2,k,j,i) * bcc(IB2,k,j,i)
                      + bcc(IB3,k,j,i) * bcc(IB3,k,j,i));
        }
        Real t_i = eint *(gamma_idx - 1.0)/rho;

        
        if(t_i > cool_t(0)){

          if(t_i > cool_t(nline-1)) t_i = cool_t(nline-1);

          int t_loc=0;
          while((t_loc < nline-2) && (cool_t(t_loc+1) < t_i) ){
            ++t_loc;
          }

          Real slope = cool_index(t_loc);
          Real coef = cool_coef(t_loc);

          Real tef = cool_tef(t_loc+1);
          if(fabs(slope-1.0) < TINY_NUMBER){
            tef += lambda_t_n*log(cool_t(t_loc+1)/t_i)/coef;
          }else{
            tef += lambda_t_n*(pow(cool_t(t_loc+1),1.0-slope) -
                            pow(t_i,1.0-slope))/(coef*(1.0-slope));
          }
	 // cout << "slope: " << slope << " coef: " << coef << "\n";

          Real new_tef = tef + rho * dt * lambda_t_n;
          // Now invert TEF to get the current temperature
          // new_tef > tef
          int tef_loc=t_loc+1;
          while((tef_loc > 0) && (new_tef > cool_tef(tef_loc))){
            --tef_loc;
          }

          Real diff_tef = (new_tef - cool_tef(tef_loc+1))/lambda_t_n;
          slope = cool_index(tef_loc);
          coef = cool_coef(tef_loc);

          Real tnew = t_i;
          if(fabs(slope-1.0) < TINY_NUMBER){
            tnew = exp(log(cool_t(tef_loc+1))-(coef*diff_tef));
          }else{
            tnew = pow(cool_t(tef_loc+1),1.0-slope)
                               - (1.0-slope) * coef * diff_tef;
            tnew = pow(tnew,1.0/(1.0-slope));
          }
         /*
	 if ((t_i*tunit < 1.E6) && (t_i*tunit > 2.E4)) {
           std::cout << " t_i: " << t_i << ", t_i (cgs) " << t_i*tunit << "\n";
           std::cout << " bounding temps: " << cool_t(t_loc)*tunit << ", " << cool_t(t_loc+1)*tunit << "\n";
	  std::cout << "What it wants tnew to be: " << tnew*tunit << "\n";
	 }
	 */
	  tnew = std::max(Tfloor/tunit,tnew);
          // output variable is dE/dt?
          cons(IEN,k,j,i) += (tnew - t_i) * rho/(gamma_idx - 1.0);


	//  std::cout << "  Temperature old: " << t_i*tunit << "  Temperature new: " << tnew*tunit << "\n";

         // pmb->user_out_var(0,k,j,i) = (tnew - t_i) * rho/(gamma_idx - 1.0);

//          tot_cooling += pmb->user_out_var(0,k,j,i) * cell_vol(i);

          //add a constant heating rate -- Chad, April 9, 2020
          Real dens2cgs = 6.85e-27;
          Real dcgs = prim(IDN,k,j,i)*dens2cgs;
	  Real gHeat = 1.e-25;
          Real heat2code = 4.81e23;
          Real mp=1.6605e-24;
	  Real nH = 0.7*dcgs/mp;
	 /* if (tnew > Tfloor/tunit) {
            cons(IEN,k,j,i) += dt*nH*gHeat*heat2code;
          } */
         // std::cout << " t_i: " << t_i << ", t_i (cgs) " << t_i*tunit << "\n";
          
	 // std::cout << "dE/dt_cool: " << (tnew - t_i) * tunit * rhounit * rho/(gamma_idx - 1.0) / (dt * time_unit) << "\n";
	  
          cons(IEN,k,j,i) += dt*nH*gHeat*heat2code;
         // pmb->user_out_var(4,k,j,i) = - (tnew - t_i) * tunit * rhounit * rho/(gamma_idx - 1.0) / (dt * time_unit); // + nH*gHeat; 
          pmb->user_out_var(4,k,j,i) = ((tnew - t_i) * rho/(gamma_idx - 1.0))/(dt*heat2code) + nH*gHeat; // + nH*gHeat; 
         // if(tnew < 20.0){ // temperature cap on heating?  Note to self
        /*  if(tnew < 0.03){ // temperature cap on heating?  Note to self
            cons(IEN,k,j,i) += dt * rho* heat_coef/(gamma_idx-1.0);

           // pmb->user_out_var(1,k,j,i) = dt * rho * heat_coef/(gamma_idx-1.0);
//            tot_heating += cell_vol(i) * pmb->user_out_var(1,k,j,i);

          }else{
           // pmb->user_out_var(1,k,j,i) = 0.0;
          } */

//          tot_vol += cell_vol(i);

        }
      }
    }
  }

 /*   
    Real global_cooling=0.0;
    Real global_heating=0.0;
    Real global_vol=0.0;
    MPI_Allreduce(&tot_heating, &global_heating, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(&tot_cooling, &global_cooling, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(&tot_vol, &global_vol, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
    
    //this is the amount of energy we need to add back per unit volume, the same for all cells
    Real diff = -(global_heating + global_cooling)/global_vol;
    
    
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
        for(int i=il; i<=iu; ++i){
            cons(IEN,k,j,i) += diff;
        }
      }
    }

*/
//    cell_vol.DeleteAthenaArray();
  // Add the heating and cooling due to overal energy balance from last step
  /*
  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
      for(int i=il; i<=iu; ++i){
          cons(IEN,k,j,i) += -dt * pmb->user_out_var(6,k,j,i);
      }
    }
  } */  //Chad -- Got rid of this April 9, 2020

}



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

       // pcr->sigma_diff(0,k,j,i) = sigmaInput;
        pcr->sigma_diff(0,k,j,i) = 1./(pcr->diffInput);
        pcr->sigma_diff(1,k,j,i) = pcr->max_opacity;
        pcr->sigma_diff(2,k,j,i) = pcr->max_opacity;  

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
	
	  /////////////////////////////////////////////////////
	  Real va1 = bcc(IB1,k,j,i) * inv_sqrt_rho;
          Real va2 = bcc(IB2,k,j,i) * inv_sqrt_rho;
          Real va3 = bcc(IB3,k,j,i) * inv_sqrt_rho;

          Real va = sqrt(pb) * inv_sqrt_rho;
          Real dpc_sign = 0.0;

          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = 1.0;
          else if(-(pcr->b_grad_pc(k,j,i)) > TINY_NUMBER) dpc_sign = -1.0;

          pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
          pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
          pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

          Real stream_boost_factor = 1.0;
          if (pcr->ion_alfven_speed > 0) { // Chad, March 8, 2020
            Real fneutral = 1.2;
            Real dens2cgs = 6.85e-27;
            Real pres2cgs = 6.54e-11;
            Real dcgs = prim(IDN,k,j,i)*dens2cgs;  // density in cgs units
            Real pcgs = prim(IEN,k,j,i)*pres2cgs;  // gas pressure in cgs units
           // Real temp = 0.6*1.67e-24*pcgs/(dcgs*1.38e-16);
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
            Real lcgs = 3.0856e21; // 1 kpc is a length unit of 1? 
            Real tempCutoff = 1.6e4;
           // Real coeff = 3.e3;
            Real coeff = 2.e3;
            fneutral = 0.5*(pcr->maxNeut)*(1.0 + tanh((-temp + tempCutoff)/(coeff)));

           // std::cout << "fneutral: " << fneutral << " maxNeut: " << pcr->maxNeut << " tanh term: " << tanh((-temp + tempCutoff)/(coeff)) << "\n";
            Real fion = 1.0 - fneutral;
            fion = std::max(fion,1.e-8);
          //  fneutral = 1.0 - fneutral/1.20;

            // Real temp = (1.377e14*1.67e-24/1.38e-16)*w(IDN,k,j,i)/w(IEN,k,j,i);
            // assumes the fiducial temp outside the cloud is 1e6 when P = 1/gamma, rho = 1, and 
            // I assume that rho = 1 corresponds to an actual rho of 1e-26
            //Real dampRate = 1.e-9 * fneutral * (temp / 
            //              1000.0)**0.5 * w(IDN,k,j,i)/1.e2; 
            Real dampRate = 1.e-9 * fneutral * sqrt(temp / 1000.0) * dcgs/1.e-24;
           // stream_boost_factor = 1.0 + (4.0*3.0e10*dampRate*0.5*pb*pres2cgs
           //                   /(4.803e-10*pow(pb*pres2cgs,1.5)*2.0*(sqrt(8*3.1415)*3.1415/dcgs)*w(IDN,k,j,i)*cr(CRE,k,j,i)*pres2cgs));
            stream_boost_factor = 1.0/sqrt(fion);

           // std::cout << "b_grad_pc: " << b_grad_pc << "\n";
           // std::cout << "b_grad_pc: " << sqrt(pow(b_grad_pc,2.)) << "\n";
            Real lcr = 1./3.*prim(IDN,k,j,i)*u_cr(CRE,k,j,i) / std::max(TINY_NUMBER,sqrt(pow(pcr->b_grad_pc(k,j,i),2.)));
           // std::cout << "lcr: " << lcr << "\n";
            lcr *= sqrt(pb);  // to get Pcr/|grad_parallel(Pcr)|


            lcr *= lcgs;
            Real damping_diff_term = (4.0*3.0e10*dampRate*0.5*pb*pres2cgs*fion*4./3.*lcr
                                /(4.803e-10*pow(pb*pres2cgs,1.5)*2.0*(sqrt(8*3.1415)*3.1415/dcgs)*prim(IDN,k,j,i)*u_cr(CRE,k,j,i)*pres2cgs));


            damping_diff_term *= va/sqrt(fion);
            damping_diff_term *= 3.155e13/(lcgs*lcgs);  // back into code units now?
	    damping_diff_term *= 1./624.15; // because gamma_L is in GeV, and 1/624 = 1 GeV in ergs

            va1 *= stream_boost_factor;  // Chad, March 28, 2020
            va2 *= stream_boost_factor;  // Chad, March 28, 2020
            va3 *= stream_boost_factor;  // Chad, March 28, 2020

            va1 = std::min(va1,pcr->max_ion_va);
            va2 = std::min(va2,pcr->max_ion_va);
            va3 = std::min(va3,pcr->max_ion_va);
            va = sqrt(pow(va1,2.) + pow(va2,2.) + pow(va3,2.));

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
                pcr->sigma_diff(0,k,j,i) = 1./(pcr->diffInput + std::min(damping_diff_term, pcr->max_in_diff));  // Chad, March 28, 2020
              
            	}
	    }
	    else {
	//	pcr->sigma_diff(0,k,j,i) = sigmaInput; // user input for sigma
		pcr->sigma_diff(0,k,j,i) = 1./(pcr->diffInput); // user input for sigma
	    }
          
	  }

          if(va > TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = fabs(pcr->b_grad_pc(k,j,i))/(sqrt(pb) * va *
                                   (4.0/3.0) * invlim * u_cr(CRE,k,j,i));
           /* if (pcr->ion_neutral_friction > 0) {
              pcr->sigma_adv(0,k,j,i) = fabs(b_grad_pc)/(sqrt(pb) * 
                                   std::min(va * stream_boost_factor,1.0) * 
                                   (4.0/3.0) * invlim * cr(CRE,k,j,i));
            } */
            pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
            pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
          }





	  //////////////////////////////////////////////////////////////

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
           }
	   else {
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

/*
	  // The added part where I calculate CR loss terms and put them in uov -- Chad
	  Real dens2cgs = 6.85e-27;
	  Real dcgs = prim(IDN,k,j,i)*dens2cgs;  // density in cgs units
          Real ne = dcgs/1.67e-24 / 1.18e0;   // me
          Real coul = 2.78e-16 * ne * u_cr(CRE,k,j,i);
          Real pion = 7.44e-16 * ne * u_cr(CRE,k,j,i);

          pmb->user_out_var(0,k,j,i) = coul;
          pmb->user_out_var(1,k,j,i) = pion;
          pmb->user_out_var(2,k,j,i) = -coul - pion; // CR loss part
          pmb->user_out_var(3,k,j,i) = (1.0/6.0)*pion + coul;
	  cout << coul << endl;
	  cout << pmb->user_out_var(0,k,j,i) << endl;
*/
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
   // Real ramptime = 3.0;
   // Real time_factor = (1.-exp(-time/ramptime))*(1.-exp((time-30.0)/ramptime)); // time-dependence of CR source
   // time_factor = std::max(time_factor, 0.);
   // Real xflux = 2.4e-5*time_factor;
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


/* void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
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
void FixMHDRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
     int ks, int ke, int ngh)
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
}

void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
     int ks, int ke, int ngh)
{

  // copy hydro variables into ghost zones, reflecting v1

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {

        prim(IDN,k,j,is-i) = prim(IDN,k,j,is);
       // prim(IVX,k,j,is-i) = -prim(IVX,k,j,is); // reflect 1-velocity
        prim(IVX,k,j,is-i) = prim(IVX,k,j,is); // reflect 1-velocity
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

void CROutflowOuterX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{

/*
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
*/
  if(CR_ENABLED){
    for (int n=0; n<(NCR); ++n) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          // enforce slight outward gradient to allow CRs to leave
          if(n==CRE) {
            u_cr(n,k,j,ie+i) = 0.99999999*u_cr(n,k,j,ie);
          } else {
            u_cr(n,k,j,ie+i) = u_cr(n,k,j,ie);
          }
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

	pmb->user_out_var(4,k,j,i) = -(nH*nH*lambda - nH*gHeat); // in cgs (dE/dt)
      }
    }
  } // end i,j,k
  return;
}

void CRLosses(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u_cr,AthenaArray<Real> &cons)
{
  Real dens2cgs = 6.85e-27;
  Real pres2cgs = 6.54e-11;
  Real heat2code = 4.81e23;

  cout << "Fix CRLosses is being called " << endl;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real dcgs = prim(IDN,k,j,i)*dens2cgs;  // density in cgs units
        Real ne = dcgs/1.67e-24 / 1.18e0;   // me
        Real dedt_coul = 2.78e-16 * ne * u_cr(CRE,k,j,i);
        Real dedt_pion = 7.44e-16 * ne * u_cr(CRE,k,j,i);
        u_cr(CRE,k,j,i) -= dt*(dedt_pion + dedt_coul);

	cons(IEN,k,j,i) += dt*(1.0/6.0)*dedt_pion + dt*dedt_coul;

	pmb->user_out_var(0,k,j,i) = dedt_coul;
	pmb->user_out_var(1,k,j,i) = dedt_pion;
	pmb->user_out_var(2,k,j,i) = dedt_coul + dedt_pion;
	pmb->user_out_var(3,k,j,i) = dt*(1.0/6.0)*dedt_pion + dt*dedt_coul;

	cout << "coulomb term: " << dedt_coul << endl;
	cout << "user variable 0: " << pmb->user_out_var(0,k,j,i) << endl;
      }
    }
  }
  return;
}
