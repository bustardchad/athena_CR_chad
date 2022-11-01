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
#include <random>

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

AthenaArray<Real> fboundary, Phii1Vec, Phii2Vec, Phii3Vec, AiiVec; // here too, so they
using namespace std;

// postshock flow variables are shared with IIB function
namespace {
Real presUp, rhoUp, velUp, B0, BAngle, shockPos;
Real L, alpha_exp;
int sim_kx1, sim_ky1, sim_kx2, sim_ky2;
bool turbClumps;
// can be used in boundary functions
} // namespace
//======================================================================================
/*! \file coolpulse.cpp
 *  \brief Inject CRs in a pulse at the left boundary, toward a cool cloud
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief CR injection with radiative cooling
//======================================================================================

static Real pi = 3.14159265358979323846;
static Real sigma=1.e30; // some large number to prevent diffusion (large enough?)
static int n_user_var=5;
// prototypes for user-defined diffusion, source terms, and boundary conditions
//void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
//        AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Real dt);

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void RadCoolCGS(MeshBlock *pmb, const Real time, const Real dt, 
     const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
     AthenaArray<Real> &cons);

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

//void FixMHDRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
//     int ks, int ke, int ngh);

//void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//     FaceField &b, Real time, Real dt, int is, int ie, int js, int je,
//     int ks, int ke, int ngh);

//void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
//    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
//    int js, int je, int ks, int ke, int ngh);

void ShockOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                       FaceField &b, Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 
    
}

// enroll user-defined boundary conditions and source terms
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserBoundaryFunction(outer_x1, ShockOuterX1);
  if(CR_ENABLED) {
    EnrollUserCRBoundaryFunction(outer_x1, CROutflowOuterX1);
  }
 // EnrollUserExplicitSourceFunction(RadCoolCGS);
 //
 
 // AthenaArray<Real> fboundary, Phii1Vec, Phii2Vec, Phii3Vec, AiiVec; // here too, so they
  
  sim_kx1 = pin->GetOrAddReal("problem","sim_kx1", 10);
  sim_ky1 = pin->GetOrAddReal("problem","sim_ky1", 10);
  sim_kx2 = pin->GetOrAddReal("problem","sim_kx2", 50);
  sim_ky2 = pin->GetOrAddReal("problem","sim_ky2", 50);
  alpha_exp = pin->GetOrAddReal("problem","alpha_exp", 0.01);
  AllocateRealUserMeshDataField(4);
 // ruser_mesh_data[0].NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);
 // ruser_mesh_data[1].NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);
 // ruser_mesh_data[2].NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);
 // ruser_mesh_data[3].NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);

  //
  Phii1Vec.NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);
  Phii2Vec.NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);
  Phii3Vec.NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);
  AiiVec.NewAthenaArray(sim_kx2,sim_ky2,sim_ky2);

  // See Bruggen 2013 for details
   for (int kx = sim_kx1; kx<=sim_kx2; ++kx) {
     for (int ky = sim_ky1; ky<=sim_ky2; ++ky) {
        for (int kz = sim_ky1; kz<=sim_ky1; ++kz) { // use the same for y and z ! Chad -- 2/11/18
          Phii1Vec(kx,ky,kz) = 2.*pi*rand() / double(RAND_MAX);
          Phii2Vec(kx,ky,kz) = 2.*pi*rand() / double(RAND_MAX);
          Phii3Vec(kx,ky,kz) = 2.*pi*rand() / double(RAND_MAX);
	  AiiVec(kx,ky,kz) = rand() / double(RAND_MAX);
  	  shockPos = pin->GetOrAddReal("problem","shockPos", 0.0);
  	  L = pin->GetOrAddReal("problem","L", 8.E23);
	
        //  ruser_mesh_data[0](kx,ky,kz) = 2.*pi*rand() / double(RAND_MAX);
        //  ruser_mesh_data[1](kx,ky,kz) = 2.*pi*rand() / double(RAND_MAX);
        //  ruser_mesh_data[2](kx,ky,kz) = 2.*pi*rand() / double(RAND_MAX);
	//  ruser_mesh_data[3](kx,ky,kz) = rand() / double(RAND_MAX);
        }
     }
  }
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
    
  Real gamma = peos->GetGamma();

  // get parameters from input file
  Real pgas  = pin->GetOrAddReal("problem","pgas", 1.);         // gas pressure
  Real rho_c = pin->GetOrAddReal("problem","rhocold",  1.);     // gas density inside cloud
 // Real rho_h = pin->GetOrAddReal("problem","rhohot", .01);      // gas density outside cloud
  Real delta_r = pin->GetOrAddReal("problem","deltar", .01);    // thickness of cloud interface
  Real r_cloud = pin->GetOrAddReal("problem","rcloud", .1);     // radius of cloud
  Real x_cloud = pin->GetOrAddReal("problem","xcloud", 1.);     // x-coordinate of cloud center
  Real y_cloud = pin->GetOrAddReal("problem","ycloud", 0.);     // y-coordinate of cloud center

  Real r_cloud2 = pin->GetOrAddReal("problem","rcloud2", .1);     // radius of cloud
  Real x_cloud2 = pin->GetOrAddReal("problem","xcloud2", 1.);     // x-coordinate of cloud center
  Real y_cloud2 = pin->GetOrAddReal("problem","ycloud2", 0.);     // y-coordinate of cloud center
 // Real pgas = .00494;         // gas pressure
 // Real rho_c = 34.31;     // gas density inside cloud
 // Real rho_h = 0.3066;      // gas density outside cloud
 // Real delta_r = .00025;    // thickness of cloud interface
 // Real r_cloud = .05;     // radius of cloud
 // Real x_cloud = 1.0;     // x-coordinate of cloud center
 // Real y_cloud = 0.0;     // y-coordinate of cloud center

  // New parameters for shock problem: (added by Chad, January 28, 2020)
  velUp = pin->GetOrAddReal("problem","velUp", 1.E6);
  rhoUp = pin->GetOrAddReal("problem","rhoUp", 1.E-26);
  Real TUp = pin->GetOrAddReal("problem","TUp", 1.E7);
  BAngle = pin->GetOrAddReal("problem","BAngle", 0.0);
  B0 = pin->GetOrAddReal("problem","B0", 0.0);
  bool magJump = pin->GetOrAddBoolean("problem","magJump", true);

  // for clumps
  bool multipleClouds = pin->GetOrAddBoolean("problem","multipleClouds", false);
  turbClumps = pin->GetOrAddBoolean("problem","turbClumps", false);

 // int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

 // f.NewAthenaArray(nc3,nc2,nc1);

  Real BDown = 0.0;
  Real velDown = 0.0;
  Real rhoDown = 0.0;
  Real presDown = 0.0;

  Real presZone = 0.0;
  Real velxZone = 0.0;
  Real rhoZone = 0.0;
  Real BZone = 0.0;
        
  // March 28, 2017 -- changed 2.E-24 to 1.66E-24
  presUp = 1.380658E-16 * TUp * rhoUp / (1.6726E-24);  // upstream pressure
  Real Mach = sqrt(rhoUp*(pow(velUp,2.0))/(gamma*presUp));  // Mach number
  // Factors that go into Rankine-Hugoniot jump conditions:
  Real factor1 = ((gamma + 1.)*pow(Mach,2.0))/((gamma-1.)*pow(Mach,2.0) + 2.);
  Real factor2 = (2.0*(gamma)*pow(Mach,2.0) - (gamma - 1.))/((gamma) + 1.);

  // !!!!!!!!!!!!!!!!!
  if (magJump == true) {
  /*
        ! using magnetized jump condition for B field parallel to shock
        ! taken from Wolfram site for solving a cubic
    !    a2 = (velUp/((Mach**2)*(sim_gamma_gas - 1.)) + &
    !            velUp*sim_gamma_gas/(sim_gamma_gas - 1.) + &
    !            (B0*cos(BAngle*pi/180.))**2*sim_gamma_gas / &
    !            (8*pi*rhoUp*velUp*(sim_gamma_gas -1.)) ) / &
    !            (1./2. - (sim_gamma_gas/(sim_gamma_gas - 1.)))
    !    print *, a2
    !    a1 = (-(B0*cos(BAngle*pi/180.))**2/(4*pi*rhoUp) + &
    !            (velUp**2)*(1./2.+1./((Mach**2)*(sim_gamma_gas - 1.)))) / &
    !            (1./2. - (sim_gamma_gas/(sim_gamma_gas - 1.)))
    !    a0 = (-((B0*cos(BAngle*pi/180.))**2/(4*pi*rhoUp))*velUp*sim_gamma_gas / &
    !            (2*(sim_gamma_gas - 1.) -1.)) / &
    !            (1./2. - (sim_gamma_gas/(sim_gamma_gas - 1.)))
    !    Q = (3*a1-a2**2)/9
    !    R = (9*a2*a1-(27*a0)-(2*a2**3))/54    
    !    D = (Q**3)+(R**2)
    !    S = (R+sqrt(D))**(1./3.)          
    !    T = (R-sqrt(D))**(1./3.) 

    !    !solutions to the cubic: Cardano's formula
    !    z1 = -a2/3. + (S+T)
    !    z2 = -a2/3. - (S+T)/2. + (sqrt(-1.)*sqrt(3.)/2)*(S-T)
    !    z3 = -a2/3. - (S+T)/2. - (sqrt(-1.)*sqrt(3.)/2)*(S-T)        
    !    print *, "velDown solved from cubic: "
    !    print *, z1
    !    print *, z2
    !    print *, z3
*/

	Real z1 = 0.0;
	Real z2 = 0.0;
	Real root = 1.;
       // !instead, use quadratic formula from UT farside website
        if (BAngle == 90.0) {
             // cout << "BAngle = 90";
        }
        else {
                Real beta = 8.0*pi*presUp/pow((B0*cos(BAngle*pi/180.)),2.0);
                Real a = 2.0*(2. - gamma);
                Real b = gamma*(2.0*(1.+beta) + (gamma - 1.)*beta*pow(Mach,2.0));
                Real c = -gamma*(gamma + 1.)*beta*pow(Mach,2.0);

                z1 = (-b + sqrt(pow(b,2.0) - 4.0*a*c))/(2.0*a);
                z2 = (-b - sqrt(pow(b,2.0) - 4.0*a*c))/(2.0*a);
		
             // !  print *, "roots solved from quadratic: "
            //  !  print *, z1
            //  !  print *, z2

                if ((z1 > 0.0) && (z1 <= 4.0)) {
                        root = z1;
		}
                else if ((z2 > 0.0) && (z2 <= 4.0)) {
                        root = z2;
		}
                else {
                      //  print *, "No valid roots"
                }
                

               // print *, "Compression ratio: "
               // print *, root
                factor1 = root;
                factor2 = 1 + (gamma)*(pow(Mach,2.0))*(1. - (1./root)) + 
                        ((1./beta)*(1-(pow(root,2.0))));
        }
        //! jump conditions:
        BDown = B0*factor1;
        rhoDown = factor1 * rhoUp;
        velDown = velUp/factor1;
        presDown = factor2 * presUp;
  }
  else if (magJump == false) {
        rhoDown = factor1 * rhoUp;
        velDown = velUp/factor1;
        presDown = factor2 * presUp;
        BDown = B0;
  }
  else  {
        rhoDown = rhoUp;
        velDown = velUp;
        presDown = presUp;
        BDown = B0;
  }

  // log-normal distribution of clumps -- Chad, January 30, 2019
  // phii1, phii2, phii3, Ai are randomly generated here for clumps
  // AiiVec, Phii1Vec, Phii2Vec, Phii3Vec hold those random values for each
  // kx,ky,kz

  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        Real xcoord = pcoord->x1v(i);
  	Real ycoord = pcoord->x2v(j);

  //	f(i,j,k) = 0.0;
	Real fnew = 0.0;
  	for (int kx = sim_kx1; kx<=sim_kx2; ++kx) {
    	  for (int ky = sim_ky1; ky<=sim_ky2; ++ky) {
      	    for (int kz = sim_ky1; kz<=sim_ky1 ;++kz) { // use the same for y and z
        	//f(i,j,k) = f(i,j,k) + AiiVec(kx,ky,kz)*sin((kx*xcoord*2.0*pi/L+Phii1Vec(kx,ky,kz)))*sin(ky*ycoord*2.0*pi/L + Phii2Vec(kx,ky,kz)); //*Sin(kz*z(k)*2.0*PI/L+Phii3Vec(kx,ky,kz))
        	fnew = fnew + AiiVec(kx,ky,kz)*sin((kx*xcoord*2.0*pi/L+Phii1Vec(kx,ky,kz)))*sin(ky*ycoord*2.0*pi/L + Phii2Vec(kx,ky,kz)); //*Sin(kz*z(k)*2.0*PI/L+Phii3Vec(kx,ky,kz))
        	// fnew = fnew + ruser_mesh_data[3](kx,ky,kz)*sin((kx*xcoord*2.0*pi/L+ruser_mesh_data[0](kx,ky,kz)))*sin(ky*ycoord*2.0*pi/L + ruser_mesh_data[1](kx,ky,kz)); //*Sin(kz*z(k)*2.0*PI/L+Phii3Vec(kx,ky,kz))
            }
          }
        }
        Real density = rhoUp;
	if (multipleClouds == true) {
          Real x1 = xcoord - x_cloud;
          Real x2 = ycoord - y_cloud;
          Real r2 = x1*x1 + x2*x2;
          Real r  = sqrt(r2);

          Real x1_b = xcoord - x_cloud2;
          Real x2_b = ycoord - y_cloud2;
          Real r2_b = x1_b*x1_b + x2_b*x2_b;
          Real r_b  = sqrt(r2_b);

          density = rhoUp + (rho_c - rhoUp) * 0.5 *
                               (1.0 - 1.0*tanh((r-r_cloud)/delta_r));
          // add in 2nd cloud
          density = density + (rho_c - rhoUp) * 0.5 *
                               (1.0 - 1.0*tanh((r_b-r_cloud2)/delta_r));
	}
	else if (turbClumps == true) {
	  density = rhoUp*exp(alpha_exp*fnew);
	}

	// Add in magnetized shock -- Chad, January 27, 2020
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //         ! put in a shock 
  	if (xcoord < shockPos) { // downstream
  		rhoZone = rhoDown;
        	presZone = presDown;
        	velxZone = -velDown;
        	BZone = BDown;
 	}
  	else { //upstream
  	  presZone = presUp;
  	  velxZone = -velUp;
          BZone = B0;
  	  rhoZone = density; // changed to density, already defined above
  	}
          // endif
  /*
	   if(sim_killdivb) then
      		facexData(MAG_FACE_VAR,i,j,k) = B0*sin(BAngle*pi/180.)
      		if(NDIM >= 2) then
        		if (xL(i) > shockPos) faceyData(MAG_FACE_VAR,i,j,k) = B0*cos(BAngle*pi/180.)
        		if (xR(i) < shockPos) faceyData(MAG_FACE_VAR,i,j,k) = BDown*cos(BAngle*pi/180.)
        		if ((xL(i) < shockPos) .AND. (xR(i) > shockPos)) then
            			faceyData(MAG_FACE_VAR,i,j,k) = (BDown+B0)*cos(BAngle*pi/180.)/2.
        		end if
      		end if
      		if(NDIM == 3) facezData(MAG_FACE_VAR,i,j,k) = 0.
	    endif

  */
        phydro->u(IDN,k,j,i) = rhoZone;
        phydro->u(IM1,k,j,i) = rhoZone*velxZone;                             // no initial momentum
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){
         // phydro->u(IEN,k,j,i) = pgas/(gamma-1.0);
          phydro->u(IEN,k,j,i) = presZone/(gamma-1.0) + 0.5*rhoZone*velxZone*velxZone;
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
         //  pcr->u_cr(CRF1,k,j,i) = -velUp*ecr0;                       // no initial CR flux
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
   // Real bx = .0349;
   // Real by = 0.;
   // Real bz = 0.;

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = B0*sin(BAngle*pi/180.);
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            if (pcoord->x1v(i) > shockPos) pfield->b.x2f(k,j,i) = B0*cos(BAngle*pi/180.);
            if (pcoord->x1v(i) < shockPos) pfield->b.x2f(k,j,i) = BDown*cos(BAngle*pi/180.);
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


/*

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
         if(pcr->stream_flag){
           vtot1 += pcr->v_adv(0,k,j,i);
           vtot2 += pcr->v_adv(1,k,j,i);
           vtot3 += pcr->v_adv(2,k,j,i);
         }

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

           // perpendicular energy source term is already added via ec_source_
           // we will still need to include this in general for diffusion case
           vtot2 = 0.0;
           vtot3 = 0.0;

         }

         Real sigma_x = pcr->sigma_diff(0,k,j,i);
         Real sigma_y = pcr->sigma_diff(1,k,j,i);
         Real sigma_z = pcr->sigma_diff(2,k,j,i);

         if(pcr->stream_flag){
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

        Real heat2code = 4.81e23; // Chad
        // Now apply the invert rotation
        if(MAGNETIC_FIELDS_ENABLED){
         // Apply rotation of the vectors
          InvRotateVec(sint_b[i],cost_b[i],sinp_b[i],cosp_b[i],
                                         newfr1,newfr2,newfr3);
          new_ec += dt * ec_source_(k,j,i);
        }
         // Add the energy source term
         if (NON_BAROTROPIC_EOS && (pcr->src_flag > 0)){
           Real new_eg = u(IEN,k,j,i) - (new_ec - ec[i]);
           if(new_eg < 0.0) new_eg = u(IEN,k,j,i);
          // u(IEN,k,j,i) = new_eg;   
          // Chad    
           u(IEN,k,j,i) = new_eg + heat2code*dt*((1.0/6.0)*pcr->dedt_pion(k,j,i) + pcr->dedt_coul(k,j,i));
         }

         if(new_ec < 0.0) new_ec = ec[i];

         // JW: comment next 3 lines out to decouple gas from CR pressure gradients
         if(pcr->src_flag > 0){
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
         pmb->user_out_var(0,k,j,i) = pcr->dedt_coul(k,j,i);
         pmb->user_out_var(1,k,j,i) = pcr->dedt_pion(k,j,i);
         pmb->user_out_var(2,k,j,i) = -pcr->dedt_coul(k,j,i) - pcr->dedt_pion(k,j,i); // CR loss part
         pmb->user_out_var(3,k,j,i) = (1.0/6.0)*pcr->dedt_pion(k,j,i) + pcr->dedt_coul(k,j,i);
        // cout << pcr->dedt_coul(k,j,i) << endl;
        // cout << pmb->user_out_var(0,k,j,i) << endl;

      }// end 

    }// end j
  }// end k
  // Add user defined source term for cosmic rays
  if(pcr->cr_source_defined)
    pcr->UserSourceTerm_(pmb, pmb->pmy_mesh->time, dt, w,bcc, u_cr);

//}

// Chad note:
// Maybe I just need to enroll a CR user source term, and I don't need the full integrator? 
// There doesn't seem to be a way to enroll the integrator function


// Detecting shocks

// from dmr.cpp: ("simple shock finder")
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;
  int k = pmb->ks;
  for (int j=pmb->js; j<=pmb->je; j++) {
    for (int i=pmb->is; i<=pmb->ie; i++) {
      Real epsr = (std::abs(w(IDN,k,j,i+1) - 2.0*w(IDN,k,j,i) + w(IDN,k,j,i-1))
                  + std::abs(w(IDN,k,j+1,i) - 2.0*w(IDN,k,j,i) + w(IDN,k,j-1,i)))
                  /w(IDN,k,j,i);
      Real epsp = (std::abs(w(IPR,k,j,i+1) - 2.0*w(IPR,k,j,i) + w(IPR,k,j,i-1))
                  + std::abs(w(IPR,k,j+1,i) - 2.0*w(IPR,k,j,i) + w(IPR,k,j-1,i)))
                  /w(IPR,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
    }
  }
  // refine : curvature > 0.01
  if (maxeps > 0.01) return 1;
  // derefinement: curvature < 0.005
  if (maxeps < 0.005) return -1;
  // otherwise, stay
  return 0;
}


// Adding CRs where (strong) shocks are detected

                ! Initial setup was a gaussian profile of CR pressure
                ! ecr is specific CR energy
                ! precursorWidth is the width of the initial CR gaussian profile

                pcrZone  = ecr(1) * rho * (sim_gamma_cr - 1.) !CR pressure
                if (CRInjectFraction == .true.) then
                   if (abs(xCenter(i+1)) > 0.0) then
                      dx = xCenter(i+1)-xCenter(i)
                   else
                      dx = xCenter(i)-xCenter(i-1)
                   end if
                   eCRAdded = (eta*rho*velUp**3)*exp(-0.5*((xCenter(i)-peak)**2)/precursorWidth**2)
                   ecrZone = ecr(1) + eCRAdded*dt
                   eintZone = ei + eCRAdded*dt
                else if (CRKEFlux == .true.) then
                  ! if ((xL(i) <= peak) .and. (xR(i) >= peak)) then
                     if (abs(xCenter(i+1)) > 0.0) then
                        dx = xCenter(i+1)-xCenter(i)
                     else
                        dx = xCenter(i)-xCenter(i-1)
                     end if
                     !pCRAdded = dt*(etaInject*rho*velUp**3)/abs(dx)*exp(-0.5*((xCenter(i)-peak)**2)/precursorWidth**2)
                     pCRAdded = rampup*dt*(etaInject*rhoUp*velUp**3)/abs(dx)*exp(-0.5*((xCenter(i)-peak)**2)/precursorWidth**2)
                    ! pCRAdded = dt*eta*0.5*(rhoUp*velUp**3)/abs(dx)
                    ! print *, "PCRAdded: "
                    ! print *, pCRAdded
                    ! print *, "PCRZone: "
                    ! print *, pcrZone

                     pcrZone = pcrZone + pCRAdded
                     ecrZone = pcrZone / (rho * (sim_gamma_cr - 1.))
                     eintZone = ei + pCRAdded/(rho*(sim_gamma_cr - 1.))
                 !  end if  
                else if (CRKeepFraction == .true.) then  ! what I've been using
                   ! Keeps the peak CR pressure at fraction (eta) of flow pressure
                   if ((xL(i) <= peak) .and. (xR(i) >= peak)) then
                      pCRAdded = eta*(rhoUp*velUp**2) * &
                        exp(-0.5*((xCenter(i)-peak)**2)/precursorWidth**2) - &
                        pcrZone
                      pcrZone = eta*(rhoUp*velUp**2) * &
                        exp(-0.5*((xCenter(i)-peak)**2)/precursorWidth**2)
                      ecrZone = pcrZone / (rho * (sim_gamma_cr - 1.))
                      eintZone = ei + pCRAdded/(rho*(sim_gamma_cr - 1.))

                   else
                      eintZone = ei
                      ecrZone = ecr(1)
                   end if

                else if (CRDynamic == .true.) then  ! follows shock
                   if (abs(xCenter(i+1)) > 0.0) then
                      dx = xCenter(i+1)-xCenter(i)
                   else
                      dx = xCenter(i)-xCenter(i-1)
                   end if
                 !  if (abs(solnData(VELX_VAR,i+1,j,k)/solnData(VELX_VAR,i-1,j,k)) > shockGrad) then
                  ! print *, okShock
                   if (okShock) then
                     pCRAdded = rampup*dt*(etaInject*rhoUp*velUp**3)/abs(dx)
                     ! pCRAdded = timeProf*eta*(solnData(DENS_VAR,i+1,j,k)*solnData(VELX_VAR,i+1,j,k)**2) - pcrZone
                     ! pcrZone = timeProf*eta*(solnData(DENS_VAR,i+1,j,k)*solnData(VELX_VAR,i+1,j,k)**2)
                     ! ecrZone = pcrZone / (rho * (sim_gamma_cr - 1.))
                     ! eintZone = ei + pCRAdded/(rho*(sim_gamma_cr - 1.))
                !     print *,"pCRAdded: ", pCRAdded
                     pcrZone = pcrZone + pCRAdded
                     ecrZone = ecr(1) + pCRAdded / (rho * (sim_gamma_cr - 1.))
                    ! ecrZone = pcrZone / (rho * (sim_gamma_cr - 1.))
                 !    print *,"eintZone: ", eintZone
                     eintZone = ei + pCRAdded/(rho*(sim_gamma_cr - 1.))
                  !   print *, ecrZone 
                   else
                      eintZone = ei
                      ecrZone = ecr(1)
                   end if
                else
                   pcrZone = pcrZone + pCRAdded * dt
                   eintZone = ei + pCRAdded*dt/(rho*(sim_gamma_cr - 1.))
                   ecrZone = pcrZone / (rho * (sim_gamma_cr - 1.))
                end if

                ekinZone = 0.5 * (velx**2 + vely**2 + velz**2)
                enerZone = eintZone + ekinZone
                gamcZone = (sim_gamma_gas*pres + sim_gamma_cr*pcrZone) / &
                        (pres + pcrZone)
                gameZone = (sim_gamma_gas*(sim_gamma_cr-1.)*pres + &
                        sim_gamma_cr*(sim_gamma_gas-1.)*pcrZone) / &
                        ((sim_gamma_cr-1.)*pres + (sim_gamma_gas-1.)*pcrZone)




} // new ending curl for adding cosmic rays function

*/





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

/*
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

*/

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

/*
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

*/

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
         // if(n==CRE) {
         //   u_cr(n,k,j,ie+i) = 0.999*u_cr(n,k,j,ie);
         // } else {
         //   u_cr(n,k,j,ie+i) = u_cr(n,k,j,ie);
         // }
	  u_cr(n,k,j,ie+i) = u_cr(n,k,j,ie);
	 // std::cout << "u_cr: " << u_cr(n,k,j,ie+i) << " ";
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
void RadCoolCGS(MeshBlock *pmb, const Real time, const Real dt, 
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
       // Real dcgs = prim(IDN,k,j,i)*dens2cgs;  // density in cgs units
       // Real pcgs = prim(IEN,k,j,i)*pres2cgs;  // gas pressure in cgs units
        Real dcgs = prim(IDN,k,j,i);  // density in cgs units
        Real pcgs = prim(IEN,k,j,i);  // gas pressure in cgs units
                                               // not that the variable IEN in the primitives vector is PRESSURE, not energy density
                                               // so no gamma factor is needed

        Real temp = mu*mp*pcgs/(dcgs*kb);      // temperature in Kelvin
        Real nH = 0.7*dcgs/mp;                 // proton density for fully ionized medium with ~25% Helium by mass

        Real tx = log10(temp)-5.;
        Real theta = 0.4*tx - 3. +5.2/(exp(tx+.08)+exp(-1.5*(tx+.08)));
        Real lambda = 1.1e-21*pow(10.,theta);
        if (temp < Tfloor) { lambda = 0.; }
        //cons(IEN,k,j,i) -= dt*(nH*nH*lambda - nH*gHeat)*heat2code;  // apply cooling
        cons(IEN,k,j,i) -= dt*(nH*nH*lambda - nH*gHeat);  // apply cooling

	pmb->user_out_var(4,k,j,i) = -(nH*nH*lambda - nH*gHeat); // in cgs
      }
    }
  } // end i,j,k
  return;
}

//----------------------------------------------------------------------------------------
//  \brief Sets boundary condition on right X boundary (iib)

void ShockOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                       FaceField &b, Real time, Real dt,
                       int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
 // Real velUp = pin->GetOrAddReal("problem","velUp", 1.E6);
 // Real rhoUp = pin->GetOrAddReal("problem","rhoUp", 1.E-26);
 // Real TUp = pin->GetOrAddReal("problem","TUp", 1.E7);
 // Real presUp = 1.380658E-16 * TUp * rhoUp / (1.6726E-24);  // upstream pressure
 
 // int nc2 = ju-jl;
 // fboundary.NewAthenaArray(ngh,nc2,1);
  Real fbound = 0.0;
  Real sim_xMax = pmb->pcoord->x1v(iu);
 // std::cout << "xmax: " << sim_xMax;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
       // prim(IDN,k,j,il+i) = rhoUp;
        if (turbClumps == true) {
  	  Real ycoord = pmb->pcoord->x2v(j);
        //  std::cout << "x2v(j): " << pmb->pcoord->x2v(j);
      	 // fboundary(i,j,1) = 0.0;
      	  fbound = 0.0;
      //  alpha_exp = 2./40.0
      //  print *, sim_xMax, velUp, L, Phii1Vec(10,10,10)
          for (int kx = sim_kx1; kx<=sim_kx2; ++kx) {
            for (int ky = sim_ky1; ky<=sim_ky2; ++ky) {
              for (int kz = sim_ky1; kz<=sim_ky1; ++kz) { // use the same for y and z
                // print *, AiiVec(kx,ky,kz), Phii1Vec(kx,ky,kz),Phii2Vec(kx,ky,kz),shockPos, L
               // fboundary(i,j,1) = fboundary(i,j,1) + 
                fbound = fbound + 
		  AiiVec(kx,ky,kz)*sin((kx*(sim_xMax + velUp*time - shockPos)*2.0*pi/L +
                  Phii1Vec(kx,ky,kz)))*sin(ky*ycoord*2.0*pi/L +
                  Phii2Vec(kx,ky,kz)); //*Sin(kz*z(k)*2.0*pi/L+Phii3Vec(kx,ky,kz))
	//	std::cout << "fbound" << fbound;
               // fboundary(i,j,1) = fboundary(i,j,1) + 
		//  ruser_mesh_data[3](kx,ky,kz)*sin((kx*(sim_xMax + velUp*time - shockPos)*2.0*pi/L+
                //  ruser_mesh_data[0](kx,ky,kz)))*sin(ky*ycoord*2.0*pi/L +
                //  ruser_mesh_data[1](kx,ky,kz)); //*Sin(kz*z(k)*2.0*pi/L+Phii3Vec(kx,ky,kz))
              }  // using yCoord instead of yCoordf because this is for density, not B fields
	    }
	  }
         // prim(IDN,k,j,il+i) = rhoUp*exp(alpha_exp*fboundary(i,j,1));
          prim(IDN,k,j,iu+i) = rhoUp*exp(alpha_exp*fbound);
	 // std::cout << prim(IDN,k,j,il+i);
	}
	else {
	  prim(IDN,k,j,iu+i) = rhoUp;
	}
	
        prim(IVX,k,j,iu+i) = -velUp;
        prim(IVY,k,j,iu+i) = 0.0;
        prim(IVZ,k,j,iu+i) = 0.0;
        prim(IPR,k,j,iu+i) = presUp;
      }
    }
  }
}
