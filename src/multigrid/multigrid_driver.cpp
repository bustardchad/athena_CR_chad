//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid_driver.cpp
//  \brief implementation of functions in class MultigridDriver

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "multigrid.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// constructor, initializes data structures and parameters

MultigridDriver::MultigridDriver(Mesh *pm, MeshBlock *pmb, MGBoundaryFunc_t *MGBoundary,
                                 ParameterInput *pin)
{
  pmy_mesh_=pm;
  if(pblock_->block_size.nx1!=pblock_->block_size.nx2
  || pblock_->block_size.nx1!=pblock_->block_size.nx3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The Multigrid solver requires logically cubic MeshBlocks." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pblock_->block_size.nx2==1 || pblock_->block_size.nx3==1 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Currently the Multigrid solver works only in 3D." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  if(pm->use_meshgen_fn_[X1DIR]==true || pm->use_meshgen_fn_[X2DIR]==true
  || pm->use_meshgen_fn_[X3DIR]==true) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "Non-uniform mesh spacing is not supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  Real dx=pblock_->pcoord->dx1f(0);
  if(dx!=pblock_->pcoord->dx2f(0) || dx!=pblock_->pcoord->dx3f(0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The cell size must be cubic." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  // count multigrid levels
  nmblevel_=0;
  for(int l=0; l<20; l++) {
    if((1<<l) == pblock_->block_size.nx1) {
      nmblevel_=l+1;
      break;
    }
  }
  if(nmblevel_==0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The MeshBlock size must be power of two." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  // count multigrid levels
  nrlx=0, nrly=0, nrlz=0;
  for(int l=0; l<20; l++) {
    if((1<<l) == pm->nrbx1)
      nrlx=l+1;
    if((1<<l) == pm->nrbx2)
      nrly=l+1;
    if((1<<l) == pm->nrbx3)
      nrlz=l+1;
  }
  if(nrlx==0 || nrly==0 || nrlz==0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MultigridDriver::MultigridDriver" << std::endl
        << "The root grid size must be power of 2." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    break;
  }
  nrootlevel_=std::min(nrlx,std::min(nrly,nrlz));
  rootsrc_.NewAthenaArray(nvar_,pm->nrbx3,pm->nrbx2,pm->nrbx1);
  ntotallevel_=nrootlevel_+nmblevel_-1;
  current_level_=ntotallevel_-1;

  fperiodic_=false;
  for(int i=0; i<6; i++)
    MGBoundaryFunction_[i]=pm->MGBoundaryFunction_[i];
  if(MGBoundaryFunction_[INNER_X1]==MGPeriodicInnerX1
  && MGBoundaryFunction_[OUTER_X1]==MGPeriodicOuterX1
  && MGBoundaryFunction_[INNER_X2]==MGPeriodicInnerX2
  && MGBoundaryFunction_[OUTER_X2]==MGPeriodicOuterX2
  && MGBoundaryFunction_[INNER_X3]==MGPeriodicInnerX3
  && MGBoundaryFunction_[OUTER_X3]==MGPeriodicOuterX3)
    fperiodic_=true;

  mode_=0; // 0: FMG+V(1,1), 1: FMG+F(0,1), 2: V(1,1)

  mgroot_=NULL;

  // Setting up the MPI information
  // *** this part should be modified when dedicate processes are allocated ***
  pblock_=pmb;
  MeshBlocks *pb=pblock_;
  nblocks_=0;
  while(pb!=NULL) {
    nblocks_++;
    pb=pb->next;
  }

  nranks_  = Globals::nranks;
  nslist_  = new int[nranks_];
  nblist_  = new int[nranks_];
  nvlist_  = new int[nranks_];
  nvslist_ = new int[nranks_];
#ifdef MPI_PARALLEL
  MPI_Comm_dup(MPI_COMM_WORKD, &MPI_COMM_MULTIGRID);
#endif
  for(int n=0; n<nranks_; n++) {
    nslist_[n]  = pm->nslist[n];
    nblist_[n]  = pm->nblist[n];
    nvlist_[n]  = nblist_[n]*nvar_;
    nvslist_[n] = nslist_[n]*nvar_;
  }
  rootbuf_=new Real[pm->nbtotal*nvar_];
}

// destructor

MultigridDriver::~MultigridDriver()
{
  delete mgroot_;
  delete [] nslist_;
  delete [] nblist_;
  delete [] nvlist_;
  delete [] nvslist_;
  delete [] rootbuf_;
  rootsrc_.DeleteAthenaArray();
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SetupMultigrid(void)
//  \brief initialize the source assuming that the source terms are already loaded

void MultigridDriver::SetupMultigrid(void)
{
  Mesh *pm=pmy_mesh_;
  MeshBlocks *pb=pblock_;

  if(fperiodic_ || mode_<=1) { // periodic or FMG
    while(pb!=NULL) {
      Multigrid *pmg=GetMultigridBlock(pb);
      for(int v=0; v<nvar_; v++)
        rootbuf_[pb->gid*nvar_+v]=pmg->CalculateTotal(v);
      pb=pb->next;
    }
#ifdef MPI_PARALLEL
    MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                   rootbuf_, nvlist, nvslist, MPI_ATHENA_REAL, MPI_COMM_GRAVITY);
#endif
  }
  if(fperiodic_) { // periodic - subtract average
    Real vol=(pm->mesh_size.(x1max-x1min)*pm->mesh_size.(x2max-x2min)
                                         *pm->mesh_size.(x3max-x3min));
    for(int v=0; v<nvar_; v++) {
      Real total=0.0;
      for(int n=0; n<Globals::nranks; n++)
        total+=rootbuf_[n*nvar_+v];
      ave=total/vol;
      for(int n=0; n<pm->nbtotal; n++)
        rootbuf_[n*nvar_+v]-=ave;
      pb=pblock_;
      while(pb!=NULL) {
        Multigrid *pmg=GetMultigridBlock(pb);
        pmg->SubtractAverage(0, v, ave);
        pb=pb->next;
      }
    }
  }
  if(mode_<=1) { // FMG
    pb=pblock_;
    while(pb!=NULL) {
      Multigrid *pmg=GetMultigridBlock(pb);
      pmg->RestrictFMGSource();
      pb=pb->next;
    }
    FillRootGridSource();
    mgroot_->RestrictFMGSource();
    current_level_=0;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::CollectSource(void)
//  \brief collect the coarsest data and store in the rootbuf_ array

void MultigridDriver::CollectSource(void)
{
  MeshBlocks *pb=pblock_;
  while(pb!=NULL) {
    Multigrid *pmg=GetMultigridBlock(pb);
    for(int v=0; v<nvar_; v++)
      rootbuf_[pb->gid*nvar_+v]=pmg->GetRootSource(v);
    pb=pb->next;
  }
#ifdef MPI_PARALLEL
  MPI_Allgatherv(MPI_IN_PLACE, nblist_[Globals::my_rank]*nvar_, MPI_ATHENA_REAL,
                 rootbuf_, nvlist, nvslist, MPI_ATHENA_REAL, MPI_COMM_GRAVITY);
#endif
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FillRootGrid(void)
//  \brief Fill the root grid using the rootbuf_ array

void MultigridDriver::FillRootGrid(void)
{
  if(pmy_mesh_->multilevel) {
    // *** implement later
  }
  else { // uniform
    for(int n=0; n<pmy_mesh_->nbtotal; n++) {
      LogicalLocation &loc=pmy_mesh_->loclist[n];
      for(int v=0; v<nvar_; v++) {
        rootsrc_(v,loc.lx3,loc.lx2,loc.lx1);
      }
    }
    mgroot_->LoadSource(rootsrc_,0,0,1.0);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::FMGProlongate(void)
//  \brief Prolongation for FMG Cycle

void MultigridDriver::FMGProlongate(void)
{
  if(current_level_ >= nrootlevel_) {
    if(current_level_==nrootlevel_)
      TransferFromRootToBlocks(true);
    else {
      mgtlist_->SetMGTaskListFMGProlongate();
      mgtlist_->DoTaskList(this);
    }
  }
  else { // root grid
    for(int n=0; n<nsmooth; n++) {
      mgroot_->ApplyPhysicalBoundary();
      mgroot_->FMGProlongate();
    }
  }
  current_level++;
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::TransferFromRootToBlock(bool fmgflag)
//  \brief Transfer the data from the root grid to the coarsest level of each MeshBlock

void MultigridDriver::TransferFromRootToBlock(bool fmgflag)
{
  MeshBlock *pb=pblock_;
  AthenaArray<Real> &src=mgroot_->GetCurrentData();
  if(pmy_mesh_->multilevel) {
    // *** implement later ***
  }
  else {
    while(pb!=NULL) {
      Multigrid *pmg=GetMultigridBlock(pb);
      LogicalLocation &loc=pb->loc;
      pmg_->SetFromRootGrid(src, loc.lx1, loc.lx2, loc.lx3);
      if(fmgflag) pmg->FMGProlongate();
      else pmg->ProlongateAndCorrect();
      pb=pb->next;
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepUpward(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToFiner(int nsmooth)
{
  if(current_level_ >= nrootlevel_) {
    mgtlist_->SetMGTaskListToFiner(nsmooth);
    mgtlist_->DoTaskList(this);
    if(current_level_==nrootlevel_) {
      CollectSource();
      FillRootGrid();
    }
  }
  else { // root grid
    mgroot_->ApplyPhysicalBoundary();
    mgroot_->ProlongateAndCorrect();
    for(int n=0; n<nsmooth; n++) {
      mgroot_->ApplyPhysicalBoundary();
      mgroot_->Smooth(0);
      mgroot_->ApplyPhysicalBoundary();
      mgroot_->Smooth(1);
    }
  }
  current_level--;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::OneStepToCoarser(int nsmooth)
//  \brief smoothing and restriction one level

void MultigridDriver::OneStepToCoarser(int nsmooth)
{
  if(current_level_ >= nrootlevel_) {
    if(current_level_==nrootlevel_)
      TransferFromRootToBlocks(0);
    else {
      mgtlist_->SetMGTaskListToCoarser(nsmooth);
      mgtlist_->DoTaskList(this);
    }
  }
  else { // root grid
    for(int n=0; n<nsmooth; n++) {
      mgroot_->ApplyPhysicalBoundary();
      mgroot_->Smooth(0);
      mgroot_->ApplyPhysicalBoundary();
      mgroot_->Smooth(1);
    }
    Restrict();
  }
  current_level++;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveVCycle(int npresmooth, int npostsmooth)
//  \brief Solve the V-cycle starting from the current level

void MultigridDriver::SolveVCycle(void)
{
  int startlevel=current_level_;
  while(current_level_>0)
    OneStepToFiner(npresmooth);
  SolveCoarsestGrid();
  while(current_level_<=startlevel)
    OneStepToFiner(npresmooth);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFCycle(int npresmooth, int npostsmooth)
//  \brief Solve the F-cycle starting from the current level

void MultigridDriver::SolveFCycle(void)
{
  int startlevel=current_level_;
  while(current_level_>0)
    OneStepToFiner(npresmooth);
  SolveCoarsestGrid();
  while(current_level_<=startlevel)
    OneStepToFiner(npresmooth);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveFMGCycle(void)
//  \brief Solve the FMG Cycle

void MultigridDriver::SolveFMGCycle(void)
{
  for(int lev=0; lev<=ntotal_level_; lev++) {
    if(mode_==0)
      SolveVCycle();
    else if(mode_==1)
      SolveFCycle();
    FMGProlongate();
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::SolveCoarsestGrid(void)
//  \brief Solve the coarsest root grid

virtual void MultigridDriver::SolveCoarsestGrid(void)
{
  int ni=std::max(pmy_mesh->nrbx1, std::max(pmy_mesh->nrbx2, pmy_mesh->nrbx3));
  if(fperiodic_ && ni==1) // trivial case - all zero
    return;
  else {
    for(int i=0; i<ni; ni++) { // iterate ni times
      mgroot_->ApplyPhysicalBoundary();
      mgroot_->Smooth(0);
      mgroot_->ApplyPhysicalBoundary();
      mgroot_->Smooth(1);
    }
    if(fperiodic_) {
      Real vol=(pm->mesh_size.(x1max-x1min)*pm->mesh_size.(x2max-x2min)
                                           *pm->mesh_size.(x3max-x3min));
      for(int v=0; v<nvar_; v++) {
        Real ave=pmg->CalculateTotal(v)/vol;
        pmg->SubtractAverage(1, v, ave);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn Real MultigridDriver::CalculateDefectNorm(int n, int nrm)
//  \brief calculate the defect norm

Real MultigridDriver::CalculateDefectNorm(int n, int nrm)
{
  MeshBlocks *pb=pblock_;
  Real norm=0.0;
  while(pb!=NULL) {
    Multigrid *pmg=GetMultigridBlock(pb);
    if(nrm==0)
      norm=std::max(norm, pmg->CauclateDefectNorm(n, nrm));
    else
      norm+=pmg->CauclateDefectNorm(n, nrm);
    pb=pb->next;
  }
#ifdef MPI_PARALLEL
  if(nrm==0)
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_MAX,MPI_COMM_GRAVITY);
  else
    MPI_Allreduce(MPI_IN_PLACE,&norm,1,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_GRAVITY);
#endif
  if(nrm==2)
    norm=std::sqrt(norm);

  return norm;
}

