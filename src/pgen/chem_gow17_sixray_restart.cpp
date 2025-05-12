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
//! \file chem_gow17_sixray_postprocessing.cpp
//  \brief problem generator, uniform mesh with chemistry
//  \python3 configure.py --prob=chem_gow17_sixray_postprocessing --chemistry=gow17 --chem_radiation=six_ray --ode_solver=cvode -fft -hdf5 -mpi --hdf5_path=/opt/cray/pe/hdf5-parallel/1.12.2.3/gnu/9.1/ --fftw_path=/opt/cray/pe/fftw/3.3.10.5/x86_milan --cvode_path=//global/homes/k/kho33/local/cvode
//======================================================================================
// c headers
#include <stdio.h>    // c style file
#include <string.h>   // strcmp()

// C++ headers
#include <algorithm>  // std::find()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../scalars/scalars.hpp"
#include "../inputs/hdf5_reader.hpp"
#include "../outputs/io_wrapper.hpp"

//User defined boundary conditions
void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
  AthenaArray<Real> G0_iang;
  Real G0, cr_rate;
} //namespace


void Mesh::InitUserMeshData(ParameterInput *pin) {

  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
#endif
  }

  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SixRayBoundaryInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SixRayBoundaryOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) { 
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, SixRayBoundaryInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, SixRayBoundaryOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, SixRayBoundaryInnerX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, SixRayBoundaryOuterX3);
  }


  G0 = pin->GetOrAddReal("chem_radiation", "G0", 0.);
  G0_iang.NewAthenaArray(6);
  G0_iang(BoundaryFace::inner_x1) = pin->GetOrAddReal("chem_radiation","G0_inner_x1",G0);
  G0_iang(BoundaryFace::inner_x2) = pin->GetOrAddReal("chem_radiation","G0_inner_x2",G0);
  G0_iang(BoundaryFace::inner_x3) = pin->GetOrAddReal("chem_radiation","G0_inner_x3",G0);
  G0_iang(BoundaryFace::outer_x1) = pin->GetOrAddReal("chem_radiation","G0_outer_x1",G0);
  G0_iang(BoundaryFace::outer_x2) = pin->GetOrAddReal("chem_radiation","G0_outer_x2",G0);
  G0_iang(BoundaryFace::outer_x3) = pin->GetOrAddReal("chem_radiation","G0_outer_x3",G0);
  cr_rate = pin->GetOrAddReal("chem_radiation", "CR", 2e-16);
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem by reading in vtk file.
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //dimensions of meshblock
  const int Nx = ie - is + 1; // sub-block size for each process
  const int Ny = je - js + 1; // sub-block size for each process
  const int Nz = ke - ks + 1; // sub-block size for each process
  //dimensions of mesh
  const int Nx_mesh = pmy_mesh->mesh_size.nx1; // total mesh size 
  const int Ny_mesh = pmy_mesh->mesh_size.nx2; // total mesh size 
  const int Nz_mesh = pmy_mesh->mesh_size.nx3; // total mesh size 
  
  //gamma-1 for hydro eos
  const Real gm1 = peos->GetGamma() - 1.0;
  //initial abundance
  const Real r_init = pin->GetReal("problem", "r_init");
  // get the user defined abundance flag
  bool user_define_abundance = pin->GetOrAddBoolean("problem", "self_define_abundance", false);
  //parse input parameters
  std::string input_filename = pin->GetString("problem", "hdffile");
  std::string dataset_prims  = pin->GetString("problem", "dataset");
  // read the hdf5 index of each phys vars
  int index_dens  = pin->GetInteger("problem", "index_dens");
  int index_vel1  = pin->GetInteger("problem", "index_vel1");
  int index_vel2  = pin->GetInteger("problem", "index_vel2");
  int index_vel3  = pin->GetInteger("problem", "index_vel3");
  int index_pres  = pin->GetInteger("problem", "index_pres");

  AthenaArray<Real> &bcc = pfield->bcc;
  FaceField &b = pfield->b;

if (MAGNETIC_FIELDS_ENABLED) {
  Real b0 = pin->GetOrAddReal("problem", "b0", 1.0);
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = 0;
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = 0.0;
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = b0;
  }}}
  pfield->CalculateCellCenteredField(b, bcc, pcoord, is,ie,js,je,ks,ke);
}

#ifdef DEBUG
  printf("Unjoined vtk file. data size = (%d, %d, %d)\n",
          Nz, Ny, Nx);
#endif

  int rank_file = 5;
  int rank_mem  = 3;

  int start_mem[3] = {0,       0,    0}; 
  int count_mem[3] = {Nx-1, Ny-1, Nz-1};

  bool collective=false;
  bool noop=false;

  //dagnostic printing of filename
  if (Globals::my_rank == 0) {
    printf("meshblock gid=%d, lx1=%lld, lx2=%lld, lx3=%lld, level=%d\n",
            gid, loc.lx1, loc.lx2, loc.lx3, loc.level);
  }

#ifdef DEBUG
if (Globals::my_rank == 0) {
  printf("Process 0: start to reading.\n");
}
#endif

  // Set conserved array selections
  int start_prims_file[5];
  start_prims_file[1] = gid;
  start_prims_file[2] = 0;
  start_prims_file[3] = 0;
  start_prims_file[4] = 0;
  int start_prims_indices[5];
  start_prims_indices[IDN] = index_dens;
  start_prims_indices[IPR] = index_vel1;
  start_prims_indices[IVX] = index_vel2;
  start_prims_indices[IVY] = index_vel3;
  start_prims_indices[IVZ] = index_pres;

  int count_prims_file[5];
  count_prims_file[0] = 1;
  count_prims_file[1] = 1;
  count_prims_file[2] = block_size.nx3;
  count_prims_file[3] = block_size.nx2;
  count_prims_file[4] = block_size.nx1;
  int start_prims_mem[4];
  start_prims_mem[1] = ks;
  start_prims_mem[2] = js;
  start_prims_mem[3] = is;
  int count_prims_mem[4];
  count_prims_mem[0] = 1;
  count_prims_mem[1] = block_size.nx3;
  count_prims_mem[2] = block_size.nx2;
  count_prims_mem[3] = block_size.nx1;
    
  // Set conserved values from file
  for (int n = 0; n < NHYDRO; ++n) {
    start_prims_file[0] = start_prims_indices[n];
    start_prims_mem[0]  = n; 
    HDF5ReadRealArray(input_filename.c_str(), dataset_prims.c_str(),
                      5, start_prims_file, count_prims_file,
                      4, start_prims_mem, count_prims_mem,
                      phydro->w, true, true);
  }
  //change primative variables to conservative variables.
  peos->PrimitiveToConserved(phydro->w, bcc, phydro->u, pcoord,
                                     is, ie, js, je, ks, ke);

#ifdef DEBUG
if (Globals::my_rank == 0) {
  printf("Process 0: Starting Additional reading if MPI is enabled.\n");
}
#endif

  // Make no-op collective reads if using MPI and ranks have unequal numbers of blocks
#ifdef MPI_PARALLEL
  {
    int num_blocks_this_rank = pmy_mesh->nblist[Globals::my_rank];
    if (lid == num_blocks_this_rank - 1) {
      int block_shortage_this_rank = 0;
      for (int rank = 0; rank < Globals::nranks; ++rank) {
        block_shortage_this_rank =
            std::max(block_shortage_this_rank,
                     pmy_mesh->nblist[rank] - num_blocks_this_rank);
      }
      for (int block = 0; block < block_shortage_this_rank; ++block) {
        for (int n = 0; n < NHYDRO; ++n) {
          start_prims_file[0] = n;
          start_prims_mem[0]  = n; 
          HDF5ReadRealArray(input_filename.c_str(), dataset_prims.c_str(), 5,
                            start_prims_file, count_prims_file, 4,
                            start_prims_mem, count_prims_mem,
                            phydro->w, true, true);
          //change primative variables to conservative variables.
          peos->PrimitiveToConserved(phydro->w, bcc, phydro->u, pcoord,
                                            is, ie, js, je, ks, ke);
        }
      }
    }
  }
#endif

#ifdef DEBUG
if (Globals::my_rank == 0) {
  printf("Process 0: Finished All Reading.\n");
}
#endif
  //intialize radiation field
  if (CHEMRADIATION_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifreq=0; ifreq < pchemrad->nfreq; ++ifreq) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              pchemrad->ir(k, j, i, ifreq * pchemrad->nang + iang) = G0_iang(iang);
            }
          }
#ifdef INCLUDE_CHEMISTRY
          for (int iang=0; iang < pchemrad->nang; ++iang) {
            //cr rate
            pchemrad->ir(k, j, i,
                pscalars->chemnet.index_cr_ * pchemrad->nang + iang) = cr_rate;
          }
#endif
        }
      }
    }
    //calculate the average radiation field for output of the initial condition
    pchemrad->pchemradintegrator->CopyToOutput();
  }

  //intialize chemical species
  if (NSPECIES > 0) {
    // inital abundance by existing file
    if (user_define_abundance){ 
      std::string dataset_chem  = "chemistry";
      for (int ispec=0; ispec < NSPECIES; ++ispec) {
        start_prims_file[0] = ispec; //index of species
        start_prims_mem[0]  = ispec; 
        HDF5ReadRealArray(input_filename.c_str(), dataset_chem.c_str(),
                          5, start_prims_file, count_prims_file,
                          4, start_prims_mem, count_prims_mem,
                          pscalars->s, true, true);
      }

#ifdef MPI_PARALLEL
    {
      int num_blocks_this_rank = pmy_mesh->nblist[Globals::my_rank];
      if (lid == num_blocks_this_rank - 1) {
        int block_shortage_this_rank = 0;
        for (int rank = 0; rank < Globals::nranks; ++rank) {
          block_shortage_this_rank =
              std::max(block_shortage_this_rank,
                      pmy_mesh->nblist[rank] - num_blocks_this_rank);
        }
        for (int block = 0; block < block_shortage_this_rank; ++block) {
          for (int ispec=0; ispec < NSPECIES; ++ispec) {
            start_prims_file[0] = ispec; //index of species
            start_prims_mem[0]  = ispec; 
            HDF5ReadRealArray(input_filename.c_str(), dataset_chem.c_str(), 5,
                              start_prims_file, count_prims_file, 4,
                              start_prims_mem, count_prims_mem,
                              pscalars->s, true, true);
          }
        }
      }
    }
#endif

    }else{
      // inital abundance by hand
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            for (int ispec=0; ispec < NSPECIES; ++ispec) {
              pscalars->s(ispec, k, j, i) = r_init * phydro->u(IDN, k, j, i);
  #ifdef INCLUDE_CHEMISTRY
              Real s_ispec = pin->GetOrAddReal("problem",
                  "r_init_"+pscalars->chemnet.species_names[ispec], -1);
              if (s_ispec >= 0.) {
                pscalars->s(ispec, k, j, i) = s_ispec * phydro->u(IDN, k, j, i);
              }
  #endif
            }
          }
        }
      }
    }
  }

  Real tinit = 0.0;
  Real dt = 100.0;
  // evole the chemistry and energy to steady state
  pscalars->odew.Initialize(pin);
  pscalars->odew.Integrate(tinit, dt);
  // sync. the primitive varibles
  peos->ConservedToPrimitive(phydro->u, phydro->w, pfield->b,
                             phydro->w, pfield->bcc,
                             pcoord, is, ie, js, je, ks, ke); 

  return;
}

void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=1; i<=ngh; ++i) {
          pmb->pscalars->s(n,k,j,il-i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,il-i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,k,jl-j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,jl-j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,kl-k,j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          prim(n,kl-k,j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=1; i<=ngh; ++i) {
          pmb->pscalars->s(n,k,j,iu+i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,iu+i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,k,ju+j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          prim(n,k,ju+j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          pmb->pscalars->s(n,ku+k,j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma simd
        for (int i=il; i<=iu; ++i) {
          prim(n,ku+k,j,i) = 0;
        }
      }
    }
  }
  return;
}