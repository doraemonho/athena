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
//! \file uniform_chem.cpp
//  \brief problem generator, uniform mesh with chemistry
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
#include "../chem_rad/radiation.hpp"
#include "../scalars/scalars.hpp"

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

Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
                   const Real E, const Real Edot);

Real CoolingTimeStep(MeshBlock *pmb);
Real MyTimeStep(MeshBlock *pmb);

//Radiation boundary
namespace {
  AthenaArray<Real> G0_iang;
  Real G0, cr_rate;
  Real dfloor, pfloor;
  Real cfl_cool_sub, user_dt;
  int nsub_max;
} //namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//========================================================================================

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

  bool is_subcycling = pin->GetOrAddBoolean("chemistry", "subcycling", false);
  if (is_subcycling){
#ifdef CVODE
    std::stringstream msg;
    msg << "### FATAL ERROR in ProblemGenerator::Chem_gow17_MHD" << std::endl
        << "user enable both sub-cycling solver and CVODE solver" << std::endl;
    throw std::runtime_error(msg.str().c_str());
#endif
    EnrollUserTimeStepFunction(CoolingTimeStep);
  }

  user_dt = pin->GetOrAddReal("time", "user_dt", 0.0);
  if ( user_dt > 0.0 && is_subcycling)
    EnrollUserTimeStepFunction(MyTimeStep);

  //EnrollUserBoundaryFunction(BoundaryFace::inner_x1, SixRayBoundaryInnerX1);
  //EnrollUserBoundaryFunction(BoundaryFace::outer_x1, SixRayBoundaryOuterX1);
  //EnrollUserBoundaryFunction(BoundaryFace::inner_x2, SixRayBoundaryInnerX2);
  //EnrollUserBoundaryFunction(BoundaryFace::outer_x2, SixRayBoundaryOuterX2);
  //EnrollUserBoundaryFunction(BoundaryFace::inner_x3, SixRayBoundaryInnerX3);
  //EnrollUserBoundaryFunction(BoundaryFace::outer_x3, SixRayBoundaryOuterX3);
  G0 = pin->GetOrAddReal("chem_radiation", "G0", 1.0);
  G0_iang.NewAthenaArray(6);
  G0_iang(BoundaryFace::inner_x1) = pin->GetOrAddReal("chem_radiation", "G0_inner_x1", G0);
  G0_iang(BoundaryFace::inner_x2) = pin->GetOrAddReal("chem_radiation", "G0_inner_x2", G0);
  G0_iang(BoundaryFace::inner_x3) = pin->GetOrAddReal("chem_radiation", "G0_inner_x3", G0);
  G0_iang(BoundaryFace::outer_x1) = pin->GetOrAddReal("chem_radiation", "G0_outer_x1", G0);
  G0_iang(BoundaryFace::outer_x2) = pin->GetOrAddReal("chem_radiation", "G0_outer_x2", G0);
  G0_iang(BoundaryFace::outer_x3) = pin->GetOrAddReal("chem_radiation", "G0_outer_x3", G0);
  cr_rate = pin->GetOrAddReal("chem_radiation", "CR", 2e-16);
  cfl_cool_sub = pin->GetOrAddReal("chemistry", "cfl_cool_sub", 0.5);
  nsub_max = pin->GetOrAddInteger("chemistry","nsub_max",1e5);
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief initialize problem by reading in vtk file.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
  //read input parameters
  const Real nH = pin->GetReal("problem", "nH"); //density
  const Real vx = pin->GetOrAddReal("problem", "vx", 0);
  const Real vy = pin->GetOrAddReal("problem", "vy", 0);
  const Real vz = pin->GetOrAddReal("problem", "vz", 0);
  const Real b0 = pin->GetOrAddReal("problem", "b0", 0);
  const Real angle = (PI/180.0)*pin->GetOrAddReal("problem","angle",0.0);
  const Real G0 = pin->GetOrAddReal("problem", "G0", 0.);

  const Real s_init = pin->GetOrAddReal("problem", "s_init", 0.);
  //mean and std of the initial gaussian profile
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = nH*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;
  Real float_min = std::numeric_limits<float>::min();
  pfloor = pin->GetOrAddReal("hydro", "pfloor",(1024*(float_min)));  
  dfloor = pin->GetOrAddReal("hydro", "dfloor",(1024*(float_min)));

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        //density
        phydro->u(IDN, k, j, i) = nH;
        //velocity, x, y, z direction
        phydro->u(IM1, k, j, i) = nH*vx;
        phydro->u(IM2, k, j, i) = nH*vy;
        phydro->u(IM3, k, j, i) = nH*vz;
        //energy
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*nH*SQR(vx) + 0.5*nH*SQR(vy) +0.5*nH*SQR(vz);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          phydro->u(IEN,k,j,i)+=0.5*b0*b0;
        }
      }
    }
  }

  //intialize magnetic field field
  if (MAGNETIC_FIELDS_ENABLED) {
    if (COORDINATE_SYSTEM == "cartesian") {

        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = 0;
        }}}
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = b0 * std::cos(angle);
        }}}
        for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = b0 * std::sin(angle);
        }}}
    }else{
      std::stringstream msg;
    msg << "### FATAL ERROR in ProblemGenerator::Chem_gow17_MHD" << std::endl
        << "Only support cartesian system!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    }
  }

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
  if (NSCALARS > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSCALARS; ++ispec) {
            pscalars->s(ispec, k, j, i) = s_init * phydro->u(IDN, k, j, i);
#ifdef INCLUDE_CHEMISTRY
            Real s_ispec = pin->GetOrAddReal("problem",
                "s_init_"+pscalars->chemnet.species_names[ispec], -1);
            if (s_ispec >= 0.) {
              pscalars->s(ispec, k, j, i) = s_ispec * phydro->u(IDN, k, j, i);
            }
#endif
          }
        }
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkInLoop(ParameterInput *pin)
//  \brief
//========================================================================================
void MeshBlock::UserWorkInLoop() {
  const Real unit_length_in_cm_  = 3.086e+18;
  const Real unit_vel_in_cms_    = 1.0e5;
  const Real unit_density_in_nH_ = 1;
  const Real unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  const Real unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  const Real  g =  peos->GetGamma();
  const Real Tfloor = 20.0;
  //set density and pressure floors
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
          Real& u_d  = phydro->u(IDN,k,j,i);
          Real& w_p  = phydro->w(IPR,k,j,i);
          Real& u_e  = phydro->u(IEN,k,j,i);
          Real& u_m1 = phydro->u(IM1,k,j,i);
          Real& u_m2 = phydro->u(IM2,k,j,i);
          Real& u_m3 = phydro->u(IM3,k,j,i);
          
          Real   nH_  = u_d*unit_density_in_nH_;
          Real   ED   = w_p/(g-1.0);
          Real E_ergs = ED * unit_E_in_cgs_ / nH_;
          Real     T  =  E_ergs / (1.5*1.381e-16);

          Real pfloor = Tfloor* (1.5*1.381e-16) * nH_/unit_E_in_cgs_*(g - 1.0);
          w_p = (T > Tfloor) ?  w_p : pfloor;
          Real di = 1.0/u_d;
          Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));

            
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          u_e = w_p/(g-1.0)+ke;
#else  // MHD:
          Real me =0.5*0.25*(SQR(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))
                           + SQR(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))
                           + SQR(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)));
         u_e = w_p/(g-1.0)+ke+me;
#endif
          
      }
    }
  }
  return;
}


void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                           FaceField &b, Real time, Real dt,
                           int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //set species and column boundary to zero
  for (int n=0; n<(NSCALARS); ++n) {
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
  for (int n=0; n<(NSCALARS); ++n) {
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
  for (int n=0; n<(NSCALARS); ++n) {
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
  for (int n=0; n<(NSCALARS); ++n) {
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
  for (int n=0; n<(NSCALARS); ++n) {
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
  for (int n=0; n<(NSCALARS); ++n) {
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

Real CoolingTimeStep(MeshBlock *pmb){

  AthenaArray<Real> &u = pmb->phydro->u;
  const Real yfloor = 1e-3;
  const Real unit_length_in_cm_  = 3.086e+18;
  const Real unit_vel_in_cms_    = 1.0e5;
  const Real unit_density_in_nH_ = 1;
  const Real unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  const Real unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  const Real  g = 5.0/3.0;

  Real    E = 0;
  Real Edot = 0;
  Real    y[NSPECIES];
  Real ydot[NSPECIES];
  
  //MHD timestep
  Real min_dt = pmb->pmy_mesh->dt;
  //code time
  Real time = pmb->pmy_mesh->time;
  Real tsub = 0.0;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {

        E = pmb->phydro->w(IPR,k,j,i)/(g-1.0);

        pmb->pscalars->chemnet.InitializeNextStep(k, j, i);
        //copy species abundance
        for (int ispec=0; ispec<=NSPECIES; ispec++) {
          y[ispec] = pmb->pscalars->s(ispec,k,j,i)/u(IDN,k,j,i);
        }
        //calculate reaction rates
        pmb->pscalars->chemnet.RHS(time, y, E, ydot);
        //calculate heating and cooling rats
        if (NON_BAROTROPIC_EOS) {
          Edot = pmb->pscalars->chemnet.Edot(time, y, E);
        }
        //get the sub-cycle dt 
        tsub = nsub_max * cfl_cool_sub * GetChemTime(y, ydot, E, Edot);
        // manually choosen timestep 1e-3 dt to stablize the system at the begining of the simulation
        tsub = std::min(1.5e-3, tsub);
        min_dt = std::min(tsub, min_dt);

      }
    }
  }

  return min_dt;
}

//----------------------------------------------------------------------------------------
//! \fn Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
//                       const Real E, const Real Edot)
//! \brief calculate chemistry timescale
Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
                 const Real E, const Real Edot) {
  const Real small_ = 1024 * std::numeric_limits<float>::min();
  //put floor in species abundance
  Real yf[NSPECIES];
  Real yfloor = 1e-3;
  for (int ispec=0; ispec<=NSPECIES; ispec++) {
    yf[ispec] = std::max(y[ispec], yfloor);
  }
  //calculate chemistry timescale
  Real tchem = std::abs( yf[0]/(ydot[0] + small_) );
  for (int ispec=1; ispec<NSPECIES; ispec++) {
    tchem = std::min( tchem, std::abs(yf[ispec]/(ydot[ispec]+small_)) );
  } 
  if (NON_BAROTROPIC_EOS) {
    tchem = std::min( tchem, std::abs(E/(Edot+small_)) );
  }
  return tchem;
}
//----------------------------------------------------------------------------------------
//! \fn MyTimeStep(MeshBlock *pmb)
//! \brief Setup self-defined dt
Real MyTimeStep(MeshBlock *pmb)
{
  Real min_user_dt = user_dt;
  return min_user_dt;
}