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
//! \file chem_SimpleMultiPhase.cpp
//  \brief problem generator, uniform mesh with simple synthetic coooling function for multiphase ISM simulation
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

Real CoolingTimeStep(MeshBlock *pmb);
Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
                   const Real E, const Real Edot);

//Radiation boundary
namespace {
  Real cfl_cool_sub;
  Real d_floor;
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

//EnrollUserTimeStepFunction(CoolingTimeStep);  
cfl_cool_sub = pin->GetOrAddReal("chemistry", "cfl_cool_sub", 0.3);
nsub_max = pin->GetOrAddInteger("chemistry","nsub_max",1e5);
d_floor = pin->GetOrAddInteger("chemistry","d_floor",1e-4);

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
  const Real vx = pin->GetOrAddReal("problem", "vx", 0); //velocity x
  const Real vy = pin->GetOrAddReal("problem", "vy", 0); //velocity x
  const Real vz = pin->GetOrAddReal("problem", "vz", 0); //velocity x
  const Real b0 = pin->GetOrAddReal("problem","b0",0.0); //velocity x
  const Real angle = (PI/180.0)*pin->GetOrAddReal("problem","angle",0.0);
  const Real G0 = pin->GetOrAddReal("problem", "G0", 0.);

  const Real s_init = pin->GetOrAddReal("problem", "s_init", 0.);
  //mean and std of the initial gaussian profile
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = nH*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;

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
          phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*nH*SQR(vx);
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
    msg << "### FATAL ERROR in ProblemGenerator::chem_SimpleMultiPhase" << std::endl
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
              pchemrad->ir(k, j, i, ifreq * pchemrad->nang + iang) = G0;
            }
          }
        }
      }
    }
  }

  //intialize chemical species
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
            pscalars->s(0, k, j, i) =   1.0*phydro->u(IDN, k, j, i);
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
  
  Real min_dt = pmb->pmy_mesh->dt; //MHD timestep
  Real time = pmb->pmy_mesh->time;//code time
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
        tsub = 1e-1 * nsub_max * cfl_cool_sub * GetChemTime(y, ydot, E, Edot);
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

void MeshBlock::UserWorkInLoop() {
  const Real unit_length_in_cm_  = 3.086e+18;
  const Real unit_vel_in_cms_    = 1.0e5;
  const Real unit_density_in_nH_ = 1;
  const Real unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  const Real unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  const Real  g =  peos->GetGamma();
  const Real Tfloor = 10.0;
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
          
          // Check if u_d < d_floor
          u_d = (u_d > d_floor) ?  u_d : d_floor;
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
