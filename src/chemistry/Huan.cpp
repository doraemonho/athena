//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file Huan.cpp
//! \brief implementation of the Huan method solver (improved Euler with 2 order accuracy)

//c header
#include <stdio.h> //c style io

//c++ header
#include <ctime> //time
#include <iostream>   // endl, ostream
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept> //throw exceptions
#include <string>

// Athena++ classes headers
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

// this class header
#include "ode_wrapper.hpp" // NOLINT

namespace {
  Real cfl_cool_sub; //cfl number for subcycling
  Real yfloor; //species abundance floor for calculating the cooling time
  int nsub_max; //maximum number of substeps
  Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
                   const Real E, const Real Edot);
  void IntegrateHalfSubstep(Real tsub, Real y[NSPECIES], const Real ydot[NSPECIES],
                           Real &E, const Real Edot);
  void IntegrateFullSubstep(Real tsub, 
                             Real y[NSPECIES], const Real ydot0[NSPECIES], const Real ydot1[NSPECIES],
                             Real &E, const Real Edot0, const Real Edot1);
} //namespace

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper constructor
ODEWrapper::ODEWrapper(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block_ = pmb;
  if (NON_BAROTROPIC_EOS) {
    dim_ = NSPECIES + 1;
  } else {
    dim_ = NSPECIES;
  }
  output_zone_sec_ = pin->GetOrAddBoolean("chemistry", "output_zone_sec", false);
  cfl_cool_sub = pin->GetOrAddReal("chemistry","cfl_cool_sub",0.1);
  yfloor = pin->GetOrAddReal("chemistry","yfloor",1e-3);
  nsub_max = pin->GetOrAddInteger("chemistry","nsub_max",1e5);
}

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper destructor
ODEWrapper::~ODEWrapper() {
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Initialize(ParameterInput *pin)
//! \brief Initialize ODE solver parameters
void ODEWrapper::Initialize(ParameterInput *pin) {
  //Note: this cannot be in the constructor, since it needs the PassiveScalars
  //class, and the ODEWrapper class is constructed in the PassiveScalars constructor.
  pmy_spec_ = pmy_block_->pscalars;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Integrate(const Real tinit, const Real dt)
//! \brief Integrate the ODE forward for time dt
void ODEWrapper::Integrate(const Real tinit, const Real dt) {
  int is = pmy_block_->is;
  int js = pmy_block_->js;
  int ks = pmy_block_->ks;
  int ie = pmy_block_->ie;
  int je = pmy_block_->je;
  int ke = pmy_block_->ke;
  int ncycle = pmy_block_->pmy_mesh->ncycle;
  clock_t tstart, tstop;
  tstart = std::clock();
  const Real scalar_floor = pmy_block_->peos->GetScalarFloor();
  //primitive conserved variables
  AthenaArray<Real> &u = pmy_block_->phydro->u;
  AthenaArray<Real> &bcc = pmy_block_->pfield->bcc;
  //chemical species and rates
  Real  y[NSPECIES];
  Real y1[NSPECIES];
  Real ydot0[NSPECIES];
  Real ydot1[NSPECIES];
  //internal energy and rates
  Real  E = 0.; 
  Real E1 = 0.;
  Real Edot0 = 0.;
  Real Edot1 = 0.;
  Real time = pmy_block_->pmy_mesh->time;//code time
  Real dt_mhd = pmy_block_->pmy_mesh->dt; //MHD timestep
  //subcycling variables
  Real tend, tsub, tnow, tleft;
  Real icount;
  //loop over all cells
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        pmy_spec_->chemnet.InitializeNextStep(k, j, i);
        //copy species abundance
        for (int ispec=0; ispec<=NSPECIES; ispec++) {
          y[ispec] = y1[ispec] = pmy_spec_->s(ispec,k,j,i)/u(IDN,k,j,i);
        }
        //assign internal energy, if not isothermal eos
        if (NON_BAROTROPIC_EOS) {
          E = u(IEN,k,j,i)
            - 0.5*( SQR(u(IM1,k,j,i)) + SQR(u(IM2,k,j,i)) + SQR(u(IM3,k,j,i))
                   )/u(IDN,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            E -= 0.5*(
                SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)) );
          }
          E1 = E;
        }
        //subcycling
        icount = 0;
        tnow = time;
        tend = time + dt_mhd;
        tleft = dt_mhd;
        while (tnow < tend) {
          
          // half step calcution using forward Euler method
          //calculate reaction rates
          pmy_spec_->chemnet.RHS(time, y, E, ydot0);
          //calculate heating and cooling rats
          if (NON_BAROTROPIC_EOS) {
            Edot0 = pmy_spec_->chemnet.Edot(time, y, E);
          }
          //get the sub-cycle dt 
          tsub = cfl_cool_sub * GetChemTime(y, ydot0, E, Edot0);
          tsub = std::min(tsub, tleft);
          
          //advance half step
          IntegrateHalfSubstep(tsub, y1, ydot0, E1, Edot0);

          //Full step calcuation
          //calculate reaction rates
          pmy_spec_->chemnet.RHS(time, y1, E1, ydot1);
          //calculate heating and cooling rats
          if (NON_BAROTROPIC_EOS) {
            Edot1 = pmy_spec_->chemnet.Edot(time, y1, E1);
          }

          // advance the full step
          IntegrateFullSubstep(tsub, y, ydot0, ydot1, E, Edot0, Edot1);

          //update timing
          tnow += tsub;
          tleft = tend - tnow;
          icount++;

          //check maximum number of steps
          if (icount > nsub_max) {
            std::stringstream msg;
            msg << "### FATAL ERROR in function ODEWrapper::Integrate: "
              << "Maximum number of substeps = " << nsub_max 
              << ", tnow = "  << tnow << ", tleft = "  << tleft << ", tsub = " << tsub
              << " exceeded for Huan solver." << std::endl;
            ATHENA_ERROR(msg);
          }
        }
        //copy species abundance back to s
        for (int ispec=0; ispec<=NSPECIES; ispec++) {
          //apply floor to passive scalar concentrations
          y[ispec] = (y[ispec] < scalar_floor) ?  scalar_floor : y[ispec];
          pmy_spec_->s(ispec,k,j,i) = y[ispec]*u(IDN,k,j,i);
        }
        //assign internal energy, if not isothermal eos
        if (NON_BAROTROPIC_EOS) {
          u(IEN,k,j,i) = E
            + 0.5*( SQR(u(IM1,k,j,i)) + SQR(u(IM2,k,j,i)) + SQR(u(IM3,k,j,i))
                )/u(IDN,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            u(IEN,k,j,i) += 0.5*(
                SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)) );
          }
        }
      }
    }
  }
  tstop = std::clock();
  if (output_zone_sec_) {
    double cpu_time = (tstop>tstart ? static_cast<double> (tstop-tstart) :
                       1.0)/static_cast<double> (CLOCKS_PER_SEC);
    std::uint64_t nzones =
      static_cast<std::uint64_t> (pmy_block_->GetNumberOfMeshBlockCells());
    double zone_sec = static_cast<double> (nzones) / cpu_time;
    printf("chemistry ODE integration: ");
    printf("ncycle = %d, total time in sec = %.2e, zone/sec=%.2e\n",
        ncycle, cpu_time, Real(nzones)/cpu_time);
  }
  return;
}

namespace {

//----------------------------------------------------------------------------------------
//! \fn Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
//                       const Real E, const Real Edot)
//! \brief calculate chemistry timescale
Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
                 const Real E, const Real Edot) {
  const Real small_ = 1024 * std::numeric_limits<float>::min();
  //put floor in species abundance
  Real yf[NSPECIES];
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
//! \fn void IntegrateHalfSubstep(Real tsub, Real y[NSPECIES], const Real ydot[NSPECIES],
//!                              Real &E, const Real Edot)
//! \brief advance the half-step of chemical abundance and energy for tsub
void IntegrateHalfSubstep(Real tsub, Real y[NSPECIES], const Real ydot[NSPECIES],
                         Real &E, const Real Edot) {
  for (int ispec=0; ispec<NSPECIES; ispec++) {
    y[ispec] += ydot[ispec] * tsub;
  }
  if (NON_BAROTROPIC_EOS) {
    E += Edot * tsub;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void IntegrateFullSubstep(Real tsub, Real y[NSPECIES], 
//!                               const Real ydot0[NSPECIES], const Real ydot1[NSPECIES],
//!                               Real &E, const Real Edot0, const Real Edot1)
//! \brief advance the full step of chemical abundance and energy for tsub
void IntegrateFullSubstep(Real tsub, 
                          Real y[NSPECIES], const Real ydot0[NSPECIES], const Real ydot1[NSPECIES],
                          Real &E, const Real Edot0, const Real Edot1) {
  for (int ispec=0; ispec<NSPECIES; ispec++) {
    y[ispec] += (ydot0[ispec] + ydot1[ispec]) * tsub * 0.5;
  }
  if (NON_BAROTROPIC_EOS) {
    E += ( Edot0 + Edot1 ) * tsub * 0.5;
  }
  return;
}

} //namespace

