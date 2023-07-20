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
#include "../fft/turbulence.hpp"

Real CoolingTimeStep(MeshBlock *pmb);
Real GetChemTime(const Real y[NSPECIES], const Real ydot[NSPECIES],
                   const Real E, const Real Edot);

//Radiation boundary
namespace {
  Real cfl_cool_sub;
  Real d_floor;
  Real v_max;
  Real Tmax;
  int nsub_max;
  int sign(Real number);

  // User defined function for history output
  Real B2overRho(MeshBlock *pmb, int iout);
  Real absrho2(MeshBlock *pmb, int iout);
  Real rhoudota_OU(MeshBlock *pmb, int iout);
  Real curlU2(MeshBlock *pmb, int iout);
  Real abspdivV(MeshBlock *pmb, int iout);
  Real absdivV2(MeshBlock *pmb, int iout);
  Real absdivV(MeshBlock *pmb, int iout);
  Real absdivB(MeshBlock *pmb, int iout);

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
d_floor = pin->GetOrAddInteger("hydro","dfloor",1e-4);
v_max = pin->GetOrAddInteger("chemistry","v_max",100.0);
Tmax = pin->GetOrAddInteger("chemistry","Tmax",5e4);

// User defined History Output
AllocateUserHistoryOutput(8);
EnrollUserHistoryOutput(0, absdivB,   "<|∇⋅B|>");
EnrollUserHistoryOutput(1, B2overRho, "<|B2/ρ|>");
EnrollUserHistoryOutput(2, absrho2,   "<|ρ2|>");
EnrollUserHistoryOutput(3, absdivV,   "<|∇⋅V|>");
EnrollUserHistoryOutput(4, absdivV2,  "<|∇⋅V|2>");
EnrollUserHistoryOutput(5, abspdivV,  "<|p∇⋅V|>");
EnrollUserHistoryOutput(6, curlU2,    "<|∇XV|2>");
EnrollUserHistoryOutput(7, rhoudota_OU, "<|ρu⋅a|>");

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
          


          // Check if u_d < d_floor && if v > vmax 
          if (u_d  < d_floor){
            Real u_1,u_2,u_3;
            u_1 = u_m1/ u_d;
            u_2 = u_m2/ u_d;
            u_3 = u_m3/ u_d;
            u_d  = d_floor;
            u_m1 = u_d*u_1;
            u_m2 = u_d*u_2;
            u_m3 = u_d*u_3;
          }

          Real   nH_  = u_d*unit_density_in_nH_;
          Real   ED   = w_p/(g-1.0);
          Real E_ergs = ED * unit_E_in_cgs_ / nH_;
          Real     T  =  E_ergs / (1.5*1.381e-16);

          Real pfloor = Tfloor* (1.5*1.381e-16) * nH_/unit_E_in_cgs_*(g - 1.0);
          Real pmax   =   Tmax* (1.5*1.381e-16) * nH_/unit_E_in_cgs_*(g - 1.0);
          w_p = (T > Tfloor) ?  w_p : pfloor;
          w_p = (T > Tmax) ?   pmax : w_p;
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
//user defined function for the problem generator
namespace {
int sign(Real number) {
  return (number > 0) - (number < 0);
}

//2. $<\nabla\cdot \vec{v}>$ -> OK
//3. $<p\nabla\cdot\vec{v}>$ -> OK
//4. $<|\nabla\cdot\vec{v}|^2>$ -> OK 
//5. $<|\nabla\times\vec{v}|^2>$  -> OK
//6. $<|\nabla\cdot\vec{B}|>$  -> OK
//7. $<\rho u\dot a> = <\rho u\dot dv_OU/dt>$ -> OK
//8. $<B^2/\rho>$ -> OK
//9. $<\rho^2>$ -> OK

Real absdivB(MeshBlock *pmb, int iout) {

  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1, face2p, face2m, face3p, face3m;
  FaceField &b = pmb->pfield->b;
  Real absdivb = 0.0;

  face1.NewAthenaArray( (ie-is)+2*NGHOST+2);
  face2p.NewAthenaArray((ie-is)+2*NGHOST+1);
  face2m.NewAthenaArray((ie-is)+2*NGHOST+1);
  face3p.NewAthenaArray((ie-is)+2*NGHOST+1);
  face3m.NewAthenaArray((ie-is)+2*NGHOST+1);

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      pmb->pcoord->Face1Area(k,   j,   is, ie+1, face1);
      pmb->pcoord->Face2Area(k,   j+1, is, ie,   face2p);
      pmb->pcoord->Face2Area(k,   j,   is, ie,   face2m);
      pmb->pcoord->Face3Area(k+1, j,   is, ie,   face3p);
      pmb->pcoord->Face3Area(k,   j,   is, ie,   face3m);
#pragma omp simd
      for(int i=is; i<=ie; i++) {
        absdivb += std::abs((face1(i+1)*b.x1f(k,j,i+1)-face1(i)*b.x1f(k,j,i)
                            +face2p(i)*b.x2f(k,j+1,i)-face2m(i)*b.x2f(k,j,i)
                            +face3p(i)*b.x3f(k+1,j,i)-face3m(i)*b.x3f(k,j,i)));
      }
    }
  }

  return absdivb/N;
}

Real absrho2(MeshBlock *pmb, int iout) {
  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  AthenaArray<Real> vol;
  vol.NewAthenaArray((ie-is)+2*NGHOST+1);
  Real absrho2 = 0.0;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
#pragma omp simd
      for(int i=is; i<=ie; i++) {
        absrho2+= vol(i)*pow(pmb->phydro->u(IDN,k,j,i),2);
      }
    }
  }

  return absrho2/N;
}

Real B2overRho(MeshBlock *pmb, int iout) {

  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  FaceField &b = pmb->pfield->b;
  Real b2overrho = 0.0;

  AthenaArray<Real> vol;
  vol.NewAthenaArray((ie-is)+2*NGHOST+1);

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
#pragma omp simd
      for(int i=is; i<=ie; i++) {          
        Real b2 = 0.5*( SQR(b.x1f(k,j,i) + b.x1f(k,j,i+1))
                      + SQR(b.x3f(k,j,i) + b.x3f(k+1,j,i))
                      + SQR(b.x2f(k,j,i) + b.x2f(k,j+1,i)));
        b2overrho += vol(i)*b2/pmb->phydro->u(IDN,k,j,i);
      }
    }
  }

  return b2overrho/N;
}


Real absdivV(MeshBlock *pmb, int iout) {
  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  int il = is - 1; int iu = ie + 1;
  int jl = js - 1; int ju = je + 1;
  int kl = ks - 1; int ku = ke + 1;

  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> x1area_, x2area_, x3area_, x2area_p1_, x3area_p1_, vol_;

  x1area_.NewAthenaArray((ie-is)+2*NGHOST+2);
  x2area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x2area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  vol_.NewAthenaArray((ie-is)+2*NGHOST+1);

  AthenaArray<Real> tmp_div;
  tmp_div.NewAthenaArray((ie-is)+2*NGHOST+1, (je-js)+2*NGHOST+1, (ke-ks)+2*NGHOST+1);
  
  Real area_p1, area;
  Real vel_p1, vel;
  Real div_vel = 0.0;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, il, iu+1, x1area_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(w(IVX,k,j,i+1) + w(IVX,k,j,i  ));
        vel     = 0.5*(w(IVX,k,j,i  ) + w(IVX,k,j,i-1));
        tmp_div(i,j,k) = area_p1*vel_p1 - area*vel;
      }
      // calculate x2-flux divergnece
      pmb->pcoord->Face2Area(k, j  , il, iu, x2area_);
      pmb->pcoord->Face2Area(k, j+1, il, iu, x2area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x2area_p1_(i);
        area    = x2area_(i);
        vel_p1  = 0.5*(w(IVY,k,j+1,i) + w(IVY,k,j  ,i));
        vel     = 0.5*(w(IVY,k,j  ,i) + w(IVY,k,j-1,i));
        tmp_div(i,j,k) += area_p1*vel_p1 - area*vel;
      }
      pmb->pcoord->Face3Area(k  , j, il, iu, x3area_);
      pmb->pcoord->Face3Area(k+1, j, il, iu, x3area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x3area_p1_(i);
        area    = x3area_(i);
        vel_p1  = 0.5*(w(IVZ,k+1,j,i) + w(IVZ, k  ,j,i));
        vel     = 0.5*(w(IVZ,k  ,j,i) + w(IVZ, k-1,j,i));
        tmp_div(i,j,k) += area_p1*vel_p1 - area*vel;
      }
    }
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->CellVolume(k,j,il,iu,vol_);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        div_vel += tmp_div(i,j,k)/vol_(i);
      }
    }
  }

  return div_vel/N;
}

// Computing |div V|^2
Real absdivV2(MeshBlock *pmb, int iout) {
  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  int il = is - 1; int iu = ie + 1;
  int jl = js - 1; int ju = je + 1;
  int kl = ks - 1; int ku = ke + 1;

  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> x1area_, x2area_, x3area_, x2area_p1_, x3area_p1_, vol_;

  x1area_.NewAthenaArray((ie-is)+2*NGHOST+2);
  x2area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x2area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  vol_.NewAthenaArray((ie-is)+2*NGHOST+1);

  AthenaArray<Real> tmp_div;
  tmp_div.NewAthenaArray((ie-is)+2*NGHOST+1, (je-js)+2*NGHOST+1, (ke-ks)+2*NGHOST+1);
  
  Real area_p1, area;
  Real vel_p1, vel;
  Real div_vel2 = 0.0;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, il, iu+1, x1area_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(w(IVX,k,j,i+1) + w(IVX,k,j,i  ));
        vel     = 0.5*(w(IVX,k,j,i  ) + w(IVX,k,j,i-1));
        tmp_div(i,j,k) = area_p1*vel_p1 - area*vel;
      }
      // calculate x2-flux divergnece
      pmb->pcoord->Face2Area(k, j  , il, iu, x2area_);
      pmb->pcoord->Face2Area(k, j+1, il, iu, x2area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x2area_p1_(i);
        area    = x2area_(i);
        vel_p1  = 0.5*(w(IVY,k,j+1,i) + w(IVY,k,j  ,i));
        vel     = 0.5*(w(IVY,k,j  ,i) + w(IVY,k,j-1,i));
        tmp_div(i,j,k) += area_p1*vel_p1 - area*vel;
      }
      pmb->pcoord->Face3Area(k  , j, il, iu, x3area_);
      pmb->pcoord->Face3Area(k+1, j, il, iu, x3area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x3area_p1_(i);
        area    = x3area_(i);
        vel_p1  = 0.5*(w(IVZ,k+1,j,i) + w(IVZ, k  ,j,i));
        vel     = 0.5*(w(IVZ,k  ,j,i) + w(IVZ, k-1,j,i));
        tmp_div(i,j,k) += area_p1*vel_p1 - area*vel;
      }
    }
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->CellVolume(k,j,il,iu,vol_);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        div_vel2 += std::pow(tmp_div(i,j,k),2)/vol_(i);
      }
    }
  }

  return div_vel2/N;
}

//Computing |p div V|
Real abspdivV(MeshBlock *pmb, int iout) {
  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  int il = is - 1; int iu = ie + 1;
  int jl = js - 1; int ju = je + 1;
  int kl = ks - 1; int ku = ke + 1;

  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> x1area_, x2area_, x3area_, x2area_p1_, x3area_p1_, vol_;

  x1area_.NewAthenaArray((ie-is)+2*NGHOST+2);
  x2area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x2area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  vol_.NewAthenaArray((ie-is)+2*NGHOST+1);

  AthenaArray<Real> tmp_div;
  tmp_div.NewAthenaArray((ie-is)+2*NGHOST+1, (je-js)+2*NGHOST+1, (ke-ks)+2*NGHOST+1);
  
  Real area_p1, area;
  Real vel_p1, vel;
  Real p_div_vel = 0.0;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, il, iu+1, x1area_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(w(IVX,k,j,i+1) + w(IVX,k,j,i  ));
        vel     = 0.5*(w(IVX,k,j,i  ) + w(IVX,k,j,i-1));
        tmp_div(i,j,k) = area_p1*vel_p1 - area*vel;
      }
      // calculate x2-flux divergnece
      pmb->pcoord->Face2Area(k, j  , il, iu, x2area_);
      pmb->pcoord->Face2Area(k, j+1, il, iu, x2area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x2area_p1_(i);
        area    = x2area_(i);
        vel_p1  = 0.5*(w(IVY,k,j+1,i) + w(IVY,k,j  ,i));
        vel     = 0.5*(w(IVY,k,j  ,i) + w(IVY,k,j-1,i));
        tmp_div(i,j,k) += area_p1*vel_p1 - area*vel;
      }
      pmb->pcoord->Face3Area(k  , j, il, iu, x3area_);
      pmb->pcoord->Face3Area(k+1, j, il, iu, x3area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x3area_p1_(i);
        area    = x3area_(i);
        vel_p1  = 0.5*(w(IVZ,k+1,j,i) + w(IVZ, k  ,j,i));
        vel     = 0.5*(w(IVZ,k  ,j,i) + w(IVZ, k-1,j,i));
        tmp_div(i,j,k) += area_p1*vel_p1 - area*vel;
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
        pmb->pcoord->CellVolume(k,j,il,iu,vol_);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        p_div_vel += w(IPR,k,j,i)*tmp_div(i,j,k)/vol_(i);
      }
    }
  }

  return p_div_vel/N;
}
// Computing |Curl U|^2
Real curlU2(MeshBlock *pmb, int iout) {
  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
 
  int il = is - 1; int iu = ie + 1;
  int jl = js - 1; int ju = je + 1;
  int kl = ks - 1; int ku = ke + 1;

  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> x1area_, x2area_, x3area_, x2area_p1_, x3area_p1_, vol_;

  x1area_.NewAthenaArray((ie-is)+2*NGHOST+2);
  x2area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x2area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  x3area_p1_.NewAthenaArray((ie-is)+2*NGHOST+1);
  vol_.NewAthenaArray((ie-is)+2*NGHOST+1);

  AthenaArray<Real> tmp_curlx1, tmp_curlx2, tmp_curlx3;
  tmp_curlx1.NewAthenaArray((ie-is)+2*NGHOST+1, (je-js)+2*NGHOST+1, (ke-ks)+2*NGHOST+1);
  tmp_curlx2.NewAthenaArray((ie-is)+2*NGHOST+1, (je-js)+2*NGHOST+1, (ke-ks)+2*NGHOST+1);
  tmp_curlx3.NewAthenaArray((ie-is)+2*NGHOST+1, (je-js)+2*NGHOST+1, (ke-ks)+2*NGHOST+1);

  Real area_p1, area;
  Real vel_p1, vel;
  Real abs_curl2 = 0.0;

// Update curl
// x1 dir : d/dx2 (v3) - d/dx3 (v2)
// x2 dir : d/dx3 (v1) - d/dx1 (v3)
// x3 dir : d/dx1 (v2) - d/dx2 (v1)
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // x1->x3 face :   d/dx1 (v2)
      // x1->x2 face : - d/dx1 (v3)
      pmb->pcoord->Face1Area(k, j, il, iu+1, x1area_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        area_p1 = x1area_(i+1);
        area    = x1area_(i);
        vel_p1  = 0.5*(w(IVY,k,j,i+1) + w(IVY,k,j,i  ));
        vel     = 0.5*(w(IVY,k,j,i  ) + w(IVY,k,j,i-1));
        tmp_curlx3(i,j,k) += area_p1*vel_p1 - area*vel;
    
        vel_p1  = 0.5*(w(IVZ,k,j,i+1) + w(IVZ,k,j,i  ));
        vel     = 0.5*(w(IVZ,k,j,i  ) + w(IVZ,k,j,i-1));
        tmp_curlx2(i,j,k) -= (area_p1*vel_p1 - area*vel);
      }

      pmb->pcoord->Face2Area(k, j  , il, iu, x2area_);
      pmb->pcoord->Face2Area(k, j+1, il, iu, x2area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        // x2->x1 dir :    d/dx2 (v3)
        // x2->x3 dir :  - d/dx2 (v1)
        area_p1 = x2area_p1_(i);
        area    = x2area_(i);
        vel_p1  = 0.5*(w(IVZ,k,j+1,i) + w(IVZ,k,j  ,i));
        vel     = 0.5*(w(IVZ,k,j  ,i) + w(IVZ,k,j-1,i));
        tmp_curlx1(i,j,k) += area_p1*vel_p1 - area*vel;

        vel_p1  = 0.5*(w(IVX,k,j+1,i) + w(IVX,k,j  ,i));
        vel     = 0.5*(w(IVX,k,j  ,i) + w(IVX,k,j-1,i));
        tmp_curlx3(i,j,k) -= area_p1*vel_p1 - area*vel;
      }
      pmb->pcoord->Face3Area(k  , j, il, iu, x3area_);
      pmb->pcoord->Face3Area(k+1, j, il, iu, x3area_p1_);
#pragma omp simd private(area_p1, area, vel_p1, vel)
      for (int i=il; i<=iu; ++i) {
        // x1 dir : - d/dx3 (v2)
        // x2 dir :   d/dx3 (v1)
        area_p1 = x3area_p1_(i);
        area    = x3area_(i);
        vel_p1  = 0.5*(w(IVY,k+1,j,i) + w(IVY, k  ,j,i));
        vel     = 0.5*(w(IVY,k  ,j,i) + w(IVY, k-1,j,i));
        tmp_curlx1(i,j,k) -= area_p1*vel_p1 - area*vel;

        vel_p1  = 0.5*(w(IVX,k+1,j,i) + w(IVX, k  ,j,i));
        vel     = 0.5*(w(IVX,k  ,j,i) + w(IVX, k-1,j,i));
        tmp_curlx2(i,j,k) += area_p1*vel_p1 - area*vel;
      }
      pmb->pcoord->CellVolume(k,j,il,iu,vol_);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real curl2 = tmp_curlx1(i,j,k)*tmp_curlx1(i,j,k) + 
                     tmp_curlx2(i,j,k)*tmp_curlx2(i,j,k) + 
                     tmp_curlx3(i,j,k)*tmp_curlx3(i,j,k);
        abs_curl2 += curl2/vol_(i);
      }
    }
  }
  return abs_curl2/N;
}

//7. $<\rho u\dot a> = <\rho u\dot dv_OU/dt>$ 
Real rhoudota_OU(MeshBlock *pmb, int iout) {
  Real N = pmb->pmy_mesh->mesh_size.nx1 * pmb->pmy_mesh->mesh_size.nx2 * pmb->pmy_mesh->mesh_size.nx3;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;

  Real rho_u_dot_aOU = 0.0;
  Real dt = pmb->pmy_mesh->dt;
  Mesh *pmesh = pmb->pmy_mesh;
  TurbulenceDriver *ptrbd = pmesh->ptrbd;
  AthenaArray<Real> &dv1 = ptrbd->vel[0], &dv2 = ptrbd->vel[1], &dv3 = ptrbd->vel[2];
  
  AthenaArray<Real> vol;
  vol.NewAthenaArray((ie-is)+2*NGHOST+1);

  for (int nb=0; nb<pmesh->nblocal; ++nb) {
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);  
#pragma omp simd
        for(int i=is; i<=ie; i++) {
          Real rho = pmb->phydro->u(IDN,k,j,i);
          Real  ux = pmb->phydro->u(IVX,k,j,i);
          Real  uy = pmb->phydro->u(IVY,k,j,i);
          Real  uz = pmb->phydro->u(IVZ,k,j,i);
          Real dvx = dv1(nb,k,j,i);
          Real dvy = dv2(nb,k,j,i);
          Real dvz = dv3(nb,k,j,i);
          rho_u_dot_aOU += vol(i)*rho*(ux*dvx + uy*dvy + uz*dvz);
        }
      }
    }
  }
  return rho_u_dot_aOU/dt/N;
}

} //namespace