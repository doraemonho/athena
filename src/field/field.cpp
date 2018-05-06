//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field.cpp
//  \brief implementation of functions in class Field

// C++ headers
#include <string>

// Athena++ headers
#include "field.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

// constructor, initializes data structures and parameters

Field::Field(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  // Allocate memory for interface fields, but only when needed.
  if (MAGNETIC_FIELDS_ENABLED) {
    int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
    if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

    //  Note the extra cell in each longitudinal dirn for interface fields
    b.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    b1.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b1.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b1.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    // Allocated 3rd registers if user-requested time integrator is type 3N or 3S*
    std::string integrator = pin->GetOrAddString("time","integrator","vl2");
    if (integrator == "ssprk5_4") {
      // future extension may add "int nregister" to Hydro class
      b2.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
      b2.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
      b2.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    }

    bcc.NewAthenaArray(NFIELD, ncells3, ncells2, ncells1);

    //-------- begin allocations for fourth-order MHD
    // TODO(kfelker): add xorder==4 switch for array allocation
    b_fc.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b_fc.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b_fc.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    bcc_center.NewAthenaArray(NFIELD, ncells3, ncells2, ncells1);

    // fourth-order UCT reconstructions at corners
    // TODO(kfelker): cut down on these temporary arrays
    // TODO(kfelker): check array limits-- add +1 for upper face?
    by_W.NewAthenaArray(ncells3, ncells2, ncells1);
    by_E.NewAthenaArray(ncells3, ncells2, ncells1);
    bx_S.NewAthenaArray(ncells3, ncells2, ncells1);
    bx_N.NewAthenaArray(ncells3, ncells2, ncells1);
    // 3D states
    bz_R1.NewAthenaArray(ncells3, ncells2, ncells1);
    bz_L1.NewAthenaArray(ncells3, ncells2, ncells1);
    bz_R2.NewAthenaArray(ncells3, ncells2, ncells1);
    bz_L2.NewAthenaArray(ncells3, ncells2, ncells1);
    by_R3.NewAthenaArray(ncells3, ncells2, ncells1);
    by_L3.NewAthenaArray(ncells3, ncells2, ncells1);
    bx_R3.NewAthenaArray(ncells3, ncells2, ncells1);
    bx_L3.NewAthenaArray(ncells3, ncells2, ncells1);

    // TODO(kfelker): expand to 3D
    v_NE.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_SE.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_NW.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_SW.NewAthenaArray(2, ncells3, ncells2, ncells1);
    // 3D states
    v_R3R2.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_R3L2.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_L3R2.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_L3L2.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_R3R1.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_R3L1.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_L3R1.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_L3L1.NewAthenaArray(2, ncells3, ncells2, ncells1);

    vl_temp_.NewAthenaArray(2, ncells3, ncells2, ncells1);
    vr_temp_.NewAthenaArray(2, ncells3, ncells2, ncells1);
    alpha_plus_x1_.NewAthenaArray(ncells3, ncells2, ncells1);
    alpha_minus_x1_.NewAthenaArray(ncells3, ncells2, ncells1);
    alpha_plus_x2_.NewAthenaArray(ncells3, ncells2, ncells1);
    alpha_minus_x2_.NewAthenaArray(ncells3, ncells2, ncells1);
    alpha_plus_x3_.NewAthenaArray(ncells3, ncells2, ncells1);
    alpha_minus_x3_.NewAthenaArray(ncells3, ncells2, ncells1);
    //-------- end allocations for fourth-order MHD

    e.x1e.NewAthenaArray((ncells3+1),(ncells2+1), ncells1   );
    e.x2e.NewAthenaArray((ncells3+1), ncells2   ,(ncells1+1));
    e.x3e.NewAthenaArray( ncells3   ,(ncells2+1),(ncells1+1));

    wght.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    wght.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    wght.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    e2_x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    e3_x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    e1_x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    e3_x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    e1_x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    e2_x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    // Allocate memory for scratch vectors
    cc_e_.NewAthenaArray(ncells3,ncells2,ncells1);

    face_area_.NewAthenaArray(ncells1);
    edge_length_.NewAthenaArray(ncells1);
    edge_length_p1_.NewAthenaArray(ncells1);
    if (GENERAL_RELATIVITY) {
      g_.NewAthenaArray(NMETRIC,ncells1);
      gi_.NewAthenaArray(NMETRIC,ncells1);
    }

    // fourth-order MHD
    // 4D scratch arrays
    // TODO(kfelker): these could all share the same 4D array, extended
    // by 1 in all directions
    scr1_nkji_cc_.NewAthenaArray(NFIELD, ncells3, ncells2, ncells1);
    scr1_kji_x1fc_.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    scr2_kji_x2fc_.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    scr3_kji_x3fc_.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
  }
}

// destructor

Field::~Field() {
  b.x1f.DeleteAthenaArray();
  b.x2f.DeleteAthenaArray();
  b.x3f.DeleteAthenaArray();
  b1.x1f.DeleteAthenaArray();
  b1.x2f.DeleteAthenaArray();
  b1.x3f.DeleteAthenaArray();
  // b2 only allocated if integrator was 3S* integrator
  b2.x1f.DeleteAthenaArray();
  b2.x2f.DeleteAthenaArray();
  b2.x3f.DeleteAthenaArray();
  bcc.DeleteAthenaArray();

  //-------- begin deallocations for fourth-order MHD
  b_fc.x1f.DeleteAthenaArray();
  b_fc.x2f.DeleteAthenaArray();
  b_fc.x3f.DeleteAthenaArray();
  bcc_center.DeleteAthenaArray();

  by_W.DeleteAthenaArray();
  by_E.DeleteAthenaArray();
  bx_N.DeleteAthenaArray();
  bx_S.DeleteAthenaArray();
  // 3D states
  bz_R1.DeleteAthenaArray();
  bz_L1.DeleteAthenaArray();
  bz_R2.DeleteAthenaArray();
  bz_L2.DeleteAthenaArray();
  by_R3.DeleteAthenaArray();
  by_L3.DeleteAthenaArray();
  bx_R3.DeleteAthenaArray();
  bx_L3.DeleteAthenaArray();

  v_NE.DeleteAthenaArray();
  v_SE.DeleteAthenaArray();
  v_NW.DeleteAthenaArray();
  v_SW.DeleteAthenaArray();
  // 3D states
  v_R3R2.DeleteAthenaArray();
  v_R3L2.DeleteAthenaArray();
  v_L3R2.DeleteAthenaArray();
  v_L3L2.DeleteAthenaArray();
  v_R3R1.DeleteAthenaArray();
  v_R3L1.DeleteAthenaArray();
  v_L3R1.DeleteAthenaArray();
  v_L3L1.DeleteAthenaArray();

  vl_temp_.DeleteAthenaArray();
  vr_temp_.DeleteAthenaArray();
  alpha_plus_x1_.DeleteAthenaArray();
  alpha_minus_x1_.DeleteAthenaArray();
  alpha_plus_x2_.DeleteAthenaArray();
  alpha_minus_x2_.DeleteAthenaArray();
  // 3D states
  alpha_plus_x3_.DeleteAthenaArray();
  alpha_minus_x3_.DeleteAthenaArray();
  // 4D scratch arrays
  scr1_nkji_cc_.DeleteAthenaArray();
  scr1_kji_x1fc_.DeleteAthenaArray();
  scr2_kji_x2fc_.DeleteAthenaArray();
  scr3_kji_x3fc_.DeleteAthenaArray();
  //----------end deallocations for fourth-order MHD

  e.x1e.DeleteAthenaArray();
  e.x2e.DeleteAthenaArray();
  e.x3e.DeleteAthenaArray();
  wght.x1f.DeleteAthenaArray();
  wght.x2f.DeleteAthenaArray();
  wght.x3f.DeleteAthenaArray();
  e2_x1f.DeleteAthenaArray();
  e3_x1f.DeleteAthenaArray();
  e1_x2f.DeleteAthenaArray();
  e3_x2f.DeleteAthenaArray();
  e1_x3f.DeleteAthenaArray();
  e2_x3f.DeleteAthenaArray();

  cc_e_.DeleteAthenaArray();
  face_area_.DeleteAthenaArray();
  edge_length_.DeleteAthenaArray();
  edge_length_p1_.DeleteAthenaArray();
  if (GENERAL_RELATIVITY) {
    g_.DeleteAthenaArray();
    gi_.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief

void Field::CalculateCellCenteredField(const FaceField &bf, AthenaArray<Real> &bc,
            Coordinates *pco, int is, int ie, int js, int je, int ks, int ke) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    // calc cell centered fields first
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        const Real& b1_i   = bf.x1f(k,j,i  );
        const Real& b1_ip1 = bf.x1f(k,j,i+1);
        const Real& b2_j   = bf.x2f(k,j  ,i);
        const Real& b2_jp1 = bf.x2f(k,j+1,i);
        const Real& b3_k   = bf.x3f(k  ,j,i);
        const Real& b3_kp1 = bf.x3f(k+1,j,i);

        Real& bcc1 = bc(IB1,k,j,i);
        Real& bcc2 = bc(IB2,k,j,i);
        Real& bcc3 = bc(IB3,k,j,i);

        // cell center B-fields are defined as spatial interpolation at the volume center
        const Real& x1f_i  = pco->x1f(i);
        const Real& x1f_ip = pco->x1f(i+1);
        const Real& x1v_i  = pco->x1v(i);
        const Real& dx1_i  = pco->dx1f(i);
        Real lw=(x1f_ip-x1v_i)/dx1_i;
        Real rw=(x1v_i -x1f_i)/dx1_i;
        bcc1 = lw*b1_i + rw*b1_ip1;
        const Real& x2f_j  = pco->x2f(j);
        const Real& x2f_jp = pco->x2f(j+1);
        const Real& x2v_j  = pco->x2v(j);
        const Real& dx2_j  = pco->dx2f(j);
        lw=(x2f_jp-x2v_j)/dx2_j;
        rw=(x2v_j -x2f_j)/dx2_j;
        bcc2 = lw*b2_j + rw*b2_jp1;
        const Real& x3f_k  = pco->x3f(k);
        const Real& x3f_kp = pco->x3f(k+1);
        const Real& x3v_k  = pco->x3v(k);
        const Real& dx3_k  = pco->dx3f(k);
        lw=(x3f_kp-x3v_k)/dx3_k;
        rw=(x3v_k -x3f_k)/dx3_k;
        bcc3 = lw*b3_k + rw*b3_kp1;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief

void Field::CalculateCellCenteredFieldFourth(const FaceField &bf_center,
                                             AthenaArray<Real> &bc_center,
                                             Coordinates *pco, int il, int iu, int jl,
                                             int ju, int kl, int ku) {
  MeshBlock *pmb = pmy_block;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
		const Real& b1_im1 = bf_center.x1f(k  , j  , i-1);
        const Real& b1_i   = bf_center.x1f(k  , j  , i  );
        const Real& b1_ip1 = bf_center.x1f(k  , j  , i+1);
        const Real& b1_ip2 = bf_center.x1f(k  , j  , i+2);
        const Real& b2_j   = bf_center.x2f(k  , j  , i);
        const Real& b2_jp1 = bf_center.x2f(k  , j+1, i);
        const Real& b3_k   = bf_center.x3f(k  , j  , i);
		const Real& b3_kp1 = bf_center.x3f(k+1, j  , i);

        Real& bcc1 = bc_center(IB1,k,j,i);
        Real& bcc2 = bc_center(IB2,k,j,i);
        Real& bcc3 = bc_center(IB3,k,j,i);

        // cell center B-fields are defined as spatial interpolation at the volume center
        const Real& x1f_i  = pco->x1f(i);
        const Real& x1f_ip = pco->x1f(i+1);
        const Real& x1v_i  = pco->x1v(i);
        const Real& dx1_i  = pco->dx1f(i);
        Real lw=(x1f_ip-x1v_i)/dx1_i;
        Real rw=(x1v_i -x1f_i)/dx1_i;
		bcc1 = -1.0/16.0*(b1_im1 + b1_ip2) + 9.0/16.0*(b1_i + b1_ip1);

        const Real& x2f_j  = pco->x2f(j);
        const Real& x2f_jp = pco->x2f(j+1);
        const Real& x2v_j  = pco->x2v(j);
        const Real& dx2_j  = pco->dx2f(j);
        lw=(x2f_jp-x2v_j)/dx2_j;
        rw=(x2v_j -x2f_j)/dx2_j;
		if (pmb->block_size.nx2 > 1) {
		  const Real& b2_jm1 = bf_center.x2f(k,j-1,i);
		  const Real& b2_jp2 = bf_center.x2f(k,j+2,i);
		  bcc2 = -1.0/16.0*(b2_jm1 + b2_jp2) + 9.0/16.0*(b2_j + b2_jp1);
		} else { // default to second-order cell-centered field reconstruction in 1D:
		  bcc2 = 0.5*(b2_j + b2_jp1);
        }

        const Real& x3f_k  = pco->x3f(k);
        const Real& x3f_kp = pco->x3f(k+1);
        const Real& x3v_k  = pco->x3v(k);
        const Real& dx3_k  = pco->dx3f(k);
        lw=(x3f_kp-x3v_k)/dx3_k;
        rw=(x3v_k -x3f_k)/dx3_k;
		if (pmb->block_size.nx3 > 1) {
		  const Real& b3_km1 = bf_center.x3f(k-1,j,i);
		  const Real& b3_kp2 = bf_center.x3f(k+2,j,i);
		  bcc3 = -1.0/16.0*(b3_km1 + b3_kp2) + 9.0/16.0*(b3_k + b3_kp1);
		} else {
		  bcc3 = 0.5*(b3_k + b3_kp1);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief

void Field::CellCenteredToAveragedField(const AthenaArray<Real> &bc_center,
                                       AthenaArray<Real> &bc, Coordinates *pco,
                                       int il, int iu, int jl, int ju, int kl, int ku) {
  MeshBlock *pmb = pmy_block;
  // Transform cell-centered B to cell-averaged <B>
  // No need to add +1 ghost to longitudinal directions here as with FaceField
  AthenaArray<Real> laplacian_cc;
  laplacian_cc.InitWithShallowCopy(scr1_nkji_cc_);

  pco->Laplacian(bc_center, laplacian_cc, il, iu, jl, ju, kl, ku, 0, 2);

  // TODO(kfelker): deal with this usage of C
  Real h = pco->dx1f(il);  // pco->dx1f(i); inside loop
  Real C = (h*h)/24.0;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        Real& bc1 = bc(IB1,k,j,i);
        Real& bc2 = bc(IB2,k,j,i);
        Real& bc3 = bc(IB3,k,j,i);
		const Real& bcc1 = bc_center(IB1,k,j,i);
        const Real& bcc2 = bc_center(IB2,k,j,i);
        const Real& bcc3 = bc_center(IB3,k,j,i);

        bc1 = bcc1 + C*laplacian_cc(0,k,j,i);
        bc2 = bcc2 + C*laplacian_cc(1,k,j,i);
        bc3 = bcc3 + C*laplacian_cc(2,k,j,i);
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief
// Function now automatically shrinks the transverse loop limits for LaplacianX(),
// correction and output by 1. E.g. jl+1: ju-1 for Laplacian1() applied to bf.x1f
// As always, the function automatically extends the longitudinal loop limits to correct
// the upper FaceField. E.g. iu+1 for bf.x1f()

void Field::CalculateFaceCenteredField(const FaceField &bf, FaceField &bf_center,
                                       Coordinates *pco, int il, int iu, int jl, int ju,
                                       int kl, int ku) {
  // Convert face-averaged magnetic field to face-centered, even for 1D runs
  MeshBlock *pmb = pmy_block;
  BoundaryValues *pbval = pmb->pbval;

  // Laplacians (in orthogonal directions) of face-centered FaceField
  AthenaArray<Real> laplacian_bx1, laplacian_bx2, laplacian_bx3;
  laplacian_bx1.InitWithShallowCopy(scr1_kji_x1fc_);
  laplacian_bx2.InitWithShallowCopy(scr2_kji_x2fc_);
  laplacian_bx3.InitWithShallowCopy(scr3_kji_x3fc_);

  // Use 1x cell per boundary edge as buffer
  int il_buf=il, iu_buf=iu, jl_buf=jl, ju_buf=ju, kl_buf=kl, ku_buf=ku;
  int nl=0, nu=0;
  // If the x1 boundaries are periodic, use 1x cell on inner/outer boundary as buffer,
  // even in 1D. This is because all cells are passed to this function in a peridoic
  // domain, but only the real cells are passed for a physical boundary
  if (pbval->nblevel[1][1][0]!=-1) il_buf+=1;
  if (pbval->nblevel[1][1][2]!=-1) iu_buf-=1;

  if (pmb->block_size.nx2 > 1) {
	if (pmb->block_size.nx3 == 1) {// 2D
	  jl_buf+=1, ju_buf-=1;
    } else { // 3D
	  jl_buf+=1, ju_buf-=1, kl_buf+=1, ku_buf-=1;
    }
  }

  // Compute and store Laplacian of cell-averaged conserved variables
  pco->LaplacianX1(bf.x1f, laplacian_bx1, il, iu+1, jl_buf, ju_buf, kl_buf, ku_buf,
                   nl, nu);
  pco->LaplacianX2(bf.x2f, laplacian_bx2, il_buf, iu_buf, jl, ju+1, kl_buf, ku_buf,
                   nl, nu);
  pco->LaplacianX3(bf.x3f, laplacian_bx3, il_buf, iu_buf, jl_buf, ju_buf, kl, ku+1,
                   nl, nu);

  // TODO(kfelker): deal with this usage of C
  Real h = pco->dx1f(il);  // pco->dx1f(i); inside loop
  Real C = (h*h)/24.0;

  // Compute fourth-order approximation to cell-centered conserved variables
  for (int k=kl_buf; k<=ku_buf; ++k) {
    for (int j=jl_buf; j<=ju_buf; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        bf_center.x1f(k,j,i) =  bf.x1f(k,j,i) - C*laplacian_bx1(k,j,i);
      }
    }
  }
  for (int k=kl_buf; k<=ku_buf; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
      for (int i=il_buf; i<=iu_buf; ++i) {
        bf_center.x2f(k,j,i) =  bf.x2f(k,j,i) - C*laplacian_bx2(k,j,i);
      }
    }
  }

  for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl_buf; j<=ju_buf; ++j) {
      for (int i=il_buf; i<=iu_buf; ++i) {
        bf_center.x3f(k,j,i) =  bf.x3f(k,j,i) - C*laplacian_bx3(k,j,i);
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief

void Field::FaceAveragedToCellAveragedField(const FaceField &bf, FaceField &bf_center,
                                            AthenaArray<Real> &bc,
                                            AthenaArray<Real> &bc_center,
                                            Coordinates *pco, int il, int iu, int jl,
                                            int ju, int kl, int ku) {
  MeshBlock *pmb = pmy_block;
  BoundaryValues *pbval = pmb->pbval;
  // Assuming all cells (ghost and real) are passed as limits:
  CalculateFaceCenteredField(bf, bf_center, pco, il, iu, jl, ju, kl, ku);
  pbval->ApplyPhysicalBoundariesFaceField(bf_center);
  // ... output shrinks by 1 in transverse directions

  // Take 1x longitudinal cell as buffer, such that all directions shrink UPON OUTPUT.

  // it is important to emphasize that the local loop limits, il, iu, etc. represent the
  // validity of the output here.

  // the function receives ALL of the calculated values, but applies to il, iu only due to
  // stencil width
  if (pbval->nblevel[1][1][0]!=-1) il+=1;
  if (pbval->nblevel[1][1][2]!=-1) iu-=1;
  if (pbval->nblevel[1][0][1]!=-1) jl+=1;
  if (pbval->nblevel[1][2][1]!=-1) ju-=1;
  if (pbval->nblevel[0][1][1]!=-1) kl+=1;
  if (pbval->nblevel[2][1][1]!=-1) ku-=1;
  // Does this make any sense for outflow conditions? Need to revamp the whole il, iu, etc
  // buffer,  boundary condition syntax
  CalculateCellCenteredFieldFourth(bf_center, bc_center, pco, il, iu, jl, ju, kl, ku);
  pbval->ApplyPhysicalBoundariesCellField(bc_center);

  // All directions shrink by 1x again for Laplacian
  if (pbval->nblevel[1][1][0]!=-1) il+=1;
  if (pbval->nblevel[1][1][2]!=-1) iu-=1;
  if (pbval->nblevel[1][0][1]!=-1) jl+=1;
  if (pbval->nblevel[1][2][1]!=-1) ju-=1;
  if (pbval->nblevel[0][1][1]!=-1) kl+=1;
  if (pbval->nblevel[2][1][1]!=-1) ku-=1;

  CellCenteredToAveragedField(bc_center, bc, pco, il, iu, jl, ju, kl, ku);
  pbval->ApplyPhysicalBoundariesCellField(bc);

  // .... leading to output that is 2x cells smaller at each real edge

  return;
}
