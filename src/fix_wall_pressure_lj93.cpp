/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include "fix_wall_pressure_lj93.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallPressureLJ93::FixWallPressureLJ93
  (LAMMPS *lmp, int narg, char **arg) :
  FixWallPressure(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixWallPressureLJ93::precompute()
{
  coeff1 = 6.0/5.0 * epsilon * pow(sigma,9.0);
  coeff2 = 3.0 * epsilon * pow(sigma,3.0);
  coeff3 = 2.0/15.0 * epsilon * pow(sigma,9.0);
  coeff4 = epsilon * pow(sigma,3.0);

  double rinv = 1.0/cutoff;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  offset = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv;
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallPressureLJ93::wall_particle(int which, double coord)
{
  double delta,rinv,r2inv,r4inv,r10inv,frc;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r10inv = r4inv*r4inv*r2inv;
      frc = side * (coeff1*r10inv - coeff2*r4inv);
      f[i][dim] -= frc;
      ewall += coeff3*r4inv*r4inv*rinv -
        coeff4*r2inv*rinv - offset;
      fwall += frc;
    }

  if (onflag) error->one(FLERR,"Particle on or inside fix wall surface");
}
