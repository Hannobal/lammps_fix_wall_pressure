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
#include "fix_wall_pressure_lj126.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallPressureLJ126::FixWallPressureLJ126
  (LAMMPS *lmp, int narg, char **arg) :
  FixWallPressure(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixWallPressureLJ126::precompute()
{
  coeff1 = 48.0 * epsilon * pow(sigma,12.0);
  coeff2 = 24.0 * epsilon * pow(sigma,6.0);
  coeff3 = 4.0 * epsilon * pow(sigma,12.0);
  coeff4 = 4.0 * epsilon * pow(sigma,6.0);

  double r2inv = 1.0/(cutoff*cutoff);
  double r6inv = r2inv*r2inv*r2inv;
  offset = r6inv*(coeff3*r6inv - coeff4);
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallPressureLJ126::wall_particle(int which, double coord)
{
  double delta,rinv,r2inv,r6inv,frc;

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
      r6inv = r2inv*r2inv*r2inv;
      frc = side * r6inv*(coeff1*r6inv - coeff2) * rinv;
      f[i][dim] -= frc;
      ewall     += r6inv*(coeff3*r6inv - coeff4) - offset;
      fwall     += frc;
    }

  if (onflag) error->one(FLERR,"Particle on or inside fix wall surface");
}
