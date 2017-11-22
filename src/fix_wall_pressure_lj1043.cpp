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

/* ----------------------------------------------------------------------
   Contributing author: Jonathan Lee (Sandia)
------------------------------------------------------------------------- */

#include <math.h>
#include "fix_wall_pressure_lj1043.h"
#include "atom.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixWallPressureLJ1043::FixWallPressureLJ1043(LAMMPS *lmp, int narg, char **arg) :
  FixWallPressure(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixWallPressureLJ1043::precompute()
{
  coeff1 = MY_2PI * 2.0/5.0 * epsilon * pow(sigma,10.0);
  coeff2 = MY_2PI * epsilon * pow(sigma,4.0);
  coeff3 = MY_2PI * pow(2.0,1/2.0) / 3 * epsilon * pow(sigma,3.0);
  coeff4 = 0.61 / pow(2.0,1/2.0) * sigma;
  coeff5 = coeff1 * 10.0;
  coeff6 = coeff2 * 4.0;
  coeff7 = coeff3 * 3.0;

  double rinv = 1.0/cutoff;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  offset = coeff1*r4inv*r4inv*r2inv - coeff2*r4inv -
	coeff3*pow(cutoff+coeff4,-3.0);
}

/* ---------------------------------------------------------------------- */

void FixWallPressureLJ1043::wall_particle(int which, double coord)
{
  double delta,rinv,r2inv,r4inv,r10inv,frc;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta <= 0.0) continue;
      if (delta > cutoff) continue;
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r10inv = r4inv*r4inv*r2inv;

      frc = side * (coeff5*r10inv*rinv - coeff6*r4inv*rinv -
	coeff7*pow(delta+coeff4,-4.0));
      f[i][dim] -= frc;
      ewall += coeff1*r10inv - coeff2*r4inv -
	coeff3*pow(delta+coeff4,-3.0) - offset;
      fwall += frc;
    }
}
