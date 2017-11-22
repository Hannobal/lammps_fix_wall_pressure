/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/pressure/harmonic,FixWallPressureHarmonic)

#else

#ifndef LMP_FIX_WALL_PRESSURE_HARMONIC_H
#define LMP_FIX_WALL_PRESSURE_HARMONIC_H

#include "fix_wall_pressure.h"

namespace LAMMPS_NS {

class FixWallPressureHarmonic : public FixWallPressure {
 public:
  FixWallPressureHarmonic(class LAMMPS *, int, char **);
  void precompute() {}
  void wall_particle(int, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Particle on or inside fix wall surface

Particles must be "exterior" to the wall in order for energy/force to
be calculated.

*/
