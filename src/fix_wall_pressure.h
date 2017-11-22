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

#ifndef LMP_FIX_WALL_PRESSURE_H
#define LMP_FIX_WALL_PRESSURE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallPressure : public Fix {
 public:
  int wallwhich;
  double xwall;

  FixWallPressure(class LAMMPS *, int, char **);
  virtual ~FixWallPressure();
  int setmask();
  virtual void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  void initial_integrate(int);
  void final_integrate();
  void reset_dt();
  void write_restart(FILE *);
  void restart(char *);
  double calc_area();
  double compute_scalar();
  double compute_vector(int);

  virtual void precompute() = 0;
  virtual void wall_particle(int, double) = 0;

 protected:
  double epsilon,sigma,cutoff,pressure,density,mass,damping;
  double vwall,fwall,ewall,area;
  double dtv,dtf,nktv2p;
  double buf1[2],buf2[2];
  int me,nprocs;
  int estyle,sstyle,pstyle,mstyle,dstyle,wstyle;
  int eindex,sindex,pindex,mindex,dindex;
  char *estr,*sstr,*pstr,*mstr,*dstr;
  int varflag; // 1 if any wall position,epsilon,sigma is a var
  int ilevel_respa;
  int fldflag;
  double *step_respa;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Wall defined twice in fix wall command

Self-explanatory.

E: Fix wall cutoff <= 0.0

Self-explanatory.

E: Cannot use fix wall zlo/zhi for a 2d simulation

Self-explanatory.

E: Cannot use fix wall in periodic dimension

Self-explanatory.

E: Variable name for fix wall does not exist

Self-explanatory.

E: Variable for fix wall is invalid style

Only equal-style variables can be used.

E: Variable evaluation in fix wall gave bad value

The returned value for epsilon or sigma < 0.0.

*/