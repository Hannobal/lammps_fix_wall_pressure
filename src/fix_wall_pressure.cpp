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

/* format:
fix <ID> <group> wall/pressure/<style> <direction (e.g. zhi)> <position>
<epsilon> <sigma> <cutoff> <presure> <density> <damping>
!!! density is in (unit of mass)/(unit of length)^2, not in density unit!!!
*/

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "stdio.h"
#include "fix_wall_pressure.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */


FixWallPressure::FixWallPressure(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  time_integrate = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  restart_global = 1;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // parse args

  int scaleflag = 1;
  fldflag = 0;
  int pbcflag = 0;
  
  for (int i = 0; i < 6; i++)
    estr = sstr = pstr = mstr = dstr = NULL;

  int iarg = 3;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+8 > narg) error->all(FLERR,"Illegal fix wall command");
      
      if      (strcmp(arg[iarg],"xlo") == 0) wallwhich = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) wallwhich = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) wallwhich = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) wallwhich = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) wallwhich = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) wallwhich = ZHI;
      
      xwall = force->numeric(FLERR,arg[iarg+1]);
      
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        estr = new char[n];
        strcpy(estr,&arg[iarg+2][2]);
        estyle = VARIABLE;
      } else {
        epsilon = force->numeric(FLERR,arg[iarg+2]);
        estyle = CONSTANT;
      }
      
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        sstr = new char[n];
        strcpy(sstr,&arg[iarg+3][2]);
        sstyle = VARIABLE;
      } else {
        sigma = force->numeric(FLERR,arg[iarg+3]);
        sstyle = CONSTANT;
      }
      
      cutoff = force->numeric(FLERR,arg[iarg+4]);
      
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
        int n = strlen(&arg[iarg+5][2]) + 1;
        pstr = new char[n];
        strcpy(pstr,&arg[iarg+5][2]);
        pstyle = VARIABLE;
      } else {
        pressure = force->numeric(FLERR,arg[iarg+5]);
        pstyle = CONSTANT;
      }
      
      if (strstr(arg[iarg+6],"v_") == arg[iarg+6]) {
        int n = strlen(&arg[iarg+6][2]) + 1;
        mstr = new char[n];
        strcpy(mstr,&arg[iarg+6][2]);
        mstyle = VARIABLE;
      } else {
        density = force->numeric(FLERR,arg[iarg+6]);
        mstyle = CONSTANT;
      }
      
      if (strstr(arg[iarg+7],"v_") == arg[iarg+7]) {
        int n = strlen(&arg[iarg+7][2]) + 1;
        dstr = new char[n];
        strcpy(dstr,&arg[iarg+7][2]);
        dstyle = VARIABLE;
      } else {
        damping = force->numeric(FLERR,arg[iarg+7]);
        dstyle = CONSTANT;
      }

      /*if(me==0) {
        std::cout << "#########READ:" << std::endl
        << "face " << wallwhich
        << " pos " << xwall
        << " eps " << epsilon
        << " sig " << sigma
        << " cut " << cutoff
        << " pres " << pressure
        << " rho " << density
        << " damp " << damping << std::endl;
      }*/
      
      iarg += 8;
      
    } else if (strcmp(arg[iarg],"fld") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"no") == 0) fldflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) fldflag = 1;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"yes") == 0) pbcflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pbcflag = 0;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall command");
  }
  
  size_vector = 4;

  // error checks

  if (cutoff <= 0.0) error->all(FLERR,"Fix wall cutoff <= 0");
  if (epsilon < 0.0) error->all(FLERR,"Fix wall epsilon < 0");
  if (sigma < 0.0) error->all(FLERR,"Fix wall sigma < 0.0");
  if (density <= 0.0) error->all(FLERR,"Fix wall density <= 0");
  if (pressure <= 0.0) error->all(FLERR,"Fix wall pressure <= 0");
  if (damping <= 0.0 || damping > 1.0)
    error->all(FLERR,"Fix wall damping <= 0 or damping >= 1");

  if ((wallwhich == ZLO || wallwhich == ZHI) && domain->dimension == 2)
    error->all(FLERR,"Cannot use fix wall zlo/zhi for a 2d simulation");

  if (!pbcflag) {
    if ((wallwhich == XLO || wallwhich == XHI) && domain->xperiodic)
      error->all(FLERR,"Cannot use fix wall in periodic dimension");
    if ((wallwhich == YLO || wallwhich == YHI) && domain->yperiodic)
      error->all(FLERR,"Cannot use fix wall in periodic dimension");
    if ((wallwhich == ZLO || wallwhich == ZHI) && domain->zperiodic)
      error->all(FLERR,"Cannot use fix wall in periodic dimension");
  }

  // set varflag if any wall parameters are variable
  // set wstyle to VARIABLE if any wall parameters are variable

  varflag = 0;
  if (estyle == VARIABLE || sstyle == VARIABLE || pstyle == VARIABLE
    || mstyle == VARIABLE || dstyle == VARIABLE) {
    varflag = 1;
    wstyle = VARIABLE;
  } else {
    wstyle = CONSTANT;
  }

  ewall = 0.0;
  fwall = 0.0;
}


/* ---------------------------------------------------------------------- */

FixWallPressure::~FixWallPressure()
{
  delete estr;
  delete sstr;
  delete pstr;
  delete mstr;
  delete dstr;
}

/* ---------------------------------------------------------------------- */

int FixWallPressure::setmask()
{
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallPressure::init()
{
  if (estyle == VARIABLE) {
    eindex = input->variable->find(estr);
    if (eindex < 0)
      error->all(FLERR,"Variable name for fix wall does not exist");
    if (!input->variable->equalstyle(eindex))
      error->all(FLERR,"Variable for fix wall is invalid style");
  }
  if (sstyle == VARIABLE) {
    sindex = input->variable->find(sstr);
    if (sindex < 0)
      error->all(FLERR,"Variable name for fix wall does not exist");
    if (!input->variable->equalstyle(sindex))
      error->all(FLERR,"Variable for fix wall is invalid style");
  }
  if (pstyle == VARIABLE) {
    pindex = input->variable->find(pstr);
    if (pindex < 0)
      error->all(FLERR,"Variable name for fix wall does not exist");
    if (!input->variable->equalstyle(pindex))
      error->all(FLERR,"Variable for fix wall is invalid style");
  }
  if (mstyle == VARIABLE) {
    mindex = input->variable->find(mstr);
    if (mindex < 0)
      error->all(FLERR,"Variable name for fix wall does not exist");
    if (!input->variable->equalstyle(mindex))
      error->all(FLERR,"Variable for fix wall is invalid style");
  }
  if (dstyle == VARIABLE) {
    dindex = input->variable->find(dstr);
    if (dindex < 0)
      error->all(FLERR,"Variable name for fix wall does not exist");
    if (!input->variable->equalstyle(dindex))
      error->all(FLERR,"Variable for fix wall is invalid style");
  }

  precompute();

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
  
  dtv    = update->dt;
  nktv2p = force->nktv2p;
  dtf    = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixWallPressure::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    if (!fldflag) post_force(vflag);
  } else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallPressure::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   only called if fldflag set, in place of post_force
------------------------------------------------------------------------- */

void FixWallPressure::pre_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallPressure::post_force(int vflag)
{

  // xwall = current position of wall
  // evaluate variables if necessary, wrap with clear/add
  // for parameters variables need to re-invoke precompute()

  if (varflag) modify->clearstep_compute();

  ewall = 0.0;
  fwall = 0.0;
  area = calc_area();
  if (wstyle == VARIABLE) {
    if (estyle == VARIABLE) {
      epsilon = input->variable->compute_equal(eindex);
      if (epsilon < 0.0)
        error->all(FLERR,"Variable evaluation in fix wall gave bad value");
    }
    if (sstyle == VARIABLE) {
      sigma = input->variable->compute_equal(sindex);
      if (sigma < 0.0)
        error->all(FLERR,"Variable evaluation in fix wall gave bad value");
    }
    if (pstyle == VARIABLE) {
      pressure = input->variable->compute_equal(pindex);
      if (pressure < 0.0)
        error->all(FLERR,"Variable evaluation in fix wall gave bad value");
    }
    if (mstyle == VARIABLE) {
      mass = input->variable->compute_equal(mindex);
      if (mass < 0.0)
        error->all(FLERR,"Variable evaluation in fix wall gave bad value");
    } else {
      mass = area*density;
    }
    if (dstyle == VARIABLE) {
      damping = input->variable->compute_equal(dindex);
      if (damping < 0.0)
        error->all(FLERR,"Variable evaluation in fix wall gave bad value");
    }
    precompute();
  } else {
    mass = area*density;
  }

  wall_particle(wallwhich,xwall);
  buf1[0]=fwall;
  buf1[1]=ewall;
  MPI_Allreduce(buf1,buf2,2,MPI_DOUBLE,MPI_SUM,world);
  fwall=buf2[0];
  ewall=buf2[1];
  if(wallwhich%2==0) // XLO YLO ZLO
    fwall += pressure*area/nktv2p;
  else // XHI YHI ZHI
    fwall -= pressure*area/nktv2p;

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixWallPressure::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallPressure::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixWallPressure::calc_area() {
  // area is the cross section of the cell perpendicular to the axis
  if(wallwhich<2) // x-xaxis
    return (domain->boxhi[1] - domain->boxlo[1]) *
    (domain->boxhi[2] - domain->boxlo[2]);
  else if(wallwhich<4) // y-axis
    return (domain->boxhi[0] - domain->boxlo[0]) *
    (domain->boxhi[2] - domain->boxlo[2]);
  else // z-axis
    return (domain->boxhi[0] - domain->boxlo[0]) *
    (domain->boxhi[1] - domain->boxlo[1]);
}

void FixWallPressure::initial_integrate(int vflag)
{
  double dtfm;
  dtfm = dtf/mass;
  vwall += dtfm * fwall;
  xwall += dtv * vwall;
  /*if(me==0 && update->ntimestep%100==0) {
    std::cout << update->ntimestep
    << " x " << xwall
    << " v " << vwall
    << " f " << fwall
    << " e " << ewall
    << " m " << mass
    << " a " << area
    << " dtv " << dtv
    << std::endl;
  }*/
}

/* ---------------------------------------------------------------------- */

void FixWallPressure::final_integrate()
{
  double area,dtfm;
    dtfm = dtf/mass;
    vwall += dtfm * fwall;
    vwall *= damping;
}

void FixWallPressure::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

void FixWallPressure::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

void FixWallPressure::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

void FixWallPressure::write_restart(FILE *fp)
{
  if (me) return;
  //std::cout << "########### writing pressure wall restart"<<std::endl;
  buf1[0]=xwall;
  buf1[1]=vwall;
  int size = 2*sizeof(double);
  fwrite(&size,sizeof(int),1,fp);
  fwrite(&buf1[0],sizeof(double),2,fp);
}

void FixWallPressure::restart(char *buf)
{
  /*if (me == 0) {
    std::cout << "########### restarting pressure wall"<<std::endl;
  }*/
  double *list = (double *)buf;
  xwall=list[0];
  vwall=list[1];
}
double FixWallPressure::compute_scalar()
{
  return ewall;
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWallPressure::compute_vector(int n)
{
  if(n==0) return ewall;
  else if(n==1) return xwall;
  else if(n==2) return vwall;
  else if(n==3) return fwall;
  else return 0.0;
}
