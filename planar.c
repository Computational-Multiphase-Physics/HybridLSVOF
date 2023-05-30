#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "curvature.h"
#include "thermal.h"
#include "tension-hybrid.h"
//#include "view.h"
//#include "vtkdaniel.h"

int MAXLEVEL = 8;
double Re = 0.01, Ma = 0.01, Ca = 0.01;
double deltaT = 1.;
double sigma0,sigmaT;
double tend = 0.1;
scalar sigma[],kappa[];
vector SurfGradSigma[];

T[top] = dirichlet(1.);
T[bottom] = dirichlet((1.+deltaT*cos(2.*pi/L0*x)));

uf.n[top] = dirichlet(0.);
uf.t[top] = dirichlet(0.);

uf.n[bottom] = dirichlet(0.);
uf.t[bottom] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.n[top] = dirichlet(0.);

u.t[bottom] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);

int main(){
  L0=8.;
  Y0 = -L0/8.;
  X0 = -L0/2.;
  mu1 = 1./Re;
  mu2 = 1./Re;
  sigma0 = 1./Re/Ca;
  sigmaT = L0/Re/deltaT;
  lambda1 = cp1/Ma;
  lambda2 = cp2/Ma;

  periodic(right);
  TOLERANCE = 1e-6;
  for(MAXLEVEL = 7; MAXLEVEL <=9; MAXLEVEL++){
  init_grid(1<<MAXLEVEL);
  run();
  }
}

event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

event init (i++){
  mask (y > 1. ? top : none);

  fraction(f,(y-((double)L0)/((double)N)/2.));

  foreach(){
    sigma[] = sigma0 + sigmaT*(T[] - 1.);
    //    T[] = 1.;
  }
  boundary((scalar *){sigma});
}


event properties (i++){
  foreach()
    sigma[] = sigma0 + sigmaT*(T[] - 1.);
  boundary((scalar *){sigma});
}


/* event movie (t += tend/90){ */

/*   char s[80]; */
/*   sprintf (s, "t = %.3f", t); */
/*   clear(); */
/*   view(fov = 20., ty = 0, quat = {0,0,0,0}, width = 1980, height = 1980); */
/*   draw_vof("f",lw=4); */
/*   squares("T", min = 1., max = 1.+deltaT, linear=true, map = cool_warm); */
/*   vectors("u",0.075,lc = {0,0,0}, lw = 4); */
  
/*   draw_string ("HOT", pos = 4, size = 30, lc = {0,0,0}, lw = 4); */
/*   draw_string ("COLD", pos = 2, size = 30, lc = {0,0,0}, lw = 4); */

/*   draw_string (s, pos = 1, size = 30, lc = {255,255,255}, lw = 4);     */
/*   char filename[60]; */
/*   sprintf(filename,"planar-%d.mp4",MAXLEVEL); */
/*   save(filename); */
/* } */

double prevu=0.;
event logs(i++){
  double uavg = normf(u.x).avg;
  fprintf(ferr,"%g %g %g %g\n",t,statsf(T).max,uavg,fabs(uavg - prevu));
  prevu = uavg;
}

vector acc[];
event end(t = tend){

  foreach(){
    foreach_dimension()
      acc.x[] = (a.x[] + a.x[1])/2.;
  }
  boundary((scalar *){acc});
  
  /* char filename[60]; */
  /* sprintf(filename,"output-%d.vtk",MAXLEVEL); */
  /* FILE * fp = fopen(filename,"w"); */
  /* output_vtk((scalar *){T,u,p,f,acc},fp,{{-10,-10},{10,10}}); */
  /* fclose(fp); */

  char fileout[60];
  sprintf(fileout,"out-%d.dat",MAXLEVEL);
  FILE *fpo = fopen(fileout,"w");
  
  for(double y = -1.; y < 1.; y+= 0.01){
    fprintf(fpo,"%g %g %g\n",y,interpolate(u.y,0.,y),interpolate(T,0.,y));
  }
  fclose(fpo);

  char namevec[80];
  sprintf (namevec, "vector-%d.dat",MAXLEVEL);
  FILE * fpv = fopen (namevec, "w");
  foreach()
    fprintf(fpv,"%g %g %g %g\n",x,y,u.x[],u.y[]);
  fclose (fpv);

}
