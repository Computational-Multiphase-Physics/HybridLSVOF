#include "run.h"
#include "timestep.h"
#include "level-set.h"
#include "fractions.h"

(const) face vector a = zerof, alpha = unityf;
double dtmax = 0.1;

#include "tension-hybridV2.h"


void levelset2vof(scalar d, scalar f)
{
  vertex scalar phi_cor[];
  face vector s_tmp[];
  foreach_vertex()
  {
    phi_cor[] = 0.25 * (d[] + d[-1] + d[-1,-1] + d[0,-1]);
  }
  fractions(phi_cor, f, s_tmp);
}

vector u[];
scalar kappa[],sigma[];
scalar f[];

int main(){
  L0=8.;
  init_grid(64);
  run();
}

event init (i++){
  dt = 0.1;
  foreach(){
    u.x[] = 5.*cos(atan2(y,(x-0.5*L0)))/sqrt(sq(x-0.5*L0) + y*y);
    u.y[] = 5.*sin(atan2(y,(x-0.5*L0)))/sqrt(sq(x-0.5*L0) + y*y);	
    d[] = sqrt(sq(x-0.5*L0) + sq(y)) - 1;
    kappa[] = curvatureLS (point, d); // I am opposite to curvature fn of basilisk
  }
  boundary((scalar *){u,kappa});

}

event advection(t += 0.1){
  advectLevelSet(d, u, dt, i);
}

event logs(i++){
   fprintf(ferr,"%g %g\n",t,dt);
}

event stability (i++) {
  dt = 0.1;
} // event stability ends

event end(t = 1000.){
  levelset2vof(d, f);
  dump("last");
    }
