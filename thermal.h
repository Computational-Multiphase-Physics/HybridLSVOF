#include "bcg.h"
#include "diffusion.h"

extern face vector uf;
extern double dt;
extern double rho1,rho2;
extern scalar f;

double lambda1 = 0,lambda2 = 0, cp1 = 1, cp2 = 1;
scalar T[];

#ifndef rhocp
# define rhocp(f) (clamp(f,0.,1.)*(rho1*cp1 - rho2*cp2) + rho2*cp2)
#endif
#ifndef lambdav
# define lambdav(f) (clamp(f,0.,1.)*(lambda1 - lambda2) + lambda2)
#endif


#if TREE
event defaults (i = 0)
{
    T.refine  = refine_linear;
    T.restriction = restriction_volume_average;
    //    T.gradient = p.gradient;
    T.dirty = true;
}
#endif // TREE

event tracer_advection(i++){
  advection ((scalar *){T}, uf, dt);
}

face vector lambdav[];
scalar rhocp[];
event tracer_diffusion(i++){
  foreach_face(){
    double ff = (f[] + f[-1])/2.;
    lambdav.x[] = fm.x[]*lambdav(ff);
  }
  foreach()
    rhocp[] = cm[]*(rhocp(f[]));

  diffusion(T,dt,lambdav,theta=rhocp);
}
