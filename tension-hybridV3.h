//This has exact lines from SPs integral.h (unlike LS it did not improve the things)
//For begining we calculate the sign distance fn at cell centers by simple averaging
#include "voftols.h"

#undef SEPS
#define SEPS 1e-30

scalar d[];
extern scalar kappa;
extern scalar sigma;

void vof2distCC(scalar c){

  vertex scalar ls[];
  vof2dist(c,ls);

  foreach(){
    int num = 0;
    d[] = 0.;
    for(int ii = 0; ii <= 1; ii++)
      for(int jj = 0; jj <= 1; jj++){
	if(fabs(ls[ii,jj]) < L0){
	  d[] += ls[ii,jj];
	  num += 1;
	}
      }
    if(num > 0)
      d[] /= num;
    else
      d[] = nodata;
  }
}

event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
  }
}

event stability (i++)
{
  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;
    }
  double rhom = (1./amin + 1./amax)/2.;

  foreach(reduction(min:dtmax)){
    double dt = sqrt (rhom*cube(dmin)/(pi*sigma[]));
    if (dt < dtmax)
      dtmax = dt;
  }
}

static inline double curvatureLS (Point point, scalar d) {

  double dx = (d[1] - d[-1])/2.;
  double dy = (d[0,1] - d[0,-1])/2.;
  double dxx = d[1] - 2.*d[] + d[-1];
  double dyy = d[0,1] - 2.*d[] + d[0,-1];
  double dxy = (d[1,1] - d[-1,1] - d[1,-1] + d[-1,-1])/4.;
  double dn = sqrt(sq(dx) + sq(dy)) + 1e-30;

  bool flagnd = true;
  foreach_dimension()
    if(fabs(d[]) < nodata || fabs(d[1]) < nodata || fabs(d[-1]) < nodata)
      flagnd = false;
  
  if(flagnd)
    return nodata;
  else
    return (sq(dx)*dyy - 2.*dx*dy*dxy + sq(dy)*dxx)/cube(dn)/Delta;
}

event properties(i++){
  vof2distCC(f);
  curvature(f, kappa);
  boundary ({kappa,d});

  foreach(){
    if(kappa[] == nodata && fabs(d[]) != nodata)
      kappa[] = curvatureLS (point, d); // Not true, why did i say earlier that "I am opposite to curvature fn of basilisk"
  }
  boundary ((scalar *){kappa});
}

event acceleration (i++)
{
  tensor S[];

  foreach (){
    foreach_dimension(){
      S.y.y[] = 0.0;
      for (int i = -1; i <= 1; i += 2){
	if(fabs(d[]) < nodata && fabs(d[1]) < nodata && fabs(d[-1]) < nodata && fabs(kappa[]) < nodata){ // FIX me
  	  if (d[]*(d[] + d[i]) < 0.) {
	    double xi = d[]/(d[] - d[i] + SEPS);
	    double nx = ((d[1] - d[-1])/2. +
			 xi*i*(d[-1] - 2.*d[] + d[1]))/Delta;
  	    double ki = kappa[];// + xi*(kappa[i] - kappa[]);
	    double sigmaxi = sigma[] + xi * (sigma[i] - sigma[]);
	    S.y.y[] += sigmaxi*(fabs(nx)/Delta - sign(d[])*ki*(0.5 - xi));
  	  }
  	}
      }
    }
  }
  
  boundary({S.x.x,S.y.y});

  foreach_vertex(){
    foreach_dimension(){
      S.x.y[] = 0.0;
      if ((d[] + d[0, -1]) * (d[-1] + d[-1, -1]) > 0.) // perhaps i can use vertex values
	S.x.y[] = 0.0;
      else{
	if(fabs(d[]) < nodata && fabs(d[-1]) < nodata && fabs(d[-1,-1]) < nodata &&fabs(d[0,-1]) < nodata){ //FIX me
	  double xi = (d[-1] + d[-1,-1])/(d[-1] + d[-1,-1] - d[] - d[0,-1]);
	  double ny = (xi*(d[] - d[-1] + d[-1,-1] - d[0,-1]) +
		       d[-1] - d[-1,-1])/Delta;
	  double sigmaxi = 0.5 * (sigma[-1] + sigma[-1,-1] + xi * (sigma[] - sigma[-1] - sigma [-1,-1] + sigma[0,-1]));
	  S.x.y[] = - sigmaxi*sign(d[] + d[0,-1])*ny/Delta;
	}
      }
    }
  }
  boundary({S.x.y,S.y.x});
  
  face vector ia = a;
  foreach_face(){
    if(fabs(kappa[]) < nodata){
      ia.x[] += alpha.x[] / fm.x[] * (S.x.x[] - S.x.x[-1] + S.x.y[0,1] - S.x.y[]) / Delta;
    }
  }

  boundary((scalar *){ia});
}


