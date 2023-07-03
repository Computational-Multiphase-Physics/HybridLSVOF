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

  //fixme: try high-order discretization for the computation of curvature
  double dx = (d[1] - d[-1])/2.;
  double dy = (d[0,1] - d[0,-1])/2.;
  double dxx = d[1] - 2.*d[] + d[-1];
  double dyy = d[0,1] - 2.*d[] + d[0,-1];
  double dxy = (d[1,1] - d[-1,1] - d[1,-1] + d[-1,-1])/4.;
  double dn = sqrt(sq(dx) + sq(dy)) + SEPS;
  double curv = (sq(dx)*dyy - 2.*dx*dy*dxy + sq(dy)*dxx)/cube(dn)/Delta;

#if AXI
  double r = y;
  curv += dy*(dy*dy + dx*dx)/r / cube(dn);
#endif

#if (AXI || dimension == 3) // !AXI
  curv = 2.0 * curv / (2.0 - d[] * curv);
#elif (dimension == 2)
  curv = 1.0 * curv / (1.0 - d[] * curv);
#endif // !AXI

  //for robustnes
  curv = min(max(curv, -2.0 / Delta), 2.0 / Delta);

  bool flagnd = true;
  foreach_dimension()
      if(fabs(d[]) < nodata || fabs(d[1]) < nodata || fabs(d[-1]) < nodata)
	flagnd = false;

  if(flagnd)
    return nodata;
  else
    return -curv;
}

event properties(i++){
  vof2distCC(f);
  curvature(f, kappa);
  boundary ({kappa,d});

  foreach(){
    if(kappa[] >= nodata)
      kappa[] = curvatureLS (point, d); // Not true, why did i say earlier that "I am opposite to curvature fn of basilisk"
  }
  boundary ((scalar *){kappa});
}

event acceleration (i++)
{
  face vector ia = a;
  vector off_diag[];
  vector diag[];
  off_diag.n[left] = neumann(0);
  off_diag.n[right] = neumann(0);
  off_diag.n[top] = neumann(0);
  off_diag.n[bottom] = neumann(0);
  diag.n[left] = neumann(0);
  diag.n[right] = neumann(0);
  diag.n[top] = neumann(0);
  diag.n[bottom] = neumann(0);

  vertex scalar t_x[];
  vertex scalar t_y[];

    foreach ()
    {
      foreach_dimension()
      {
        diag.x[] = 0.0;
        for (int k = -1; k <= 1; k += 2)
        {
          double phi_c = d[];
          double phi_ni = d[0, k];
	  if(fabs(d[]) < nodata && fabs(d[0,1]) < nodata && fabs(d[0,-1]) < nodata && fabs(kappa[]) < nodata){ // FIX me
	    if (phi_c * (phi_c + phi_ni) <= 0)
	      {
		double xi = phi_c / (phi_c - phi_ni + SEPS);
		double tan = 1.0 / Delta * ((d[0, 1] - d[0, -1]) / 2.0
					    + k * xi * (d[0, -1] - 2.0 * d[] + d[0, 1]));
		double ka = kappa[];
		double sigmaxi = sigma[] + xi * (sigma[k] - sigma[]);
		diag.x[] += sigmaxi * (fabs(tan) / Delta - sign(phi_c) * ka * (0.5 - xi));

	      }
	  }
	}
      }
    }

    boundary((scalar *){diag});

     foreach_vertex()
     {
      foreach_dimension()
      {
        off_diag.y[] = 0.0;
        if ((d[] + d[0, -1]) * (d[-1] + d[-1, -1]) > 0.) // perhaps i can use vertex values
          off_diag.y[] = 0.0;
        else
        {
	  if(fabs(d[]) < nodata && fabs(d[-1]) < nodata && fabs(d[-1,-1]) < nodata &&fabs(d[0,-1]) < nodata){ //FIX me
	    double xi = (d[-1] + d[-1, -1]) / (d[-1] + d[-1, -1] - d[] - d[0, -1]);
	    double ny = (xi * (d[] - d[-1] + d[-1, -1] - d[0, -1]) +
			 d[-1] - d[-1, -1]) / Delta;
	    double sigmaxi = 0.5 * (sigma[-1] + sigma[-1,-1] + xi * (sigma[] - sigma[-1] - sigma [-1,-1] + sigma[0,-1]));
	    off_diag.y[] = - sigmaxi * sign(d[] + d[0, -1]) * ny / Delta;
	  }
	}
      }
     }

    foreach_face()
    {
      if(fabs(kappa[]) < nodata)
	ia.x[] += alpha.x[] / fm.x[] * (diag.x[] - diag.x[-1] + off_diag.y[0, 1] - off_diag.y[]) / Delta;
    }
  boundary((scalar *){ia});
}


