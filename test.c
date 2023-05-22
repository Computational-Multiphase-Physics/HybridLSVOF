/* #include "utils.h" */
/* #include "navier-stokes/centered.h" */
#include "tension-hybrid.h"

scalar f[];
#define Drop(x,y) (sq(x - 1.) + sq(y))
scalar kappa[],sigma[];
vector off_diag[];
vector diag[];
vector diagact[];

scalar dact[];
face vector ia[];

int main(){
  init_grid (64);
  L0 = 2;
  fraction(f,sq(0.75)-(Drop(x,y)));

  vof2distCC(f);

  foreach(){
    dact[] = sq(0.75)-(Drop(x,y));
  }
  
  /* foreach(){ */
  /*   kappa[] = curvatureLS (point, d); // I am opposite to curvature fn of basilisk */
  /*   sigma[] = 1. + x; */
  /* } */
  /* boundary ((scalar *){kappa,sigma}); */
  
  foreach ()
    {
      foreach_dimension()
  	{
  	  diag.x[] = 0.0;
	  diagact.x[] = 0.0;
	  
  	  for (int k = -1; k <= 1; k += 2)
  	    {
  	      double phi_c = d[];
  	      double phi_ni = d[0, k];
  	      if(phi_c < nodata && phi_ni < nodata){
  		if (phi_c * (phi_c + phi_ni) <= 0)
  		  {
  		    diag.x[] = 1.;
  		  }
  	      }
  	      else
  	      	diag.x[] = 0.;

	      double phi_c2 = d[];
  	      double phi_ni2 = d[0, k];
  	      if(phi_c2 < nodata && phi_ni2 < nodata){
  		if (phi_c2 * (phi_c2 + phi_ni2) <= 0)
  		  {
  		    diagact.x[] = 1.;
  		  }
  	      }
  	      else
  	      	diagact.x[] = 0.;
  	    }
  	}
    }
  
  /* boundary((scalar *){diag}); */
  
  /* foreach_vertex() */
  /*   { */
  /*     foreach_dimension() */
  /* 	{ */
  /* 	  off_diag.y[] = 0.0; */
  /* 	  if ((d[] + d[0, -1]) * (d[-1] + d[-1, -1]) > 0.) */
  /* 	    off_diag.y[] = 0.0; */
  /* 	  else */
  /* 	    { */
  /* 	      double xi = (d[-1] + d[-1, -1]) / (d[-1] + d[-1, -1] - d[] - d[0, -1]); */
  /* 	      double ny = (xi * (d[] - d[-1] + d[-1, -1] - d[0, -1]) + */
  /* 			   d[-1] - d[-1, -1]) / Delta; */
  /* 	      double sigmaxi = 0.5 * (sigma[-1] + sigma[-1,-1] + xi * (sigma[] - sigma[-1] - sigma [-1,-1] + sigma[0,-1])); */
  /* 	      off_diag.y[] = - sigmaxi * sign(d[] + d[0, -1]) * ny / Delta; */
  /* 	    } */
  /* 	} */
  /*   } */
  
  /* foreach_face() */
  /*   { */
  /*     ia.x[] += (diag.x[] - diag.x[-1] + off_diag.y[0, 1] - off_diag.y[]) / Delta; */
  /*   } */
  /* boundary((scalar *){ia}); */
  
  dump("init");
}
