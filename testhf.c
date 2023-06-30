#include "utils.h"
#include "navier-stokes/centered.h"

scalar f[];

#include "fractions.h"
#include "curvature.h"
#include "tension-hybridV4.h"

scalar kappa[],sigma[];

vector hf[];
vector tang[];
scalar nthf[];
scalar ntsd[];
vector nrm[];
vector tang_d[];


void get_tangent(vector n, vector t)
{
  foreach(){
    t.x[] = n.y[] != nodata ? -n.y[] : nodata;
    t.y[] = n.y[] != nodata ? n.x[] : nodata;
  }
}
  
int main(){
  init_grid (64);
  L0 = 8;
  Y0 = -L0/8.;
  X0 = -L0/2.;
  fraction(f,(y-((double)L0)/((double)N)/2. + 0.0*sin(2*pi*x/8)));
  //fraction(f,sq(0.5)-(sq(x) + sq(y)));
  heights(f,hf);
  /* foreach(){ */
  /*   hf.x[] = hf.x[] == nodata ? 0. : hf.x[]; */
  /*   hf.y[] = hf.y[] == nodata ? 0. : hf.y[]; */
  /*   fprintf(ferr,"%g %g %g %g\n",x,y,hf.x[],hf.y[]); */
  /* } */
    
  foreach(){
    coord nrmf = height_normal(point,f,hf);
    foreach_dimension()
      nrm.x[] = nrmf.x;
  }
  get_tangent(nrm,tang);

  vof2distCC(f);
  foreach(){
    foreach_dimension()
      tang_d.y[] = (d[1] - d[-1])/2./Delta;    
  }

  foreach(){
    if(tang.x[] != nodata && tang.y[] != nodata){
      nthf[] = norm(tang);
      ntsd[] = norm(tang_d);
    }
  }
  
  hf.x.nodump = true;
  dump("toto");
}
