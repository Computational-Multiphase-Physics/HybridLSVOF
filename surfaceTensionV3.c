#define GRAVITY 0 // make it 1 if you want to add gravity (parallel to axis)
#define FILTERED 0
//#include "axi.h"
//#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension-hybrid.h"

int MAXLEVEL = 9;
#define MINLEVEL (MAXLEVEL-4)

#define g_accel 0.0
#define R0 1.0
#define Ldomain (16.0*R0)
#define centerx (0.5*Ldomain)
#define centery (0.0*Ldomain)
#define Drop(x,y) (sq(x-centerx) + sq(y-centery))

#define Re 0.1
#define Ca 0.1

#define RHO1 (1.0)
#define RHO2 (1.0)
#define MU1 (1.0/(Re))
#define MU2 (1.0/(Re))

#define SIGMA0 (1.0/(Re*Ca))
#define SIGMAT (1.0/(Re))

double tmax = 6.;
double tsnap = 0.010;
double DTMAX = 0.00005;

/* void careful_refinement(){ */
/*   for(int ii = MINLEVEL ; ii <= MAXLEVEL; ii++){ */
/*     refine(level < ii && sqrt(Drop(x,y)) > (R0 - 4.*sqrt(2.)*Ldomain/(1<<(ii-1))) && sqrt(Drop(x,y)) < (R0 + 4.*sqrt(2.)*Ldomain/(1<<(ii-1)))); */
/*   } */
/* } */

uf.n[left] = neumann(0.0);
uf.n[right] = neumann(0.0);

uf.t[left] = neumann(0.0);
uf.t[right] = neumann(0.0);

p[left] = dirichlet(0.0);
p[right] = dirichlet(0.0);

int main(){
  dtmax = DTMAX;
  TOLERANCE = 1e-6;
  mu1=1./Re;
  mu2=1./Re;
  rho1=1.;
  rho2=1.;
  L0=Ldomain;
  origin(0.0,0.0);
  for(MAXLEVEL = 7; MAXLEVEL <= 7; MAXLEVEL++){
    N=1<<MAXLEVEL;
    run();
    //    system("gnuplot plots");
  }
}

scalar sigma[], kappa[];

event init(i=0){
  if(!restore (file = "restore")){
    //    careful_refinement();
    fraction(f,sq(R0)-(Drop(x,y)));
    foreach(){
      sigma[] = SIGMA0 + SIGMAT*x;
    }
    boundary({sigma});
  }
  dump("init");
}

/* event adapt(i++){ */
/*   adapt_wavelet((scalar *){u,f,KAPPA},(double []){0.00001,0.00001,0.01,0.001},maxlevel = MAXLEVEL, minlevel = MINLEVEL); */
/* } */

event logWriting (i++) {
  double ux = 0;
  double xc = 0;
  double weight = 0;
  foreach(reduction(+:xc) reduction(+:ux) reduction(+:weight)){
    if (f[] >= 0.5){
      weight += f[]*dv();
      ux += f[]*u.x[]*dv();
      xc += f[]*x*dv();
    }
  }
  if (weight > 0){
    xc /= (L0*weight);
    ux /= weight;
  }

  if (i == 0) {
    fprintf (ferr, "#i dt t xc yc U LEVEL = %d\n",MAXLEVEL);
  } else {
    fprintf (ferr, "%d %g %g %g %g %g %d\n", i, dt, t, xc, ux, -2./15.*SIGMAT/mu1, MAXLEVEL);
  }
}

event end (t = tmax){
  /* FILE * fp1; */
  /* char nameGrad[80]; */
  /* sprintf (nameGrad, "surfGrad-%d-%4.3f.dat",MAXLEVEL,t); */
  /* fp1 = fopen (nameGrad, "w"); */
  /* fprintf (fp1, "#x y surfGrad.x surfGrad.y\n"); */
  /* fclose(fp1); */
  /* foreach(){ */
  /*   if(interfacial(point,f)){ */
  /*     fp1 = fopen (nameGrad, "a"); */
  /*     fprintf (fp1, "%g %g %g %g\n", x, y, SurfGradSigma.x[], SurfGradSigma.y[]); */
  /*     fclose(fp1); */
  /*   } */
  /* } */

  char dumpname[80];
  sprintf(dumpname,"dump-%d",MAXLEVEL);
  dump(dumpname);
  
  char namefacet[80];
  sprintf (namefacet, "facet-%d-%d.dat",pid(),MAXLEVEL);
  FILE * fp2 = fopen (namefacet, "w");
  output_facets (f, fp2);
  fclose (fp2);

  char namevec[80];
  sprintf (namevec, "vector-%d-%d.dat",pid(),MAXLEVEL);
  FILE * fpv = fopen (namevec, "w");
  foreach()
    fprintf(fpv,"%g %g %g %g\n",x,y,u.x[],u.y[]);
  fclose (fpv);

}

/**

~~~gnuplot
reset
set term pngcairo
set output "xcvst.png"
set xlabel "Xc"
set ylabel "t"
plot "log" every 10:10 u 3:4:6 w p palette,'' u 3:5

**/
