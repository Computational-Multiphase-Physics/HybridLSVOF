/**
  this header file is to provid the level-set based multiphase solver
*/

const double width = 1.; //interfacial region 

struct LS_reinit {
  scalar dist;
  double dt;
  int it_max;
};

double sign2 (double x)
{
  return(x > 0. ? 1. : x<0 ? -1. : 0.);
}

double my_minmod(double a,double b){
  if(a*b>0){
    if(fabs(a) < fabs(b))
      return a;
    if(fabs(a) > fabs(b))
      return b;
  }
  return 0;
}

double delta_finest ()
{
  return L0 / pow(2,grid->maxdepth);
}

double vofLS(double phi, double epsilon)
{
  double H;
#if 1
  //smeared heavside function
	if (phi < -epsilon)
		H = 0.0;
	else if (phi > epsilon)
		H = 1.0;
	else
		H = 0.5 * (1.0 + phi / epsilon + 1.0 / pi * sin(pi * phi / epsilon));
#else
  //fixme: the one used in the paper of AL Saud. There is something wrong with my implementation.
  //I leave it for now as I think it's not important that much
  H = 1.0 - 0.5 * (1.0 - erf(phi / (epsilon)));
#endif

	return H;
}

foreach_dimension()
static inline double WENOdiff_x(Point point, scalar s, int i){
  double s1 = (s[2*i,0,0] + s[] - 2*s[i,0,0])/Delta; 
  double s2 = (s[1,0,0] + s[-1,0,0] - 2*s[])/Delta;
  return i*((s[i,0,0] - s[])/Delta -my_minmod(s1,s2)/2.);
}

void prehamil(Point point, coord  * grapl, coord * gramin, scalar s){
  foreach_dimension(){
    grapl->x  = WENOdiff_x(point,s,1);
    gramin->x = WENOdiff_x(point,s,-1);
  }
}

foreach_dimension()
static inline double root_x(Point point, scalar s, double eps, int dir)
{
  // dir == 1 or -1 offsets the position of the interface
  double phixx = my_minmod(s[2*dir] + s[] - 2*s[dir], 
                           s[1] + s[-1] - 2*s[]);
  if(fabs(phixx) > eps){
    double D = sq(phixx/2.-s[] - s[dir])-4*s[]*s[dir];
    // fprintf(stderr, "%g %g %g\n", D, phixx, eps);
    return 1/2.+( s[] - s[dir] - sign2(s[] - s[dir])*sqrt(D))/phixx;
  }
  else{
    return s[]/(s[]- s[dir]);
  }
}

double hamiltonian (Point point, scalar s0, coord grapl, coord  gramin)
{
  double hamil = 0;
  if(s0[] > 0){
    foreach_dimension(){
      double a = min(0.,grapl.x); 
      double b = max(0.,gramin.x);
      hamil += max(sq(a),sq(b));
    }
    return sqrt(hamil);
  }
  else{
    foreach_dimension(){      
      double a = max(0.,grapl.x);
      double b = min(0.,gramin.x);
      hamil += max(sq(a),sq(b));
    }
  }
  return sqrt(hamil);
}

double ForwardEuler(scalar dist, scalar temp, scalar dist0, double dt){
  double res=0.;

  foreach(reduction(max:res)){
    double delt =0.;

    double flag = 1.;
    foreach_dimension(){
      flag = min (flag, dist0[-1]*dist0[]);
      flag = min (flag, dist0[ 1]*dist0[]);
    }

    coord grapl, gramin;
    prehamil(point, &grapl, &gramin, temp);

    if(flag < 0.){ // the cell contains the interface
      double size = 1.e10;
      foreach_dimension(){
        if(dist0[]*dist0[1]<0){
          double dx = Delta*root_x(point, dist0, 1.e-10, 1);
          double sxx1 = (temp[2] + temp[] - 2*temp[1])/sq(Delta);
          double sxx2 = (temp[1] + temp[-1] - 2*temp[])/sq(Delta);
          if(dx !=0.)
            grapl.x = -temp[]/dx - dx* my_minmod(sxx1,sxx2)/2.;
          else 
            grapl.x = 0.;
          size = min(size, dx);
        }
        if(dist0[]*dist0[-1]<0){
          double dx = Delta*root_x(point, dist0, 1.e-10, -1);
          // if(dx>10.){
            // fprintf(stderr, "%g %g %g\n", dist0[],dist0[-1,0],dx);
          // }
          double sxx2 = (temp[1] + temp[-1] - 2*temp[])/sq(Delta);
          double sxx3 = (temp[-2] + temp[0] - 2*temp[-1])/sq(Delta);
          if(dx!=0.)
            gramin.x = temp[]/dx + dx* my_minmod(sxx3,sxx2)/2.;
          else 
            gramin.x = 0.;
          size = min(size, dx);
        }
      }
      delt = sign2(dist0[]) * min(dt,fabs(size)/2.) *
      (hamiltonian(point, dist0, grapl,gramin) - 1);
      dist[] -= delt;
    }
    else{ 
      delt = sign2(dist0[]) * 
      (hamiltonian(point, dist0, grapl, gramin)- 1);
      dist[] -= dt*delt;
    }
    res = max (res,fabs(delt));
  }


  boundary({dist});
  restriction({dist});

  return res;
}

/**
## LS_reinit() function
*/
int LS_reinit(struct LS_reinit p){
  scalar dist = p.dist; 

  double dt   = p.dt;     // standard timestep (0.5*Delta)
  int it_max  = p.it_max;// maximum number of iteration (100)

/**
In 2D, if no specific timestep is set up by the user, we take the most
restrictive one
with regards to the CFL condition :
$$
\Delta t = 0.5 * \Delta x
$$
*/
  
  if(dt == 0) dt = 0.5 * L0/(1 << grid->maxdepth);


/**
Default number of iterations is 20 times, which is sufficient to have the first
10 neighbor cells to the 0-level-set properly redistanced.
*/

  if(it_max == 0)it_max = 5;

  vector gr_LS[];
  int i ;

/**
Convergence is attained is residual is below $dt\times 10^{-6}$
*/  
  double eps = dt*1.e-6;

/**
We create `dist0[]` which will be a copy of the initial level-set function
before the iterations and `temp[]` which will be $\phi^{n}$ used for the
iterations.
*/
  scalar dist0[];
  foreach(){
    dist0[] = dist[] ;
  }
  boundary({dist0});

/**
 Time integration iteration loop.

One can choose between Runge Kutta 2 and Forward Euler temporal integration.
*/
  for (i = 1; i<=it_max ; i++){
    double res = 0;

/**

## RK3
We use a Runge Kutta 3 compact version taken from [Shu and Osher](#Shu1988)
made of 3 backward Euler steps:

* Step1-2
$$
\frac{\widetilde{\phi}^{n+1}  - \phi^n}{\Delta t}  = \text{RHS}^n\\
\dfrac{\widetilde{\phi}^{n+2}  - \widetilde{\phi}^{n+1}}{\Delta t}  = \widetilde{RHS}^{n+1} 
$$
with :
$$
RHS =  sgn (\phi_{ij}^0)\cdot \left[ H_G\left( D_x^+\phi_{ij}^n, D_x^-\phi_{ij}^n, D_y^+\phi_{ij}^n,
D_y^-\phi_{ij}^n \right)\right]
$$
*/
    scalar temp[],temp1[], temp2[];
    foreach(){
      temp[] = dist[] ;
      temp1[] = dist[] ;
    }
    boundary({temp,temp1});
    ForwardEuler(temp1,temp,dist0,dt);
    foreach(){
      temp2[] = temp1[] ;
    }
    boundary({temp2});
    ForwardEuler(temp2,temp1,dist0,dt);
/**
* Intermediate value
$$
\widetilde{\phi}^{n+1/2}  = \dfrac{3}{4}\widetilde{\phi}^{n} + \dfrac{1}{4}\widetilde{\phi}^{n+2}
$$
*/
    foreach(){
      temp1[] = 3./4*dist[] + temp2[]/4.;
      temp2[] = temp1[];
    }
    boundary({temp1,temp2});

/**
* Step 3
$$
\widetilde{\phi}^{n+3/2} - \widetilde{\phi}^{n+1/2} = \widetilde{RHS}^{n+1/2}
$$
*/
    ForwardEuler(temp2,temp1,dist0,dt);
/**
* Final Value
$$
\widetilde{\phi}^{n+1} = \widetilde{\phi}^{n} + \dfrac{2}{3}\widetilde{\phi}^{n+3/2}
$$
*/
    foreach(reduction(max:res)){
      res = max(res, 2./3.*fabs(dist[] - temp2[]));
      dist[] = dist[]/3. + temp2[]*2./3.;
    }
    boundary({dist});
    restriction({dist});
/**
Iterations are stopped when $L_1 = max(|\phi_i^{n+1}-\phi_i^n|) < eps$
*/
    if(res<eps){
      return i;
    }
  }
  return it_max;
}


foreach_dimension() 
static double WENO5_x(Point point, scalar q, int i)
{
    static const double coeff[3][3] = {
        {1. / 3., -7. / 6., 11. / 6.},
        {-1 / 6., 5. / 6., 1. / 3.},
        {1. / 3., 5. / 6., -1. / 6.}};
    static const double wIS = 13. / 12;
    static const double weights[3] = {0.1, 0.6, 0.3};
    static const double eps = 1.e-12;

    double ENO3[3], alpha[3], IS[3];
    ENO3[0] = q[i * 2] * coeff[0][0] + q[i] * coeff[0][1] + q[] * coeff[0][2];
    ENO3[1] = q[i] * coeff[1][0] + q[] * coeff[2][1] + q[-i] * coeff[1][2];
    ENO3[2] = q[] * coeff[2][0] + q[-i] * coeff[2][1] + q[-2 * i] * coeff[2][2];

    IS[0] = wIS * sq(q[2 * i] - 2 * q[i] + q[]) + sq(q[2 * i] - 4 * q[i] + 3 * q[]) * 0.25;
    IS[1] = wIS * sq(q[i] - 2 * q[] + q[-i]) + sq(q[i] - q[-i]) * 0.25;
    IS[2] = wIS * sq(q[] - 2 * q[-i] + q[-2 * i]) + sq(3 * q[] - 4 * q[-i] + q[-2 * i]) * 0.25;

    alpha[0] = weights[0] / sq(eps + IS[0]);
    alpha[1] = weights[1] / sq(eps + IS[1]);
    alpha[2] = weights[2] / sq(eps + IS[2]);

    double sum = alpha[0] + alpha[1] + alpha[2];

    return (alpha[0] * ENO3[0] + alpha[1] * ENO3[1] + alpha[2] * ENO3[2]) / sum;
}

void advectLevelSet(scalar phi, vector uadv, double dt, int i)
{
    vector upfluxp[], upfluxm[];
    foreach(){
        foreach_dimension(){
        upfluxp.x[] = (phi[1] - phi[]);
        upfluxm.x[] = (phi[] - phi[-1]);
        }
    }
    boundary((scalar *){upfluxm,upfluxp});
    scalar dvar[];
    foreach()
    {
        dvar[] = 0.0;
        foreach_dimension()
        {
            double dphi = 0.0;
            double vel = uadv.x[];
            dphi = vel > 0 ? WENO5_x(point, upfluxm.x, -1) : WENO5_x(point,upfluxp.x, 1);
            dphi = dphi / Delta * vel * dt;
            dvar[] += dphi ;
        }
    }

    foreach()
    {
        phi[] -= dvar[];
    }

#if NOREINIT == 1

#else
    if(i % 10 == 0)
    {
      LS_reinit(phi, it_max = 10);
    }
#endif
    
    boundary({phi});
}
