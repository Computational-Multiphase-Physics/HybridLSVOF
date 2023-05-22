#include "fractions.h"
#include "curvature.h"
#include "distance.h"

#define FRONTTHR 1e-6
#define BGHOSTS 2

// interface normal computed with height function
coord interface_normal2(Point point, scalar c)
{
  coord n;
  if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata)
    n = mycs (point, c);
  else {
    double nn = 0.;
    foreach_dimension()
      nn += fabs(n.x);
    foreach_dimension()
      n.x /= nn;
  }
  return n;
}

// reconstruction function from fractions but use interface normal computed with height function

void reconstruction2(const scalar c, vector n, scalar alpha)
{
  foreach() {

    /**
    If the cell is empty or full, we set $\mathbf{n}$ and $\alpha$ only to
    avoid using uninitialised values in `alpha_refine()`. */

    if (c[] <= 0. || c[] >= 1.) {
      alpha[] = 0.;
      foreach_dimension()
	n.x[] = 0.;
    }
    else {

      /**
      Otherwise, we compute the interface normal using the
      Mixed-Youngs-Centered scheme, copy the result into the normal field
      and compute the intercept $\alpha$ using our predefined function. */

      coord m = interface_normal2(point, c);
      foreach_dimension()
	n.x[] = m.x;
      alpha[] = plane_alpha (c[], m);
    }
  }

#if TREE

  /**
  On a tree grid, for the normal to the interface, we don't use
  any interpolation from coarse to fine i.e. we use straight
  "injection". */

  foreach_dimension()
    n.x.refine = n.x.prolongation = refine_injection;

  /**
  We set our refinement function for *alpha*. */

  alpha.n = n;
  alpha.refine = alpha.prolongation = alpha_refine;
#endif

  /**
  Finally we apply the boundary conditions to define $\mathbf{n}$ and
  $\alpha$ everywhere (using the prolongation functions when necessary
  on tree grids). */

  boundary ({n, alpha});
}

void vof2dist(scalar c, vertex scalar phi){
  //if phi.height exists compute height function
  if (c.height.x.i)
    heights (c, c.height);

  // compute interface reconstruction with height function normal
  scalar alpha_front[];
  vector n_front[];

  reconstruction2(c, n_front, alpha_front);


  //set magnitude of signed distance to be greater than maximum domain size
  foreach_vertex()
    phi[] = L0*sqrt(2)+1;

  boundary({phi});

#if dimension == 2
  foreach_vertex(){
    double vertex_signed_distance = phi[];
    coord vertex_c = (coord){x,y,0};
    //copied from foreach_neighbor looks at 4 by 4 stencil centered around vertex
    int _s = 2;
    int _nn = _s + 0 ? _s + 0 : GHOSTS;
    int _i = point.i, _j = point.j;
    for (int _k = - _nn; _k <= _nn-1; _k++) {
      point.i = _i + _k;
      for (int _l = - _nn; _l <= _nn-1; _l++) {
        point.j = _j + _l;
        POINT_VARIABLES;
        if (fabs(c[] - 0.5) <= (0.5-FRONTTHR)){
          //fprintf (fout, "interfacial at %g, %g\n ", x, y);
          coord n = {n_front.x[], n_front.y[], 0.0};
          //compute alpha and facets
          coord segment[2];
          coord temp1;
          coord temp2;
          if (facets (n, alpha_front[], segment) == 2){
            //define segments in absolute location
            segment[0] = (coord){x + segment[0].x*Delta, y + segment[0].y*Delta, 0.0};
            segment[1] = (coord){x + segment[1].x*Delta, y + segment[1].y*Delta, 0.0};
            //define segment such that normal matches
            if ((sign(segment[1].y -segment[0].y) == sign(n.x)) && (sign(segment[0].x -segment[1].x) == sign(n.y))){
               temp1 = segment[1];
               temp2 = segment[0];
            }
            else{
              temp1 = segment[0];
              temp2 = segment[1];
            }
            coord r;
            double s, d2 = PointSegmentDistance (&vertex_c, &(temp1), &(temp2), &r, &s);
            if (sqrt(d2) < fabs(vertex_signed_distance)){
              vertex_signed_distance = sqrt(d2)*((double)PointSegmentOrientation(&vertex_c, &(temp1), &(temp2)));}
          }
        }
      }
    }
    point.i = _i; point.j = _j;
    POINT_VARIABLES;
    phi[] = vertex_signed_distance;
  }
#else //dimension ==3


foreach_vertex(){
  double vertex_signed_distance = phi[];
  coord vertex_c = (coord){x,y,z};
  //copied from foreach_neighbor looks at 4 by 4 stencil centered around vertex
  int _s = 2;
  int _nn = _s + 0 ? _s + 0 : GHOSTS;
  int _i = point.i, _j = point.j, _k = point.k;
  for (int _a = - _nn; _a <= _nn-1; _a++) {
    point.i = _i + _a;
    for (int _b = - _nn; _b <= _nn-1; _b++) {
      point.j = _j + _b;
      for (int _c = - _nn; _c <= _nn-1; _c++) {
        point.k = _k + _c;
        POINT_VARIABLES;
        if (fabs(c[] - 0.5) <= (0.5-FRONTTHR)){
          //fprintf (fout, "interfacial at %g, %g\n ", x, y);
              coord n = {n_front.x[], n_front.y[], n_front.z[]};
              coord v[12];
              int m = facets (n, alpha_front[], v, 1.); //m is the number of coordinates
              for (int j = 0; j < m - 2; j++) {
                coord temp1 = {x + v[0].x*Delta  , y + v[0].y*Delta  , z + v[0].z*Delta};
                coord temp2 = {x + v[j+1].x*Delta, y + v[j+1].y*Delta, z + v[j+1].z*Delta};
                coord temp3 = {x + v[j+2].x*Delta, y + v[j+2].y*Delta, z + v[j+2].z*Delta};

                double s, t, d2 = PointTriangleDistance (&vertex_c, &(temp1), &(temp2), &(temp3), &s, &t);
                if (sqrt(d2) < fabs(vertex_signed_distance)){
                  vertex_signed_distance = sqrt(d2)*((double)PointTriangleOrientation(&vertex_c, &(temp1), &(temp2), &(temp3)));}
                }
            }
          }
        }
      }
  point.i = _i; point.j = _j;  point.k = _k;
  POINT_VARIABLES;
  phi[] = vertex_signed_distance;
  }
#endif

  boundary({phi});

  /*heuristic check to confirm sign of non interfacial cells
  foreach vertex check vof values of adjacent cells at most refined level.
  If any are non interfacial set the sign of the vertex such that it is consistent
  */

  #if dimension == 2
  foreach_vertex()
    for (int ii = -1; ii<=0; ii++)
      for (int jj = -1; jj<=0; jj++){
        if is_refined(neighbor(ii,jj)){
          if (fabs(fine(c,ii,jj) - 0.5) > (0.5-FRONTTHR)) //if !interfacial
            phi[] = sign(fine(c,ii,jj) - 0.5)*fabs(phi[]);}
        else{
          if (fabs(c[ii,jj] - 0.5) > (0.5-FRONTTHR)) //if !interfacial
            phi[] = sign(c[ii,jj] - 0.5)*fabs(phi[]);}
      }

  #else //dimension == 3
  foreach_vertex()
    for (int ii = -1; ii<=0; ii++)
      for (int jj = -1; jj<=0; jj++)
        for (int kk = -1; kk<=0; kk++){
          if is_refined(neighbor(ii,jj,kk)){
            if (fabs(fine(c,ii,jj,kk) - 0.5) > (0.5-FRONTTHR)) //if !interfacial
              phi[] = sign(fine(c,ii,jj,kk) - 0.5)*fabs(phi[]);}
          else{
            if (fabs(c[ii,jj,kk] - 0.5) > (0.5-FRONTTHR)) //if !interfacial
              phi[] = sign(c[ii,jj,kk] - 0.5)*fabs(phi[]);}
        }

  #endif
  boundary({phi});
}

/*function that modifies signed distance where the cell fractions and face
fractions would produce a degenerate case due the film thickness being thinner
than the grid size. In these cases the function either perforates the film or thickens
it depending on the VOF fraction */

void dist_cleanup(scalar c, vertex scalar phi){

  scalar num_faces[];

  num_faces.prolongation = num_faces.refine = refine_injection;
  num_faces.restriction = no_restriction;

  foreach(){
    num_faces[] = 0;
    num_faces[] += (((phi[0,0])*(phi[1,0]) < 0.) ? 1 : 0);
    num_faces[] += (((phi[0,0])*(phi[0,1]) < 0.) ? 1 : 0);
    num_faces[] += (((phi[1,0])*(phi[1,1]) < 0.) ? 1 : 0);
    num_faces[] += (((phi[0,1])*(phi[1,1]) < 0.) ? 1 : 0);
  }

  boundary({num_faces});

  foreach_vertex(){
    bool to_change = false;
    double f_cell;
    for (int ii = -1; ii<=0; ii++)
      for (int jj = -1; jj<=0; jj++){
        if is_refined(neighbor(ii,jj)){
          if (fine(num_faces,ii,jj) > 2){
            to_change = true;
            f_cell = fine(c,ii,jj);}
        }
        else{
        if (num_faces[ii,jj] > 2){
          to_change  = true;
          f_cell = c[ii,jj];}
        }
      }
    if (to_change) {
      if (f_cell < 0.5){
        phi[] = (phi[] >= 0.0 ? -0.00001 : phi[]);}
      else{
        phi[] = (phi[] <= 0.0 ? 0.00001 : phi[]);}
    }
  }
  boundary({phi});
}
