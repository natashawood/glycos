// Function written by B. Lachele Foley, 2007 and then Oliver Grant 2011
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* this turns a coordinate set into a vector */
vectormag_3D coord2_to_vec(coord_3D c1, coord_3D c2){
vectormag_3D nv ;

nv.i=(c2.i-c1.i);
nv.j=(c2.j-c1.j);
nv.k=(c2.k-c1.k);
// don't trust entry in 'd' -- normalize now to be certain
nv.d=sqrt(nv.i*nv.i+nv.j*nv.j+nv.k*nv.k);

return nv;
}

