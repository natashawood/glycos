#include <glylib.h>
#include <mylib.h>
#include <molecules.h>

void change_torsion_by(double delta_phi, dihedral_coord_set *DCS, int lna, int tors, coord_3D **brnch, int *nab){

int i=tors; // tors is torsion index. 1 for phi, 2 for psi and so on

/*********Get vector that describes bond to be rotated **************/
vectormag_3D rot_bond;
rot_bond=get_vector_from_coords(*DCS->X[i+1], *DCS->X[i+2]);
//dprint_vectormag_3D(&rot_bond);

/*********Rotate atoms in brnch around the bond ********************/ 
rotate_coords_about_axis_dp_list(brnch, *nab, *DCS->X[i+1], rot_bond, delta_phi);
//printf("brnch[0].X.i=%f\n",brnch[0].X[1].i);

return;
}
