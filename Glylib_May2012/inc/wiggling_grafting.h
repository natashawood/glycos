/** \file wiggling_grafting.h
\brief  Header file for wiggling and grafting

	by Oliver Grant. */

#include <geometries.h>
#include <parameter_sets.h>
#include <molecules.h>
#include <mylib.h>

void set_nmb_for_atoms_in_assembly(assembly *A);
void change_torsion_by(double delta_phi, dihedral_coord_set *DCS, int lna, int tors, coord_3D **brnch, int *nab);
void generate_torsion_window_permutations(int window_size, assembly *A, assembly *B, int tors, int *total_struct, double *lowest_clash, int linkage_resid, dihedral_coord_set *DCS, int lna, int linked2);
int find_vdw_clashes_return_clashiest_resid(assembly *A, assembly *B, int nbr);
double find_vdw_clashes_return_total(assembly *A, assembly *B, int nbr);
double find_vdw_clashes_return_resid_total(assembly *A, assembly *B, int *current_resid, int nbr);
int find_next_resid(assembly *A, int current_resid, int *core_ri);
void find_path_to_core(assembly *A, int *core_ri ,int current_resid, int *branch);
void keep_best_structure(assembly *A, assembly *B, double *total_clash);
double wiggler(assembly *A, assembly *B, int *core_ri);
void explore_branch_for_core(int *core_ri, int c_ri, int c_mi, int pri, assembly *A, int *found_core);
int find_connection_atoms(int linkage, int linked2, dihedral_coord_set *DCS, assembly *A);
void find_branch_to_rotate(assembly *A, int tors, coord_3D **brnch, int *nab, dihedral_coord_set *DCS, int lna, int linked2);
int find_vdw_clashes_return_self_clashes(assembly *A); // not used
int find_bonds_meaning_clash(assembly *A, assembly *B);
int find_new_bonds_meaning_clash(assembly *A);

//grafting
void mergeAsmbly(assembly* merged, assembly* A, assembly *B);
void appendAsmblPDB(assembly* asmbl, char* file_name);
void superimpose4atoms(assembly *A, assembly *B);
void outputResPDB(assembly* asmbl, char* file_name, int* residues, int outrn);
vectormag_3D coord2_to_vec(coord_3D c1, coord_3D c2);
void rename_superimpose(assembly *A, assembly *B, char *core_aN, char *branch_aN, int core_rn, int branch_rn, int Tcore_rn, int Tbranch_rn);
void rename_core(assembly *A, assembly *B, int TnumRes, int LGnumRes, int z, int fTtoLG[][40], int *core_ri);
double find_vdw_clashes_pairwize_between_Assemblies(assembly *A, assembly *B, char *fileprefix, int *core_ri, char wiggle);
int check_surface_plane(assembly *A, assembly *B);
int graft_spacer(assembly *A, assembly *B, char *fileprefix);
int find_vdw_clashes_spacers(assembly *A, assembly *B, char *spacername, char *spaceroutname, int *fitterCnt);
