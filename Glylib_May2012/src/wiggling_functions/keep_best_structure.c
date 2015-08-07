#include <molecules.h>
#include <glylib.h>
#include <mylib.h>

/* Load in best_structure and if it clashes less than current one keep it */

void keep_best_structure(assembly *A, assembly *B, double *total_clash){
double prev_total_clash=0;
assembly C;
char filename[20];
strcpy(filename, "best_structure.pdb");
load_pdb(filename,&C);
int ri;

prev_total_clash=*total_clash;
*total_clash=find_vdw_clashes_return_total(&C,B,1);
if (*total_clash<=prev_total_clash){ // if best structure found for resolving resid clash does not increase total_clash
	printf("Keeping best structure\n");
	//deallocateAssembly(A);
	load_pdb(filename,A); // replace A with it and continue
	for (ri=0;ri<(*A).m[0][0].nr;ri++){
       		if ((*A).m[0][0].r[ri].na>1){
               		set_smallest_rings_from_residue_atom_nodes(&(*A).m[0][0].r[ri]);
       		}
	}
	//set_nbonds_for_atoms_in_assembly(A);
	//set_nmb_for_atoms_in_assembly(A);
}
//deallocateAssembly(&C);
return;
}
