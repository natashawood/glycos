#include <glylib.h> // for load_pdb function
#include <molecules.h>
#include <mylib.h>

// takes in assembly, name to call the file, list of residue names and number of residues in the list
void mergeAsmbly(assembly* merged, assembly* A, assembly *B){
//Tells the user what's going on
printf("Entered Merge Assembly function...\n");

char file_name[50];
strcpy(file_name,"temp.pdb");

outputMolPDB_NoRenameResid(&(*A).m[0][0], file_name); //Glycan should be all in one molecule from now on. outputAsmbl now adds TER cards and needs to to work
//outputAsmblPDB(A, file_name); // overwrites temp file and put A in there
appendAsmblPDB(B, file_name); // adds B onto the bottom of A. Written by Oliver so no TER cards need to be added for it to work.
load_pdb(file_name,merged); // loads in this pdb file into one assembly
int ri=0;
for (ri=0;ri<(*merged).m[0][0].nr;ri++){
        if ((*merged).m[0][0].r[ri].na>1){
                set_smallest_rings_from_residue_atom_nodes(&(*merged).m[0][0].r[ri]);
        }
}
return;
}
