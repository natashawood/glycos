#include <glylib.h> // required for load_pdb function
#include <mylib.h>
#include <molecules.h>

//int main (void){
double wiggler(assembly *A, assembly *B, int *core_ri){

int window_size=5; // fineness of search
int branch[(*A).nr]; // stores which residues are in the branch
branch[0]=0; // stores number of residues in branch
int linkage_resid=0,linked2=0,current_resid=0; // linkage is resid being wiggled. current_resid is one with clash being resolved
double resid_clash=1000; 
double total_clash=1000; 
int struct_total=10000; // id for structures generated so far. loads into vmd in order when starting at 10000
int i=1; // 
int y=0;
int lna;
int tors=0;// dummy for the first call, set in generate_torsion...
int stop=0;

dihedral_coord_set *DCS;
DCS=(dihedral_coord_set*)calloc(8, sizeof(dihedral_coord_set));
//DCS=(dihedral_coord_set*)realloc(DCS, 8*sizeof(dihedral_coord_set)); // space for 8


while ( (total_clash>36) && (struct_total<11000) && (stop==0) ) {
	resid_clash=1000; //reset after resolved each one
	total_clash=find_vdw_clashes_return_total(A,B,1);
	printf("total_clash=%2f\n",total_clash);
	current_resid=find_vdw_clashes_return_clashiest_resid(A,B,1); // SORT OUT THAT 1! Not doing anything at the minute.
	for(y=2;y<=core_ri[0];y++){ //The converted ROH residue will be in core_ri[1]. Any residue attached there needs to be wiggled.
                printf("core_ri[%d]=%d\n",y,core_ri[y]);
		if(current_resid==core_ri[y]){
			stop=1;}} // if current_resid has already become a core then stop
	find_path_to_core(A, core_ri, current_resid, branch);
	i=1;
  //      printf("branch[%d]=%d,branch[%d+1]=%d\n",i,branch[i],i+1,branch[i+1]);
	while ( (resid_clash>1) && (i<=branch[0]) && (stop==0) ) {
		printf("Starting wiggling\n");
//		printf("branch[%d]=%d,branch[%d+1]=%d\n",i,branch[i],i+1,branch[i+1]);
		linkage_resid=branch[i+1]; // this will be rotated branch
		linked2=branch[i]; // this is the core

		printf("Will find which atoms are involved in linkage\n");
		lna=find_connection_atoms(linkage_resid, linked2, DCS, A);
		tors=(lna-3);

                // if tors>0 so won't try to wiggle OME or SO3. Can't atm. 
                if (tors>0){
                     printf("Will generate torsions\n");
                     generate_torsion_window_permutations(window_size,A,B,tors,&struct_total,&resid_clash,current_resid,DCS,lna,linked2);
                     printf("Finished generating torsions for resid %d\n",current_resid);   
                     keep_best_structure(A,B,&total_clash); // otherwise move onto next linkage and see if we can fix it from there.   
                }
		i++;
                // Add wiggled bits to core. They won't be wiggled again.
		core_ri[0]=((core_ri[0])+1);
		core_ri[core_ri[0]]=linkage_resid;
		for (y=1;y<core_ri[0];y++){
			printf("core_ri[%d]=%d\n",y,core_ri[y]);
		}
	}
}
return total_clash;
}
