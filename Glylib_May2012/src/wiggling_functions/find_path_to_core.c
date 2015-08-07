#include <glylib.h>
#include <mylib.h>
#include <molecules.h>

// Function Declaration
void find_branch_to_core(int *core_ri, int c_ri, int c_mi, int pri, assembly *A, int *found_core, int *branch);

// CODE START
void find_path_to_core(assembly *A, int *core_ri, int current_resid, int *branch){
int mi=0,ri=0; // indexes
int c_ri=0,c_mi=0; // Connected residue index, connected atom index
int i=0, found_core=0;
//int core_dir=0; // direction from the linkage that the core lies

branch[0]=0; // reset each time function is called.
// first find the mb which leads to the core
for(mi=0;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
                if(current_resid==(*A).m[mi][0].r[ri].n){
/*DEBUG*/              printf("Found linkage residue %d\n",current_resid); 
                        if((*A).m[mi][0].r[ri].nrb>0){ // if it has residues connected (It will)
/*DEBUG*/                       printf("It has %d connected residues\n",(*A).m[mi][0].r[ri].nrb);
                                for(i=0;i<(*A).m[mi][0].r[ri].nrb;i++){ //for each residue connected
                                        c_mi=(*A).m[mi][0].r[ri].rb[i].t.m;
                                        c_ri=(*A).m[mi][0].r[ri].rb[i].t.r;
                                        current_resid=(*A).m[mi][0].r[ri].n;
//*DEBUG*/                               printf("BACK OUT TO FIRST CALL\nResid connected to linkage=%d\n",(*A).m[c_mi][0].r[c_ri].n);
	                                find_branch_to_core(core_ri, c_ri, c_mi, current_resid, A, &found_core, branch); // self calling function
                                        if(found_core==1){
/*DBG*/                                       	printf("Resid %d leads to the core\n",(*A).m[c_mi][0].r[c_ri].n);
						branch[0]=(branch[0]+1);
						branch[branch[0]]=(*A).m[mi][0].r[ri].n;
                                                found_core=0; //reset
                                        }
                                }
                        }
                }
        }
}

for(i=0;i<=branch[0];i++){
	printf("branch[%d]=%d\n",i,branch[i]);
}

return;
}



// From the current residue c_ri find the direction which leads to the core.
void find_branch_to_core(int *core_ri, int c_ri, int c_mi, int pri, assembly *A, int *found_core, int *branch){
//pri is previous residue id
//c_ri is the current residue id
int i=0;
int cc_ri=0, cc_mi=0; // connected to current indexes
int nrb; //readability. Number of residue bonds to other residues
int cri; // current residue id

cri=(*A).m[c_mi][0].r[c_ri].n;
nrb=(*A).m[c_mi][0].r[c_ri].nrb;

//DEBUG printf("\nExploring resid %d which has %d neighbours\n",cri,nrb);

for (i=2;i<=core_ri[0];i++){ // the array size is stored in [0]
//*DBG*/	printf("cri=%d, core_ri[%d]=%d\n",cri,i,core_ri[i]);
        if (cri==core_ri[i]){ //FOUND_CORE!
                *found_core=1;
/*DEBUG*/       printf("Found core, resid=%d\n",core_ri[i]);
		branch[0]=((branch[0])+1);
                branch[branch[0]]=(*A).m[c_mi][0].r[c_ri].n;
        }
}

// if this isn't the core and it has connections which aren't the previous residue
if ( (*found_core==0) && (nrb>1) ){ //cc_ri can equal c_ri as long as mi is different
//*DBG*/printf("Exploring resid %d's %d neighbours\n",cri,nrb);
        for (i=0;i<nrb;i++){ // for each neighbour
                cc_ri=(*A).m[c_mi][0].r[c_ri].rb[i].t.r;
                cc_mi=(*A).m[c_mi][0].r[c_ri].rb[i].t.m;
//*DEBUG*/       printf("loop %d, Currently checking resid=%d\n",i,(*A).m[cc_mi][0].r[cc_ri].n);
                if (pri!=(*A).m[cc_mi][0].r[cc_ri].n) {
//*DEBUG*/               printf("Not the same residue as previous:%d\n ... so will check it out\n",pri);
                        if (*found_core==0){
				find_branch_to_core(core_ri, cc_ri, cc_mi, cri, A, found_core, branch);
			}
                } else {
//*DEBUG*/               printf("But it's the same residue as before=%d\n",pri);
                }
        }
	if((*found_core)==1){ // we are dropping out of the functions and have found the core
 	      	branch[0]=((branch[0])+1);
        	branch[branch[0]]=(*A).m[c_mi][0].r[c_ri].n;
//*DBG*/      	printf("Resid %d leads to the core\n",(*A).m[c_mi][0].r[c_ri].n);
		for (i=0;i<=branch[0];i++){
//*DBG*/       		printf("branch[%d]=%d\n",i,branch[i]);
		}
	}
}

//DEBUG printf("Returning found core=%d\n",*found_core);
return;
}
