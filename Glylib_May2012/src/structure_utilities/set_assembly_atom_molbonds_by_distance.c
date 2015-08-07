/*
Oliver Grant 7Feb11
This function was written as to fill in number of bonds info when a pdb was read in.
Finds the number of bonds an atom in a structure has based on distance

*/
#include <mylib.h>
#include <molecules.h>
#define DIST 1.7
void set_assembly_atom_molbonds_by_distance(assembly *A){ 
//printf("Setting nmb based on distance between atoms\n");
// revieved a pointer to an assembly need to dereference 

double x,y,z,d; // xyz are holders to make code more readable. d is distance between two atoms
//compare dist between each atom and all other atoms in that residue
(*A).nb=0; // Set this as zero intially. Using it as total number of bonds in the molecule

/* Found a strange instance where nmb wasn't initialized to zero. This fixed it. OG
for(mi=0;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
		(*A).m[mi][0].r[ri].nrb=0;
		for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){
			(*A).m[mi][0].r[ri].a[ai].nmb=0;
		}
	}
}
*/

int ai=0,aii=0,nmb=0,bi=0;

for(ai=0;ai<(*A).na;ai++){// ai is now current atom to check every other atom against
	for(aii=0;aii<(*A).na;aii++){ //for current atom check dist to every other atom
        	x=((*A).a[aii][0].x.i-(*A).a[ai][0].x.i);
                y=((*A).a[aii][0].x.j-(*A).a[ai][0].x.j);
                z=((*A).a[aii][0].x.k-(*A).a[ai][0].x.k);
                d=sqrt((x*x)+(y*y)+(z*z)); //pythagorus square on hyp in 3D
                //printf("\ndistance is %f\n",d);
                if (d<=DIST && d>0.001){ // If not the same atom but close enough for a bond
                	//increment number of bonds counter if atom within DIST but not itself (d>0.001)
                        (*A).a[ai][0].nmb++; /* shorten this for readability */
                        (*A).nb++;
			nmb=(*A).a[ai][0].nmb;  /* so the next bits don't get really long */
			/* allocate or reallocate space to accommodate the new bond info */
			if(nmb==1){(*A).a[ai][0].mb=(molbond*)calloc(1, sizeof(molbond));}
			else{(*A).a[ai][0].mb=(molbond*)realloc((*A).a[ai][0].mb,nmb*sizeof(molbond));}
			//for(bi=0;bi<nmb;bi++){
			//	A[0].a[ai][0].mb[bi].s.a=ai;
					
			bi=(nmb-1); // always one less than number of bonds so that index starts at 0
                        //set source atom info, redundant in this case but hey it's there and I can fill it in 
			A[0].a[ai][0].mb[bi].s=A[0].a[ai][0].moli;
                        A[0].a[ai][0].mb[bi].t=A[0].a[aii][0].moli;
			/*if(A[0].a[ai][0].mb[bi].s.r!=A[0].a[ai][0].mb[bi].t.r){
				printf("bonding ri=%d,a=%s to ri=%d,a=%s\n",A[0].a[ai][0].mb[bi].s.r,A[0].a[ai][0].N,A[0].a[ai][0].mb[bi].t.r,A[0].a[aii][0].N);
			}*/ //good for debugging residue bonds
		}
	}
}
return;
}
