/*
Oliver Grant 7Feb11
This function was written as to fill in number of bonds info when a pdb was read in.
Finds the number of bonds an atom in a structure has based on distance
Only finds bonds to atoms within the residue ri. This should change for more general applications but it works well as is for carbohydrate grafting

*/
#include <mylib.h>
#include <molecules.h>
#define DIST 1.7
void set_nbonds_for_atoms_in_assembly(assembly *A){ 

// revieved a pointer to an assembly need to dereference 

int mi=0,ri=0,ai=0,aii=0,bi=0;//loop counters


double x,y,z,d; // xyz are holders to make code more readable. d is distance between two atoms
//compare dist between each atom and all other atoms in that residue
for(mi=0;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
                for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){// ai is now current atom to check every other atom against
                       // printf("\nChecking %s.%d.%s against\n",(*A).m[mi][0].r[ri].N,(*A).m[mi][0].r[ri].n,(*A).m[mi][0].r[ri].a[ai].N);
			for(aii=0;aii<(*A).m[mi][0].r[ri].na;aii++){//for current atom check dist to every other atom within the residue ri
                                x=((*A).m[mi][0].r[ri].a[aii].x.i-(*A).m[mi][0].r[ri].a[ai].x.i);
                                y=((*A).m[mi][0].r[ri].a[aii].x.j-(*A).m[mi][0].r[ri].a[ai].x.j);
                                z=((*A).m[mi][0].r[ri].a[aii].x.k-(*A).m[mi][0].r[ri].a[ai].x.k);
                                d=sqrt((x*x)+(y*y)+(z*z)); //pythagorus square on hyp in 3D
                                //printf("\ndist is %f\n",d);
                         //       printf("%s=%.1f, ",(*A).m[mi][0].r[ri].a[aii].N,d);
                                if (d<=DIST && d>0.001){
                                        //increment number of bonds counter if atom within DIST but not itself (d>0.001)
                                        (*A).m[mi][0].r[ri].a[ai].nb++; /* shorten this for readability */
					bi=(*A).m[mi][0].r[ri].a[ai].nb;  /* so the next bits don't get really long */
					
					/* allocate or reallocate space to accommodate the new bond info */
					if(bi==1){(*A).m[mi][0].r[ri].a[ai].b=(bond*)calloc(1, sizeof(bond));}
					else{(*A).m[mi][0].r[ri].a[ai].b=(bond*)realloc((*A).m[mi][0].r[ri].a[ai].b,bi*sizeof(bond));}
					bi--; // always one less than number of bonds so that index starts at 0
                                        //set source atom info, redundant in this case but hey it's there and I can fill it in 
					//causes a segfault after 1st loop. I don't understand why.
                                        //(*A).m[mi][0].r[ri].a[ai].b[bi].s.r=ri;
					//(*A).m[mi][0].r[ri].a[ai].b[bi].s.a=ai;
					//(*A).m[mi][0].r[ri].a[ai].b[bi].s.m=mi;
					/* set target atom info, mi=mii and ri=rii */
					(*A).m[mi][0].r[ri].a[ai].b[bi].t.a=aii;
					//printf("mi=%d ri=%d ai=%d bi=%d t.a=%d\n",mi,ri,ai,bi,((*A).m[mi][0].r[ri].a[ai].b[bi].t.a));
					//printf("\nBond: %s.%d.%s.nb=%d,\n",(*A).m[mi][0].r[ri].N,(*A).m[mi][0].r[ri].n,(*A).m[mi][0].r[ri].a[ai].N,(*A).m[mi][0].r[ri].a[ai].nb);
                                        (*A).m[mi][0].r[ri].a[ai].b[bi].t.r=ri;
					(*A).m[mi][0].r[ri].a[ai].b[bi].t.m=mi; 
                                }
                        }
                }
        }
}
return;
}
