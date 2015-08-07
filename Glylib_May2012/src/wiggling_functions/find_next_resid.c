#include <glylib.h>
#include <mylib.h>
#include <molecules.h>

/* What this very long function does is find the resid of the residue connected to the reducing terminal of current_resid */

int find_next_resid(assembly *A, int current_resid, int *core_ri){
int mi=0,ri=0,ai=0,i=0,j=0,mii=0,rii=0; // indexes
int next_resid=0;
for(mi=0;mi<(*A).nm;mi++){
	for(ri=0;ri<(*A).m[mi][0].nr;ri++){
		if(current_resid==(*A).m[mi][0].r[ri].n){
			for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){
				if ( (strcmp("C1",(*A).m[mi][0].r[ri].a[ai].N)==0) && ((*A).m[mi][0].r[ri].a[ai].nmb>0) ){
					mii=(*A).m[mi][0].r[ri].a[ai].mb[0].t.m;
					rii=(*A).m[mi][0].r[ri].a[ai].mb[0].t.r;
					next_resid=(*A).m[mii][0].r[rii].n;
					printf("Next resid=%d\n",next_resid);
				}
				else if ( (strcmp("C1",(*A).m[mi][0].r[ri].a[ai].N)==0) && ((*A).m[mi][0].r[ri].a[ai].nmb==0) ){ // C1 isn't connected to another residue
					for(i=0;i<(*A).m[mi][0].r[ri].na;i++){
						if ( (strcmp("C2",(*A).m[mi][0].r[ri].a[i].N)==0)){
							mii=(*A).m[mi][0].r[ri].a[i].mb[0].t.m;
	                                   	 	rii=(*A).m[mi][0].r[ri].a[i].mb[0].t.r;
        	                                	next_resid=(*A).m[mii][0].r[rii].n;
                	                        	printf("Next resid=%d\n",next_resid);
						}
					}
				}
			}
		}
	}
}
for (j=1;j<=core_ri[0];j++){
	if (next_resid==core_ri[j]){
		next_resid=0;
		printf("\nOkay hit the core, we are done here\n");
	}
}
return next_resid;
}
