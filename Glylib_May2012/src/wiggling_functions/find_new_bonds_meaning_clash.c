/* OG 12Jun2012
 * Function expects two identical assemblies with different orientations of molecules. 
 * Can send the same assembly twice if you have rotated a bond or whatever.
 * If a new internal bond has formed then the new orientation of the molecule has caused a clash.
 */
#include <mylib.h>
#include <molecules.h>
#include <glylib.h>
#define DIST 1.7
int find_new_bonds_meaning_clash(assembly *A){ 
//printf("in find new bonds function\n");
double x,y,z,d;
int ai=0,aii=0,bi=0;
int is_new_bond=0;
int found_old_bond=0;

int ri,rii;

for(ai=0;ai<(*A).na;ai++){// ai is now current atom to check every other atom against
        for(aii=0;aii<(*A).na;aii++){ //for current atom check dist to every other atom
                x=((*A).a[aii][0].x.i-(*A).a[ai][0].x.i);
                y=((*A).a[aii][0].x.j-(*A).a[ai][0].x.j);
                z=((*A).a[aii][0].x.k-(*A).a[ai][0].x.k);
                d=sqrt((x*x)+(y*y)+(z*z)); //pythagorus square on hyp in 3D
		ri=A[0].a[ai][0].mb[0].s.r;
		rii=A[0].a[aii][0].mb[0].s.r;
		//if ( (A[0].r[ri][0].n==11) && (A[0].r[rii][0].n==8) && (strcmp(A[0].a[ai][0].N,"C5N")==0) && (strcmp(A[0].a[aii][0].N,"CME")==0) ){
		//if ( ri!=rii){
		//printf("Checking ri=%d,ai=%s rii=%d,aii=%s, dist=%.1f\n",A[0].r[ri][0].n,A[0].a[ai][0].N,A[0].r[rii][0].n,A[0].a[aii][0].N,d);
		//}
                if ( ri!=rii && d<=DIST && d>0.001){ // If not the same atom but close enough for a bond
		//	printf("Checking r=%d,a=%s r=%d,a=%s, dist=%.1f\n",A[0].r[ri][0].n,A[0].a[ai][0].N,A[0].r[rii][0].n,A[0].a[aii][0].N,d);
			//printf("Checking r=%d,a=%s r=%d,a=%s\n",A[0].r[ri][0].n,A[0].a[ai][0].N,A[0].r[rii][0].n,A[0].a[aii][0].N);
			for (bi=0;bi<A[0].a[ai][0].nmb;bi++){ // Check all mb in B to see if it's a new bond or not.
			//	printf("ai=%d,a=%d,N=%s,N=%s\n",ai,A[0].a[ai][0].mb[0].s.a,A[0].a[ai][0].N);
			//	printf("A[0].a[ai][0].mb[bi].t.a=%d,aii=%d\n",A[0].a[ai][0].mb[bi].t.a,A[0].a[aii][0].mb[0].s.a);
				if( (A[0].a[ai][0].mb[bi].t.a==A[0].a[aii][0].mb[0].s.a) && (A[0].a[ai][0].mb[bi].t.r==A[0].a[aii][0].mb[0].s.r) ){
					found_old_bond=1;
					bi=A[0].a[ai][0].nmb; // finish up faster if we've found it
		//			if ( (A[0].r[ri][0].n==11) && (A[0].r[rii][0].n==8) && (strcmp(A[0].a[ai][0].N,"C5N")==0) && (strcmp(A[0].a[aii][0].N,"CME")==0) ){
		//				printf("Found old bond\n");
		//			}
				}
			}
			if (found_old_bond==0){
				is_new_bond=1;
		//		printf("Found new internal bond!\n");
				ai=(*A).na;  // If find even one new bond quit
				aii=(*A).na; // If find even one new bond quit
			}
			found_old_bond=0;// reset
		}
	}
}

return is_new_bond;
}
