/* OG 12Jun2012
 * Function expects two identical assemblies with different orientations of molecules. 
 * Can send the same assembly twice if you have rotated a bond or whatever.
 * If a new internal bond has formed then the new orientation of the molecule has caused a clash.
 */
#include <mylib.h>
#include <molecules.h>
#include <glylib.h>
#define DIST 1.7
int find_bonds_meaning_clash(assembly *A, assembly *B){ 
//printf("in find new bonds function\n");
double x,y,z,d;
int ai=0,aii=0;
int Anmb=0,Bnmb=0; // number of molecule bonds in A and B

int ri,rii;
int is_new_bond=0;

//printf("Assembly A:\n");
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
			Anmb++;
			//printf("%s:%s-%s:%s\n",(*A).m[0][0].r[ri].N,(*A).a[ai][0].N,(*A).m[0][0].r[rii].N,(*A).a[aii][0].N);
		}
	}
}
//printf("Assembly B:\n");
for(ai=0;ai<(*B).na;ai++){// ai is now current atom to check every other atom against
        for(aii=0;aii<(*B).na;aii++){ //for current atom check dist to every other atom
                x=((*B).a[aii][0].x.i-(*B).a[ai][0].x.i);
                y=((*B).a[aii][0].x.j-(*B).a[ai][0].x.j);
                z=((*B).a[aii][0].x.k-(*B).a[ai][0].x.k);
                d=sqrt((x*x)+(y*y)+(z*z)); //pythagorus square on hyp in 3D
                ri=B[0].a[ai][0].mb[0].s.r;
                rii=B[0].a[aii][0].mb[0].s.r;
		//printf("ai=%d,aii=%d,ri=%d,rii=%d,dist=%.2f\n",ai,aii,ri,rii,d);
                if ( ri!=rii && d<=DIST && d>0.001){ // If not the same atom but close enough for a bond
                        Bnmb++;
			//printf("%s:%s-%s:%s\n",(*B).m[0][0].r[ri].N,(*B).a[ai][0].N,(*B).m[0][0].r[rii].N,(*B).a[aii][0].N);
                }
        }
}

if (Anmb!=Bnmb){
	printf("Anmb=%d,Bnmb=%d\n",Anmb,Bnmb);
	printf("*******************************\nAre there TER cards in your library glycans? This isn't leap! Remove them!\n\n\n");
	is_new_bond=1;
}

return is_new_bond;
}
