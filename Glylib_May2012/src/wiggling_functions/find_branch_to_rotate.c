#include <glylib.h>
#include <mylib.h>
#include <molecules.h>

/* Based on linkage residue this function will find which linkage leads to the core.*/
// function declarations:
//void follow_branch_from_residue(int core_dir, int c_ri, int c_mi, int crn, assembly *A, coord_3D **brnch, int *nab);
void follow_atom_tree(assembly *A, int mi, int ri, int ai, coord_3D **brnch, int *nab);
char check_if_visited(residue *r, int ai, coord_3D **brnch, int *nab);
void find_ring_atom(assembly *A, dihedral_coord_set *DCS, int lna, atom **Pring_bat, atom **Pprev_at);
void add_surrounding_atoms(assembly *A, int i, dihedral_coord_set *DCS, int lna, coord_3D **brnch, int *nab);

// START HERE:
void find_branch_to_rotate(assembly *A, int tors, coord_3D **brnch, int *nab, dihedral_coord_set *DCS, int lna, int linked2){

/* get atom (Aa) in linked_resid that is connected to atom (Ab) in core_dir */
/* Go to atom connected to ring atom in core_dir. From here travel to all it's neighbours bar where we came from saving co_cord. Move to next atom. Repeat.ti*/

int c_mi=0,c_ri=0,c_ai=0; // Connected residue index, connected atom index
int i=0;
char is_found='N';

atom **Pprev_at,**Pring_bat;
atom *ring_bat,*prev_at;
//Pring_bat=(atom **)calloc(1, sizeof(atom*));
//Pprev_at=(atom **)calloc(1, sizeof(atom*));
//ring_bat=(atom *)calloc(1, sizeof(atom));
//prev_at=(atom *)calloc(1, sizeof(atom));
Pring_bat=&ring_bat;
Pprev_at=&prev_at;

/************Find ring atom in branch which is in the linkage to core ****************/
find_ring_atom(A,DCS,lna,Pring_bat,Pprev_at);

/************************Include the ring atom in rotation****************************/
//	printf("Saving coord of atom %s which is ring_bat\n",ring_bat->N);
	brnch[*nab]=&(*ring_bat).x;
	*nab=(*nab+1);

/********find all atoms in branch residue and beyond we need to move******************/
for (i=0;i<(*ring_bat).nmb;i++){ // for each atom connected to ring atom
	c_mi=(*ring_bat).mb[i].t.m;
	c_ri=(*ring_bat).mb[i].t.r;
	c_ai=(*ring_bat).mb[i].t.a;
	//printf("c_ai=%d\n",c_ai);
	is_found=check_if_visited(&(*A).m[c_mi][0].r[c_ri], c_ai, brnch, nab); // have we already found it in our travels?
	if ( ((*A).m[c_mi][0].r[c_ri].a[c_ai].n!=prev_at->n) && (is_found=='N' ) ){ // Make sure we are not going down the linkage
	//	printf("%d is not equal to %d\n",(*A).m[c_mi][0].r[c_ri].a[c_ai].n,prev_at->n);
		follow_atom_tree(A, c_mi, c_ri, c_ai, brnch, nab); // recursively finds all the atoms to be rotated
	}
}

/**************add atoms in linkage that we need to move******************/
// For 3rd bond we need to move 2nd atom in DCS, for 4th bond we need to move 2nd and 3rd. 
for (i=tors;i>2;i--) { 
//	printf("As tors=%d,saving coord of atom which is in DCS->[%d]\n",tors,i);
	brnch[*nab]=DCS->X[i];
	//printf("brnch[%d]=%p,DCS->X[%d]=%p\n",(*nab),brnch[*nab],tors,DCS->X[tors]);
	//printf("x=%f in brnch[%d]\n",(*brnch[*nab]).i,*nab);
	*nab=(*nab+1);
	add_surrounding_atoms(A,i,DCS,lna,brnch,nab);
}


//printf("number of atoms in branch=%d\n",*nab);

//free(Pring_bat);
//free(Pprev_at);
//free(ring_bat);
//free(prev_at);


return;
}

void add_surrounding_atoms(assembly *A, int i, dihedral_coord_set *DCS, int lna, coord_3D **brnch, int *nab){
// For 2-8 linkages with atoms hanging off torsion atoms. Need to rotate.
int ai=0,ri=0,mi=0;
int r=0,m=0;
double x,y,z;
char atom[4];
// Find neighbours of current atom we are adding to brnch
for(mi=0;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
                for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){
                        x=(*A).m[mi][0].r[ri].a[ai].x.i;
                        y=(*A).m[mi][0].r[ri].a[ai].x.j;
                        z=(*A).m[mi][0].r[ri].a[ai].x.k;
                        if( (x==(*(DCS->X[i])).i) && (y==(*(DCS->X[i])).j) && (z==(*(DCS->X[i])).k) ){
				r=ri;
				m=mi;
				strcpy(atom,(*A).m[mi][0].r[ri].a[ai].N);
				// finish looping
//				ai=(*A).m[mi][0].r[ri].na;
//				ri=(*A).m[mi][0].nr;
//				mi=(*A).nm;
                        }
                }
        }
}
// Let's just do this ugly and only for 2-8 linkages.
// When the atom is C8, add C9 and O9. When it's C7 add O7.
if (strcmp(atom,"C8")==0){
	for(ai=0;ai<(*A).m[m][0].r[r].na;ai++){
		if ( (strcmp((*A).m[m][0].r[r].a[ai].N,"C9")==0) || (strcmp((*A).m[m][0].r[r].a[ai].N,"O9")==0) ){
			brnch[*nab]=&(*A).m[m][0].r[r].a[ai].x;
                	*nab=(*nab+1);
		}
	}
}
if (strcmp(atom,"C7")==0){
	for(ai=0;ai<(*A).m[m][0].r[r].na;ai++){
		if (strcmp((*A).m[m][0].r[r].a[ai].N,"O7")==0){
			brnch[*nab]=&(*A).m[m][0].r[r].a[ai].x;
                        *nab=(*nab+1);
                }
        }
}
	
return;
}



void find_ring_atom(assembly *A, dihedral_coord_set *DCS, int lna, atom **Pring_bat, atom **Pprev_at){
int ai=0,ri=0,mi=0;
double x,y,z;
for(mi=0;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
                for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){
			x=(*A).m[mi][0].r[ri].a[ai].x.i;
			y=(*A).m[mi][0].r[ri].a[ai].x.j;
			z=(*A).m[mi][0].r[ri].a[ai].x.k;
			if( (x==(*(DCS->X[2])).i) && (y==(*(DCS->X[2])).j) && (z==(*(DCS->X[2])).k) ){
				(*Pring_bat)=&(*A).m[mi][0].r[ri].a[ai];
			//	printf("ring_bat.N=%s,ai is %d\n",(*Pring_bat)->N,(*Pring_bat)->mb[0].s.a);
			}
			if( (x==(*(DCS->X[3])).i) && (y==(*(DCS->X[3])).j) && (z==(*(DCS->X[3])).k) ){
                                (*Pprev_at)=&(*A).m[mi][0].r[ri].a[ai];
			//	printf("prev_at.N=%s,ai is %d\n",(*Pprev_at)->N,(*Pprev_at)->mb[0].s.a);
                        }
		}
	}
}
return;
}

void follow_atom_tree(assembly *A, int mi, int ri, int ai, coord_3D **brnch, int *nab){
int i=0,nmi=0,nri=0,nai=0;
char is_found='N';
//printf("ai=%d\n",ai);
//printf("in r x=%f,y=%f,z=%f\n",r->a[ai].x.i,r->a[ai].x.j,r->a[ai].x.k);
//printf("r->a[%d].x=%p\n",ai,&r->a[ai].x);
//printf("Saving coord of atom %s in residue.n=%d\n",(*A).m[mi][0].r[ri].a[ai].N,(*A).m[mi][0].r[ri].n);
brnch[(*nab)]=&((*A).m[mi][0].r[ri].a[ai].x); // save coord
//printf("brnch[%d]=%p\n",*nab,brnch[*nab]);
//printf("r.n=%d,a.n=%d\n",r->n,r->a[ai].n);
//printf("nab=%d,x=%f,y=%f,z=%f\n",*nab,(*brnch[*nab]).i,(*brnch[*nab]).j,(*brnch[*nab]).k);
*nab=(*nab+1);
//printf("which has %d mb\n",(*A).m[mi][0].r[ri].a[ai].nmb);
for (i=0;i<(*A).m[mi][0].r[ri].a[ai].nmb;i++){ // for each atom connected to this atom
        nai=(*A).m[mi][0].r[ri].a[ai].mb[i].t.a;
	nri=(*A).m[mi][0].r[ri].a[ai].mb[i].t.r;
	nmi=(*A).m[mi][0].r[ri].a[ai].mb[i].t.m;
	//printf("%d in resi %d is bonded to %d in resi %d\n",ai,ri,nai,nri);
	is_found=check_if_visited(&(*A).m[nmi][0].r[nri], nai, brnch, nab);
	if (is_found=='N'){follow_atom_tree(A, nmi, nri, nai, brnch, nab);}
	else if (is_found=='Y'){ /*printf("Already visited this one\n");*/}
}

return;
}

char check_if_visited(residue *r, int ai, coord_3D **brnch, int *nab){
// checks if have been at this atom before and saved it
// would have been simpler to save and check resnames but this reduces variables used in code above
// could have used a system within the structure to indicate visited too. cleaner code.
int i=0;
double x,y,z;
x=r->a[ai].x.i;
y=r->a[ai].x.j;
z=r->a[ai].x.k;
char is_found='N';

for (i=0;i<(*nab);i++){
//	printf("x=%f and found x=%f\n",x,(*brnch[i]).i);
	if ( ((*brnch[i]).i==x) && ((*brnch[i]).j==y) && ((*brnch[i]).k==z) ){
		is_found='Y';
	}
}
//printf("is_found=%c\n",is_found);
return is_found;
}	
