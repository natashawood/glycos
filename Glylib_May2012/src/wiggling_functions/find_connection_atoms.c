#include <glylib.h>
#include <mylib.h>
#include <molecules.h>
int explore_atom_tree(residue *r, atom *a, dihedral_coord_set *DCS, int *lna, atom *pa, int level);
void reverse_DCS_order(dihedral_coord_set *DCS,int lna);

int find_connection_atoms(int linkage, int linked2, dihedral_coord_set *DCS, assembly *A){
// linkage defines the residue id. We want the co_ord of atoms involved in it's reducing linkage to the next residue.

int mi=0,ri=0,i=0,b_ai=0,c_ai=0;
//int l_mi=0,l_ri=0,l_ai=0;i
//char N1,N2; // type of atom
int lna=0; //linkage number of atoms
int outcome=0;
int linked2_ri=-1;
//*DBG*/printf("linkage=%d, linked2=%d\n",linkage,linked2);

residue *br,*cr;
br=(residue *)calloc(1, sizeof(residue)); //bonded  residue
cr=(residue *)calloc(1, sizeof(residue)); //current residue

atom *bat,*cat;
cat=(atom *)calloc(1, sizeof(atom));
bat=(atom *)calloc(1, sizeof(atom));



// Find the residues involved in the linkage
for(mi=0;mi<(*A).nm;mi++){ // loop through each molecule
	for (ri=0;ri<(*A).m[mi][0].nr;ri++){ // loop through each residue
//		printf("mi=%d,ri=%d,n=%d\n",mi,ri,(*A).m[mi][0].r[ri].n);
		if ((*A).m[mi][0].r[ri].n==linkage){ // if residue is the one we want
//			printf("Found linkage,ri=%d\n",ri);
			br=&(*A).m[mi][0].r[ri];
		}
		if ((*A).m[mi][0].r[ri].n==linked2){
//			printf("Found linked2,ri=%d\n",ri);
			linked2_ri=ri;
                        cr=&(*A).m[mi][0].r[ri];
                }
	}
}

//printf("(*br).nrb=%d\n",(*br).nrb);

// Loop through each residue bond of the branch residue
for (i=0;i<(*br).nrb;i++){ 
	if((*br).rb[i].t.r==linked2_ri){ // if this is the residue bond to linked2
		b_ai=(*br).rb[i].s.a; 
		c_ai=(*br).rb[i].t.a;
		bat=&(*br).a[b_ai];    
        	cat=&(*cr).a[c_ai];
		printf("bat=%s in residue %d\n",(*bat).N,(*br).n);
		printf("cat=%s in residue %d\n",(*cat).N,(*cr).n);
	}
}


// explore in bat direction
//printf("Now for bat\n");
outcome=explore_atom_tree(br,bat,DCS,&lna,bat,0);

//printf("lna=%d\n",lna);
// reverse order of DCS
reverse_DCS_order(DCS,lna);

// explore in cat direction
outcome=explore_atom_tree(cr,cat,DCS,&lna,cat,0);


//free(br);
//free(cr);

return lna;
}

int explore_atom_tree(residue *r, atom *a, dihedral_coord_set *DCS, int *lna, atom *pa, int level){
int i=0;
int ai=0; // atom index
int outcome=0;
int level_lna=*lna;
atom *na;
na=(atom *)calloc(1, sizeof(atom));
printf("Entered explore_atom_tree\n");

for (i=0;i<(*a).nmb;i++){
//	printf("%s has %d neighbours\n",(*a).N,(*a).nmb);
	if((*a).mb[i].t.r==(*a).mb[i].s.r){ // if neighbour in same residue
		ai=(*a).mb[i].t.a;
//		printf("%s(ring=%c) is connected to atom %s(ring=%c)\n",(*a).N,(*a).R,(*r).a[ai].N,(*r).a[ai].R);
		na=&(*r).a[ai];
		// If atom in ring, save any connected ring atom
		if ((*a).R=='Y'){
			if (na->R=='Y'){
				(*lna)=((*lna)+1);
				(*DCS).X[(*lna)+level]=&(*a).x;
//				printf("Dropping out having found %s in resid %d is before end of path, DCS->X[%d]\n",(*a).N,(*a).mb[0].s.r,(*lna+level));
				(*lna)=((*lna)+1);
				(*DCS).X[(*lna)+level]=&(*na).x;
//				printf("Dropping out having found %s in resid %d is end of path, DCS->X[%d]\n",(*na).N,(*na).mb[0].s.r,(*lna+level));
				i=(*a).nmb;// stop looping through connections
				outcome=1;
			}
		}
		else if ((*a).R!='Y'){
			if ( (na->mb[0].s.a!=pa->mb[0].s.a) && (na->mb[0].s.r==pa->mb[0].s.r) ) {// if not previous atom and within same residue
				outcome=explore_atom_tree(r,na,DCS,lna,a,level+1);
				if (outcome==1){
					(*lna)=((*lna)+1);
					(*DCS).X[level_lna+level+1]=&(*a).x;
//					printf("Dropping out found %s in resid %d is in path, DCS->X[%d]\n",(*a).N,(*a).mb[0].s.r,(level_lna+level+1));
					i=(*a).nmb;
				}
			}
		}
		else { printf("Something strange in explore_atom_tree\n");}
	}
}
//free(na);
return outcome;
}

void reverse_DCS_order(dihedral_coord_set *DCS,int lna){
int i;
dihedral_coord_set *tmp;
tmp=(dihedral_coord_set*)calloc(8, sizeof(dihedral_coord_set));
//printf("Entered reverse_DCS_order,lna=%d\n",lna);

for (i=1;i<=lna;i++){	
//	printf("DCS[0].X[%d]=%1f\n",i,(*(DCS->X[i])).i);
	tmp->X[lna-i]=DCS->X[i];
}
//printf("Reordering\n");
for (i=1;i<=lna;i++){
	DCS->X[i]=tmp->X[i-1];
//	printf("DCS[0].X[%d]=%1f\n",i,(*(DCS->X[i])).i);
}
return;
}

