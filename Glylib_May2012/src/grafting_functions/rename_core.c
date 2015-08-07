#include <glylib.h> // required for load_pdb function
#include <mylib.h>
#include <molecules.h>

// Renames Core residues originally from ligand.pdb to what they would be called in LG. This way names are consistent (not ROH yet)

//function prototype
void rename_core(assembly *A, assembly *B, int TnumRes, int LGnumRes, int z, int fTtoLG[][40], int *core_ri){
//z is budcore number
printf("************Renaming core*************\n");
//printf("TnumRes=%d\n",TnumRes);
int mi,ri,Ari,Bri;
int a;

for(ri=0;ri<(*A).m[0][0].nr;ri++){ // for each residue in molecule 0 
	for (a=1;a<=(TnumRes+1);a++){ 
		//printf("a=%d,(*A).m[0][0].r[%d].n=%d\n",a,ri,(*A).m[0][0].r[ri].n);
		if ((*A).m[0][0].r[ri].n==a){ // if the residue number is above 1
			printf("Converting core resid %d to %d\na=%d\n",(*A).m[0][0].r[ri].n,fTtoLG[a][z],a);
			(*A).m[0][0].r[ri].n=fTtoLG[a][z]; // convert it to the equivalent core residue number in libary glycan
			Ari=ri;
			core_ri[a]=fTtoLG[a][z];  // added OG 19June2012 change core_ri to reflect changes
			//printf("MARK2\nAri=%d",Ari);

			//Ari=fTtoLG[a][z];
			for(mi=0;mi<(*B).nm;mi++){ // go through library glycan 
				for(Bri=0;Bri<(*B).m[mi][0].nr;Bri++){
					//printf("MARK3\n Aresname=%s,Bresname=%s\n",(*A).m[0][0].r[Ari].N,(*B).m[mi][0].r[Bri].N);
					if (fTtoLG[a][z]==(*B).m[mi][0].r[Bri].n){ // if number matches one that we just changed
						strcpy((*A).m[0][0].r[Ari].N,(*B).m[mi][0].r[Bri].N); // copy LG residue name into template. 
						//printf("MARK4\n");
					}
				}
			}
			ri++; //skip on
		}
	}
}

return;
}
