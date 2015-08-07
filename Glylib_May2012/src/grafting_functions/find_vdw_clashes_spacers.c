/*
Oliver Grant 10Jan13
Finds clashes between atoms in assemblies in a pairwize fashion. 
*/
#include <glylib.h>
#include <mylib.h> 
#include <molecules.h>
// Bondi definitions of VDW in Angstrom HERE AS A REFERENCE, SET IN CODE BELOW
//#define H 	1.09 // Rowland and Taylor modification
//#define C	1.70
//#define N 	1.55
//#define O	1.52
//#define F 	1.47
//#define P	1.80
//#define S	1.80
//#define Cl	1.75

// Carbon Surface Area for normalizing results to the area of a carbon atom. Can ask the question: How many C atom equivalents are inside the protein? 
#define CSA 36.31681103 

// nbr is number of molecules in ligand.pdb. So skip checking first nbr molecules in A
int find_vdw_clashes_spacers(assembly *A, assembly *B, char *spacername, char *spaceroutname, int *fitterCnt){ // A is Spacer, B is receptor 
int mi=0,mii=0,ri=0,rii=0,ai=0,aii=0;//loop counters for A(i) and B(ii) assemblies
double soap, total_soap=0.0; // Sphere Overlap Area P? (Can't think of a better one)
double rA,rB; // radii
double x,y,z,d; // xyz are holders to make code more readable. d is distance between two atoms
FILE * pFile; // Clash results out file
char filename[500]; // Name of out file
int numAtm=0;
int swtch=0;
double ovrlap_eqv; //Equivalent overlap. So can say overlap is equivalent to X occluded atoms
int recClashRi[500]; // 500 is arbitrary. [0] will be length
recClashRi[0]=0;
int i=0;
int save=1;
char res[10]; // holds residue name and number as string

//compare dist between each atom in A assembly and all atoms in B assembly
printf("Entered the vdw FUNction for spacers\n");
printf("spacername=%s\nspaceroutname=%s\n",spacername,spaceroutname);

strcpy(filename,spaceroutname);
//strcpy(filename,"Results.txt");

for(mi=0;mi<(*A).nm;mi++){
    for(ri=0;ri<(*A).m[mi][0].nr;ri++){
        if( (strcmp((*A).m[mi][0].r[ri].N,"0YB")!=0) && (strcmp((*A).m[mi][0].r[ri].N,"0VA")!=0) ){// Don't check sugar part         
            for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){// ai is now current atom to check against every atom in assembly B
                for(mii=0;mii<(*B).nm;mii++){
                    for(rii=0;rii<(*B).m[mii][0].nr;rii++){
                        for(aii=0;aii<(*B).m[mii][0].r[rii].na;aii++){//for current atom check dist to every other atom within the assembly B
                            //printf("Comparing Values\n rA=%f\n",rA);
                            x=((*A).m[mi][0].r[ri].a[ai].x.i-(*B).m[mii][0].r[rii].a[aii].x.i);
                            y=((*A).m[mi][0].r[ri].a[ai].x.j-(*B).m[mii][0].r[rii].a[aii].x.j);
                            z=((*A).m[mi][0].r[ri].a[ai].x.k-(*B).m[mii][0].r[rii].a[aii].x.k);
                            d=sqrt((x*x)+(y*y)+(z*z)); //pythagorus square on hyp in 3D
                            //printf("\ndist is %f\n",d);
			    // set radii of the atoms
			    // looking back I'm not sure if I confirmed that N[0] would always be the element type... seems to be
			    if ((*A).m[mi][0].r[ri].a[ai].N[0]=='C'){rA=1.70;} 
			    if ((*A).m[mi][0].r[ri].a[ai].N[0]=='O'){rA=1.52;} 
			    if ((*A).m[mi][0].r[ri].a[ai].N[0]=='N'){rA=1.55;} 
			    if ((*A).m[mi][0].r[ri].a[ai].N[0]=='S'){rA=1.80;} 

			    if ((*B).m[mii][0].r[rii].a[aii].N[0]=='C'){rB=1.70;}
                            if ((*B).m[mii][0].r[rii].a[aii].N[0]=='O'){rB=1.52;}
			    if ((*B).m[mii][0].r[rii].a[aii].N[0]=='N'){rB=1.55;}
			    if ((*B).m[mii][0].r[rii].a[aii].N[0]=='S'){rB=1.80;}
						
			    //printf("Dist is %f, rA is %f, rB is %f \n\n",d,rA,rB);
			    // if the sum of the radii is greater than the distance between them.
			    if (rA + rB > d + 0.6){ // 0.6 overlap is deemed acceptable. (Copying chimera:) 
			        swtch=1;
  				soap=2*PI*rA*(rA-d/2-(((rA*rA)-(rB*rB))/(2*d))); 
				// Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours...
				// NOTE this eqn does account for double overlaps. Each atom is counted against each atom. 
				// NOTE so overlap will be double counted and a higher overlap value will result
				// NOTE the math is availabe in that paper for correcting this but 
				// NOTE I want to switch to a LJ type eqn anyway.
				// NOTE Do I? No, I don't. It's an estimate.
				total_soap=total_soap + soap;
				save=1;
				for(i=1;i<=recClashRi[0];i++){
				    if(rii==recClashRi[i]){save=0;}
				}
				if(save==1){
				    recClashRi[0]=(recClashRi[0]+1); //increment
				    recClashRi[recClashRi[0]]=rii;
				    printf("Saving %d\n",rii);
                                }
                            }
                        }
                    }
                }
                if (swtch==1){
                    numAtm++; // then this atom (should be a sugar atom) ai had some contacts 
                }
		swtch=0; // reset for next atom
            }
	//resid_soap=(total_soap - prev_soap);
	//strcpy(filename,fileprefix);
	//strcat(filename, "_Spacer_Clash_results_resid.txt");
	//pFile = fopen (filename,"a+");
	//fprintf(pFile,"Clash Area For Resid No.%d is %5.1f\n",(*A).m[mi][0].r[ri].n,resid_soap);
	//fclose (pFile);
	//prev_soap=total_soap;
        }
    }
}
ovrlap_eqv=(total_soap/CSA);
//total_soap=(total_soap / numAtm);
total_soap=(total_soap+0.5); // so can round off by truncation

int res_below=0;
res_below=check_surface_plane(A,B);

if ((ovrlap_eqv<=1.0) && (res_below==0)) {
	*fitterCnt=((*fitterCnt)+1);
}

pFile = fopen (filename,"a+");
fprintf(pFile,"%c%c%c%c%c%c|%6.1f | %3d | ",spacername[3],spacername[4],spacername[5],spacername[6],spacername[7],spacername[8],ovrlap_eqv,res_below);
for(i=1;i<=recClashRi[0];i++){ 
	sprintf(res,"%s%3d",(*B).m[0][0].r[recClashRi[i]].N,(*B).m[0][0].r[recClashRi[i]].n);
	fprintf(pFile,"%6s, ",res);
	//fprintf(pFile,"%3s%3d,",(*B).m[0][0].r[recClashRi[i]].N,(*B).m[0][0].r[recClashRi[i]].n);  //Shit, assuming mii is 0.
}
fprintf(pFile,"\n");
fclose (pFile);

/* Let's not wiggle the Spacers for now
// Need some distance from core based metric
if (ovrlap_eqv>=1.0){
	if (wiggle=='Y') {
		total_soap=wiggler(A,B,core_ri);
		ovrlap_eqv=(total_soap/CSA);
	
		if (ovrlap_eqv<=1.0){
			strcpy(filename,fileprefix);
        		strcat(filename, "_BINDER_Clash_results_branch.txt");
			pFile = fopen (filename,"a+");
        		fprintf(pFile,"\n\tPost wiggling atom_equivalent=%.1f residues below plane=%d",ovrlap_eqv,res_below);
        		fclose (pFile);
		}
	}
	if (ovrlap_eqv>=1.0){
        	strcpy(filename,fileprefix);
        	strcat(filename, "_CLASHER_Clash_results_branch.txt");
        	pFile = fopen (filename,"a+");
        	fprintf(pFile,"\n\tPost wiggling atom_equivalent=%.1f",ovrlap_eqv);
        	fclose (pFile);
	}
}
*/

return total_soap;
}
