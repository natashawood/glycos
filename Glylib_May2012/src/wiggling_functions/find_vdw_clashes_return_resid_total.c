/*
Oliver Grant 14Mar11
Finds clashes between atoms in assemblies in a pairwize fashion. 
NOTE: xsA and xsB refer to the co-ord sets to use. As this is a clash checker I'm not sure of the usefulness of checking clashes over a traj.
*/
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
double find_vdw_clashes_return_resid_total(assembly *A, assembly *B, int *current_resid, int nbr){ 
int mi=0,mii=0,ri=0,rii=0,ai=0,aii=0;//loop counters for A(i) and B(ii) assemblies
double soap=0.0, total_soap=0.0; // Sphere Overlap Area P? (Can't think of a better one)
double rA,rB; // radii
double x,y,z,d; // xyz are holders to make code more readable. d is distance between two atoms
//FILE * pFile; // Clash results out file
//char filename[500]; // Name of out file
int numAtm=0;
int swtch=0;
//double ovrlap_eqv; //Equivalent overlap. So can say overlap is equivalent to X occluded (carbon) atoms
//double prev_resid_clash=0;
//int clashiest_resid=0;


//compare dist between each atom in A assembly and all atoms in B assembly
//printf("Entered the vdw FUNction\n");

/*printf("Skipping:\n");
for(mi=0;mi<=nbr;mi++){
	for(ri=0;ri<(*A).m[mi][0].nr;ri++){
		printf("%d:%s, ",(*A).m[mi][0].r[ri].n,(*A).m[mi][0].r[ri].N);
	}
}

printf("\nChecking:\n");
for(mi=nbr+1;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
                printf("%d:%s, ",(*A).m[mi][0].r[ri].n,(*A).m[mi][0].r[ri].N);
        }
}
printf("\n");
*/

//for(mi=nbr+1;mi<(*A).nm;mi++){
for(mi=0;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
		if ((*A).m[mi][0].r[ri].n==*current_resid){
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
							// looking back I'm not sure if I confirmed that N[0] would always be the element type...
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
								//printf("Soap=%f,rA=%f,rB=%f,d=%f\n",soap,rA,rB,d); 
									// Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours...
									// NOTE this eqn does account for double overlaps. Each atom is counted against each atom. 
									// NOTE so overlap will be double counted and a higher overlap value will result
									// NOTE the math is availabe in that paper for correcting this but 
									// NOTE I am keeping it simple for now.
								total_soap=total_soap + soap;
							/*	strcpy(filename,fileprefix);
								strcat(filename, "_Clash_results_atoms.txt");
								pFile = fopen (filename,"a+"); // open file and append
								//if (pFile=NULL){ // this is always true... later: pFile==NULL should be. 
								//	printf("File open FAIL!\n");
								//	break ;
								fprintf(pFile,"Clash Area For Resid No.%d Atom_Name:%s is %5.1f \n",(*A).m[mi][0].r[ri].n,(*A).m[mi][0].r[ri].a[ai].N,soap);
								fprintf(pFile," ...Clashing with Resid No.%d Atom_Name:%s in the receptor.\n",(*B).m[mii][0].r[rii].n,(*B).m[mii][0].r[rii].a[aii].N);
								fclose (pFile);
							*/
							}
						}
					}
                        	}
				if (swtch==1){
					numAtm++; // then this atom (should be a sugar atom) ai had some contacts 
				}
				swtch=0; // reset for next atom
			}
			//printf("Total_soap=%f, mi=%d,ri=%d,resid=%d,res=%s\n",total_soap,mi,ri,(*A).m[mi][0].r[ri].n,(*A).m[mi][0].r[ri].N);
			//resid_soap=(total_soap - prev_soap);
			/*if (resid_soap > prev_resid_clash) {
				clashiest_resid=(*A).m[mi][0].r[ri].n;
				printf("Clashiest resid is now %d\n",clashiest_resid);
				prev_resid_clash=resid_soap;
			}
			*/
			//printf("resid_soap=%f\n",resid_soap);
			//prev_soap=total_soap;
        	}	
	}
}
/*
ovrlap_eqv=(total_soap/CSA);
//total_soap=(total_soap / numAtm);
total_soap=(total_soap+0.5); // so can round off by truncation
int total_soap_int;
total_soap_int=total_soap;


if (ovrlap_eqv<=1.0){
	strcpy(filename,fileprefix);
	strcat(filename, "_BINDER_Clash_results_branch.txt");
	pFile = fopen (filename,"a+");
	fprintf(pFile,"\n\tatom_equivalent=%.1f",ovrlap_eqv);
	fclose (pFile);	
}
if (ovrlap_eqv>=1.0){
        strcpy(filename,fileprefix);
        strcat(filename, "_CLASHER_Clash_results_branch.txt");
        pFile = fopen (filename,"a+");
        fprintf(pFile,"\n\tatom_equivalent=%.1f",ovrlap_eqv);
        fclose (pFile);
}
*/
return total_soap; // 
}
