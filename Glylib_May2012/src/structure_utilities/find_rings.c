/** \file find_rings.c
 *  \brief Functions for finding rings.
 *
 *  begin 20110421 BLFoley 
 */

/**  set_smallest_rings_from_residue_atom_nodes finds rings in a residue.
 *
 *  This function might not work properly unless the bonding trees were 
 *  set according to the internal "from bonds" functions.  However, even 
 *  if the function thinks this might be the case, it tries anyway.
 *
 *  If there is position information in either the main coordinate or 
 *  the first alternate coordinate set, the function will set the total
 *  length (s) based on that information.  If position information is 
 *  not available, it will be set to zero.
 */

/* Edited by Oliver Grant. As of 28May2012 it can handle carbohydrates
 * Residues with multiple or fused ring systems will fail
 * Some code is redundant/overly complicated as we each had different ideas
 * The how-to was worked out as I coded, so the codes messy.
 */
#include <stdlib.h>
#include <molecules.h>
#include <mylib.h>
void walk_along_paths(int **path, int pathn, int ci, residue *r);


void set_smallest_rings_from_residue_atom_nodes(residue *r)
{
int i=0,nri=0,ai=0,nrings=0;
int nr_localmax=0,*nstarts,nrstarts=0,niimax=0,npath=0;
//int ni1=0,ni2=0;
char isorigin='n';

//molindex_set *RING;
/*
struct
    {
    int *steps;
    ensindex *parent;
    } path_follower;
path_follower *pf;

struct
	{
	int *steps;
	int *child;
	}
ring_follower *rf;
*/
/*
    0.  Check that the nodes are set and seem sane.
        If the origin has one or more incoming, there 
        might be a problem.
*/

//printf("\nResidue is %d\n",r[0].n);
if(r[0].aT==NULL){mywhine("r[0].aT==NULL in set_smallest_rings_from_residue_atom_nodes.");}
if(r[0].aT[0].isorigin!='Y'){mywhine("r[0].aT[0] is not origin in set_smallest_rings_from_residue_atom_nodes.");}

if(r[0].na<3)
    {  /* If there aren't enough atoms to make a ring. */
    r[0].nring=0;
    return;
    }

isorigin='n';
for(ai=0;ai<r[0].na;ai++)
    {
    if(r[0].aT[ai].ni<0){mywhine("Unexpected value of r[0].aT[ai].ni in set_smallest_rings_from_residue_atom_nodes.");}
    if(r[0].aT[ai].no<0){mywhine("Unexpected value of r[0].aT[ai].no in set_smallest_rings_from_residue_atom_nodes.");}
    if(r[0].aT[ai].isorigin=='Y') isorigin = 'y';
    }
if(isorigin!='y')
	{
	printf("\nWARNING: No atom in tree set as origin in \n\tset_smallest_rings_from_residue_atom_nodes.\n");
	printf("This will almost certainly cause trouble.\n");
	printf("Setting first atom as origin and trying anyway.\n");
	printf("If other things fall apart, this is a possible reason.\n");
	}

/*
        If there are incoming bonds to the first atom, the atom nodes probably 
        aren't set properly for this function.  But, give it a go anyhow.  
*/
if(r[0].aT[0].ni>0)
    {
    printf("\nWARNING: r[0].aT[0].ni>0 in set_smallest_rings_from_residue_atom_nodes.\n");
    printf("So,the atom_node bonding probably isn't right for this function.\n");
    printf("Ignoring that and forging on anyway.\n");
    }
if(r[0].aT[0].isorigin!='Y')
    {
    printf("\nWARNING: r[0].aT[0] is not origin in set_smallest_rings_from_residue_atom_nodes.\n");
    printf("So,the atom_node bonding probably isn't right for this function.\n");
    printf("Ignoring that and forging on anyway.\n");
    }

/*
        Check number of rings.  Record start locations.
        Record lots of other bits of information.
*/

/***************Find atoms with more than one incoming bond*****************/
//printf("r[0].aT[0].isorigin=%c\n",r[0].aT[0].isorigin);
nrings=0;
nstarts=(int*)calloc(r[0].na,sizeof(int));
int a=0;
for(ai=0;ai<r[0].na;ai++) 
    { 
//    if(r[0].aT[ai].isorigin=='Y') nstarts[0]= ( r[0].aT[0].ni+1 ) * ( r[0].aT[0].ni ) / 2;
//    else nstarts[ai] = ( r[0].aT[ai].ni-1 ) * ( r[0].aT[ai].ni-2 ) / 2;
    a=r[0].aT[ai].ID.a; 
//    printf("atom %s in residue %s has %d incoming\n",r[0].a[a].N, r[0].N, r[0].aT[ai].ni);
    nstarts[ai] = ( r[0].aT[ai].ni-1 ); // Oliver
    if(nstarts[ai]>0)
        {
        nrstarts++;
        nrings += nstarts[ai];
        if(nstarts[ai]>nr_localmax) nr_localmax=nstarts[ai];
        if(r[0].aT[ai].ni>niimax) niimax=r[0].aT[ai].ni;
        }
    }
int start[nrstarts];
//start=(int*)calloc(nrstarts,sizeof(int));
nri=0;
for(ai=0;ai<r[0].na;ai++)
    { 
    if(nstarts[ai]>0)
        {
	start[npath]=ai;
	npath++; // not doing anything right now
	nri+=nstarts[ai]; /* a quick code check */ 
	}
    }
if(nri!=nrings){mywhine("nri!=nrings in set_smallest_rings_from_residue_atom_nodes.");}

/*******************For now quit if there will be more than one ring per residue********************/
if(nrstarts>1){mywhine("Can't handle more than one ring per residue in set_smallest_rings_from_residue_atom_nodes.");}

/***************Lachele's data structure idea for posterity**************/
/*
        Assign accounting spaces.

*/

/*rf = (ring_follower*) calloc(r[0].na, sizeof(ring_follower));
for(ai=0;ai<r[0].na;ai++)
    { 
    rf[ai].steps=(int*)calloc(niimax,sizeof(int));
    rf[ai].child=(int*)calloc(niimax,sizeof(int));
    }
*/ // Lachele's original code. I'm keeping it here just in case I need to see the idea again.

/*
// I'm doing this different to Lachele's original idea for a rf
rf = (ring_follower*) calloc(nstarts, sizeof(ring_follower));
for (ai=0;ai<r[0].na;ai++) // could just be a ring of atoms as molecule. So max steps would be all the atoms.
    {
    rf[i].steps=(int*)calloc(r[0].na,sizeof(int)); 
    rf[i].child=(int*)calloc(npath,sizeof(int)); // number of child paths would be equal to number of incoming bonds beyond (total number of atoms-1) -1 as origin will not have incoming.
    }
//RING = (molindex_set*) calloc(nrings, sizeof(molindex_set));
*/

/***************Create space to record paths************************/
int npaths=5; // Do better than this Oly, dynamically expand in function

int **path=(int **)malloc(npaths*sizeof(int *));
      for(i=0;i<npaths;i++)path[i]=(int *)malloc((r[0].na)*sizeof(int));
int j=0,k=0;

for (j=0;j<npaths;j++){
	path[j][0]=0;
	for (k=1;k<r[0].na;k++){
		path[j][k]=-1;
	//	printf("path[%d][%d]=%d\n",j,k,path[j][k]);
	}
}


//printf("About to enter path walker\nr[0].na=%d,nrstarts=%d,start[0]=%d\n\n",r[0].na,nrstarts,start[0]);

/*
for (i=0;i<r[0].na;i++){
	printf("r[0].aT[%d].isorigin=%c\n",i,r[0].aT[i].isorigin);
}
*/

/**************Find paths from start points*******************/
for (i=0;i<nrstarts;i++){
    	path[i][0]=0; //reset each loop
    	walk_along_paths(path, i, (start[i]), r);
}
/*
// print out the paths that were found
for (i=0;i<npaths;i++){
	for (j=1;j<=path[i][0];j++){
		printf("path[%d][%d]=%d\n",i,j,path[i][j]);
	}
	printf("\n");
}
*/
/**********************Sort paths into rings********************/
// Just doing it for simple carbohydrates for now
int current=0;
//int ring[r[0].na];
//ring[0]=0;
int end1=path[0][path[0][0]];
int end2=path[1][path[1][0]];

// Trim the ends off the paths until the numbers are unique
while (end1==end2){
	path[0][0]=((path[0][0])-1);
        path[1][0]=((path[1][0])-1);
	end1=path[0][path[0][0]];
	end2=path[1][path[1][0]];
	//printf("end1=%d,end2=%d\n",end1,end2);
} // once that's over

// include one atom before that matched
path[1][0]=((path[1][0])+1);

for (i=0;i<=1;i++){
	for (j=1;j<=path[i][0];j++){		
		current=path[i][j];
		r[0].a[current].R='Y';
		//printf("%s is part of ring\n",r[0].a[current].N);
	}
}
free(nstarts);
for(i=0;i<npaths;i++){free(path[i]);}
free(path);
return;
}

void walk_along_paths(int **path, int pathn, int ci, residue *r){ // current atom general index
int ai=0,i=0;
int stop=0;
int ns=0;
int new_ci=0;
//printf("Have just entered walk_along_paths\n");

// Outcome1: Hit origin?
if (r[0].aT[ci].isorigin=='Y') {
	//printf("Hit the origin\n");
	ns=path[pathn][0]=((path[pathn][0])+1);
        path[pathn][ns]=ci;
	//printf("path[%d][%d]=%d\n",pathn,ns,path[pathn][ns]);
	stop=1;
}

// Outcome2: Hit start?
//if (ci==start) (stop=1);
// Outcome3: Hit something already in path or another path?
if (stop!=1){
	for (ai=1;ai<=path[pathn][0];ai++){
    	//	printf("path[%d][%d]=%d,ci=%d\n",pathn,ai,path[pathn][ai],ci);
    		if (ci==path[pathn][ai]) {
			stop=1;
	//		printf("Have hit something already in path\n");
		}
	}
}
    
if (stop!=1){
    	for (i=0;i<r[0].aT[ci].ni;i++){ // for each incoming
		// add this residue to path
		pathn=pathn+i; // add a child if there is more than 1 incoming
		ns=path[pathn][0]=((path[pathn][0])+1);
        	path[pathn][ns]=ci;
		//printf("ci=%d,ni=%d\n",ci,r[0].aT[ci].ni);
       	 	new_ci=r[0].aT[ci].i[i][0].a; 
		//printf("Will now walk along path from new_ci=%d\n",new_ci);
        	walk_along_paths(path, pathn, new_ci, r);
        }
}
return;
}
