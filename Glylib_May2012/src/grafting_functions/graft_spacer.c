// Function written by Oliver Grant August 2012
#include <glylib.h>
#include <mylib.h>
#include <molecules.h>
#include <dirent.h>
int graft_spacer(assembly *A, assembly *B, char *fileprefix){ // A is carb. B is receptor

assembly Sp;
char filename[100];
char cplx_a1[3], cplx_a2[3], cplx_a3[3], cplx_a4[3];
int c_ai=0,c_ri=0,c_mi=0; // core ai, ri and mi
int i=0, j=0, a4_ai=0, a3_ai=0, a1_ri=0, a1_mi=0, prev=0;
int a3=0,a4=0; // for sorting out which atom is a3 and which is a4 in each scenario
int ai=0,ri=0,mi=0;
int bbai=0,bri=0,bai=0,bmi=0;
char ROHatom[3];
int fitterCnt=0;

// results output
char spaceroutname[300];
sprintf(spaceroutname,"%s_spacer.out",fileprefix);
printf("spaceroutname=%s\n",spaceroutname);

// pdb output file
char pdboutname[300];
int count=10000; // set to 10000 so vmd reads sequentially

set_nbonds_for_atoms_in_assembly(A);

atom *A_a1, *A_a2, *A_a3, *A_a4, *Sp_a1, *Sp_a2, *Sp_a3, *Sp_a4; // pointers to atoms
A_a1=(atom *)calloc(1, sizeof(atom));
A_a2=(atom *)calloc(1, sizeof(atom));
A_a3=(atom *)calloc(1, sizeof(atom));
A_a4=(atom *)calloc(1, sizeof(atom));
Sp_a1=(atom *)calloc(1, sizeof(atom));
Sp_a2=(atom *)calloc(1, sizeof(atom));
Sp_a3=(atom *)calloc(1, sizeof(atom));
Sp_a4=(atom *)calloc(1, sizeof(atom));
//coord_3D *pt1,*pt2,*pt3;
//pt1=(coord_3D*)calloc(1, sizeof(coord_3D));
//pt2=(coord_3D*)calloc(1, sizeof(coord_3D));
//pt3=(coord_3D*)calloc(1, sizeof(coord_3D));
char ano='A'; // Anomeric type for linker. Gets set below.

//In carb: rename atoms so they can be superimposed. Find anomeric type for linker.
for(mi=0;mi<(*A).nm;mi++){ // for each molecule
    for(ri=0;ri<(*A).m[mi][0].nr;ri++){ // for each residue in molecule
        if( (strcmp((*A).m[mi][0].r[ri].N,"ROH")==0) || (strcmp((*A).m[mi][0].r[ri].N,"OME")==0) ){ // if ROH or OME
            printf("Found ROH atom in A\n");
            for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){ // for each atom
                if( (strcmp((*A).m[mi][0].r[ri].a[ai].N,"O1")==0) || (strcmp((*A).m[mi][0].r[ri].a[ai].N,"O")==0) ){ // If O1 or O
                    strcpy(ROHatom,(*A).m[mi][0].r[ri].a[ai].N);
                    A_a2=&(*A).m[mi][0].r[ri].a[ai]; // pointer to a2 atom
                    // find atom (C1) in other residue that this one is connected to and the connect two atoms (O5,C2)
                }
            }
        }
    }
}
for(i=0;i<A_a2->nmb;i++){
    if(A_a2->mb[i].t.r!=A_a2->mb[i].s.r){
                c_ai=A_a2->mb[0].t.a;
        c_ri=A_a2->mb[0].t.r;
        c_mi=A_a2->mb[0].t.m;
        A_a1=&(*A).m[c_mi][0].r[c_ri].a[c_ai];
        printf("residue connected to ROH is %s\n",(*A).m[c_mi][0].r[c_ri].N);
	ano=(*A).m[c_mi][0].r[c_ri].N[2];
        //printf("Anomeric configuration is %c\n",ano);
        strcpy(cplx_a1,A_a1->N);
        strcpy(cplx_a2,A_a2->N);
        strcpy(A_a1->N,"a1");
        strcpy(A_a2->N,"a2");
        for (i=0;i<A_a1->nb;i++){ // for each atom bonded to a1
            bbai=A_a1->b[i].t.a; // readability: save ai of bonded atom
            printf("A a1 number bonds=%d, bbai=%d, c_ai=%d\n",A_a1->nb, bbai, c_ai);
            if ((*A).m[c_mi][0].r[c_ri].a[bbai].N[1]>a4){ // if the atom number is > number in a4
                a3=a4; // push what was a4 down to a3. These are just for atom number comparison.
                a3_ai=prev;
                a4=(*A).m[c_mi][0].r[c_ri].a[bbai].N[1]; // make a4 = the current atom number
                a4_ai=bbai;
                prev=bbai;
            }
            else if ((*A).m[c_mi][0].r[c_ri].a[bbai].N[1]>a3){ // if !> a4 check if > a3
                a3=(*A).m[c_mi][0].r[c_ri].a[bbai].N[1]; // make a3 equal the atom number
                a3_ai=bbai;
            }
        }
        A_a3=&(*A).m[c_mi][0].r[c_ri].a[a3_ai];
        A_a4=&(*A).m[c_mi][0].r[c_ri].a[a4_ai];
        strcpy(cplx_a3, A_a3->N); // SAVE ATOM NAME FOR RENAMING COMPLEX
        strcpy(A_a3->N,"a3"); // rename that atom a3
        strcpy(cplx_a4, A_a4->N); // SAVE ATOM NAME FOR RENAMING COMPLEX
        strcpy(A_a4->N,"a4"); // rename that atom a4
    }
}
strcpy(filename,"checknamesA.pdb");
outputAsmblPDB(A, filename);
a4=0; //reset

// Get ID number from fileprefix
i=25; //Need to start here because of ./RESULTS_library_search/ in name.
char id[20];
while( (fileprefix[i] != '-') && (i < 45) ){
    id[i-25]=fileprefix[i];
    i++;
}
id[i-25]='\0';
printf("Glycan Id=%s\n",id);

//Create a list of IDs that were grafted. For summary.
FILE *file;
file = fopen("grafted-IDs.txt", "a+");
fprintf(file,"%s\n",id);
fclose(file);


//Read through spacer_id.txt intil we find the spacer used by the ID.
char spacer[20], buffer[20], line[80];
int stop=0;
file = fopen("spacer_id.txt", "rt");
while( (fgets(line, 80, file) !=NULL) && (stop!=1) ){
    sscanf(line,"%s",buffer);
    if(strcmp(id,buffer)==0){
        stop=1;
        fgets(line, 80, file);
        sscanf(line,"%s",spacer);
        printf("Success %s=%s\n",id,spacer);
    }
    //printf("%s\n",buffer);
}
fclose(file);

// Put header into Results.txt
FILE * pFile;
pFile = fopen (spaceroutname,"a+"); // open file and append
fprintf(pFile,"Grafting %s\nID-Pop|  vdW  | RBP | Clashing Receptor Residues\n",spacer);
fclose (pFile);

// Go into Spacers folder. Find current spacer and create a list of all structural files associated with it.
char dir[40];
//char spacerfile[50];
sprintf(dir, "/home/oliver/work/linkerLib/%s%c/",spacer,ano);
//strcpy (dir,"/home/oliver/work/linkerLib/");
//strcat (dir,spacer);
//strcat (dir,"/");
//int len=strlen(dir);
//dir[len-1]= '/';
printf("dir=%s\n",dir);

int numfiles=0;
struct dirent **namelist;
int n;
char (*d_list)[200];
d_list = malloc(1000 * sizeof *d_list);
n = scandir(dir, &namelist, 0, alphasort);
if (n < 0) { perror("scandir");}
else {
    numfiles=n; //remove the "." and ".." names
    while (n--) {
        printf("%s\n", namelist[n]->d_name);
        sprintf(d_list[n],"%s%s",dir,namelist[n]->d_name);
        //free(namelist[n]);
        //numfiles++;
    }
    //free(namelist);
}
//n=numfiles;

/*DIR *mydir = opendir(dir);
struct dirent *entry = NULL;
int numfiles=0;
i=0;
char (*d_list)[200];
d_list = malloc(1000 * sizeof *d_list);
while((entry = readdir(mydir))){ // If we get EOF, the expression is 0 and the loop stops. 
    if ( (strcmp(entry->d_name,".")!=0) && (strcmp(entry->d_name,"..")!=0) ){ // skip . and ..
        printf("%s\n", entry->d_name);
        strcpy(d_list[i], dir);
        strcat(d_list[i], entry->d_name);
        i++;
    }
/}
//closedir(mydir);
numfiles=i;
*/

//for (i=0;i<(numfiles);i++){
//        printf("d_list[%d]=%s\n",i,d_list[i]);
//}

//printf("numfiles=%d\n",numfiles);


// For each spacer file in list
//load into B and rename atoms for superimpose4atoms.
for(i=2;i<(numfiles);i++){ // start at 2 to skip . and .. directories
    //snprintf(spacerfile, "%s%s",dir,namelist[i]->d_name);
    //printf("spacerfile=%s\n",spacerfile);
    load_pdb(d_list[i],&Sp);
    set_nbonds_for_atoms_in_assembly(&Sp);
    for (ri=0;ri<Sp.m[0][0].nr;ri++){
        if (Sp.m[0][0].r[ri].na>1){
            set_smallest_rings_from_residue_atom_nodes(&Sp.m[0][0].r[ri]);
        }
    }
// superimpose 4 atoms
    for(mi=0;mi<Sp.nm;mi++){ 
        for(ri=0;ri<Sp.m[mi][0].nr;ri++){ 
            printf("Sp=%s,Residue=%s\n",spacer,Sp.m[mi][0].r[ri].N);
            if( (strcmp(Sp.m[mi][0].r[ri].N,spacer)==0) || (strcmp(Sp.m[mi][0].r[ri].N,"NLN")==0) ) {
                printf("Sp contains spacer %s\n",Sp.m[mi][0].r[ri].N);
                for(ai=0;ai<Sp.m[mi][0].r[ri].na;ai++){ // for each atom in spacer
                    for(j=0;j<Sp.m[mi][0].r[ri].a[ai].nmb;j++){
                        bai=Sp.m[mi][0].r[ri].a[ai].mb[j].t.a;
                        bri=Sp.m[mi][0].r[ri].a[ai].mb[j].t.r;
                        bmi=Sp.m[mi][0].r[ri].a[ai].mb[j].t.m;
                        //printf("%s.%d.%s.nmb=%d\n",Sp.m[mi][0].r[ri].N,Sp.m[mi][0].r[ri].n,Sp.m[mi][0].r[ri].a[ai].N,Sp.m[mi][0].r[ri].a[ai].nmb);
                        //printf("ai=%d,ai.N=%s,ri=%d,bri=%d,bri.N=%s,nmb=%d\n",ai,Sp.m[mi][0].r[ri].a[ai].N,ri,bri,Sp.m[bmi][0].r[bri].N,Sp.m[mi][0].r[ri].a[ai].nmb);
                        // If not in same residue and not LNK residue must be sugar :(
                        //if( (ri!=bri) && (strcmp(Sp.m[bmi][0].r[bri].N,"LNK")!=0) ) {
                        // If not in same residue and atom is a ring atom must be sugar :(
                        if( (ri!=bri) && (Sp.m[bmi][0].r[bri].a[bai].R=='Y') ) {
                            Sp_a1=&Sp.m[bmi][0].r[bri].a[bai];
                            printf("%s.%d.%s=a1\n",Sp.m[bmi][0].r[bri].N,Sp.m[bmi][0].r[bri].n,Sp.m[bmi][0].r[bri].a[bai].N);
                            Sp_a2=&Sp.m[mi][0].r[ri].a[ai];
                            printf("Sp_a1->N=%s,Sp_a2->N=%s\n",Sp_a1->N,Sp_a2->N);
                            strcpy(Sp_a1->N, "a1");
                            strcpy(Sp_a2->N, "a2");
                            //a1_ai=bai;
                            a1_ri=bri;
                            a1_mi=bmi;
                            //Now stop searching by maxing counters
//                            mi=Sp.nm;
//                            ri=Sp.m[mi][0].nr;
//                            ai=Sp.m[mi][0].r[ri].na;
//                            j=Sp.m[mi][0].r[ri].a[ai].nmb;
                        }
                    }
                }
            }
        }
    }
    a4=0; //reset
    prev=0; //reset
    //int tmp=0;
    printf("Will now check for a3 and a4\n");
    for(j=0;j<Sp_a1->nb;j++){ // for each intra-residue(nb) atom bonded to a1
        bbai=Sp_a1->b[j].t.a; // readability: save ai of bonded atom
        //tmp=Sp_a1->b[j].t.r;
        //printf("nb=%d,bbri=%d,bbai=%d,bbai.N=%s\n",Sp_a1->nb,tmp,bbai,Sp.m[a1_mi][0].r[a1_ri].a[bbai].N);
        printf("Sp_a1.nb=%d,Sp_a1->b[%d]=%s\n",Sp_a1->nb,j,Sp.m[a1_mi][0].r[a1_ri].a[bbai].N);
        if( Sp.m[a1_mi][0].r[a1_ri].a[bbai].N[1]>a4 ){ // if the atom number is greater than number in a4
            a3=a4; // push what was a4 down to a3. These variables are just for atom number comparisons.
            a3_ai=prev;
            a4=Sp.m[a1_mi][0].r[a1_ri].a[bbai].N[1];
            a4_ai=bbai;
            prev=bbai;
        }
        else if ( Sp.m[a1_mi][0].r[a1_ri].a[bbai].N[1]>a3 ){ // if not greater than a4 but > a3
            a3=Sp.m[a1_mi][0].r[a1_ri].a[bbai].N[1];
            a3_ai=bbai;
        } 
    }
    Sp_a3=&Sp.m[a1_mi][0].r[a1_ri].a[a3_ai];
    Sp_a4=&Sp.m[a1_mi][0].r[a1_ri].a[a4_ai];
    printf("Linker a3=%s, a4=%s\n",Sp_a3->N,Sp_a4->N);

    //strcpy(cplx_a3, Sp_a3->N); // SAVE ATOM NAME FOR RENAMING
    strcpy(Sp_a3->N, "a3"); // rename that atom a3
    //strcpy(cplx_a4, Sp_a4->N); // SAVE ATOM NAME FOR RENAMING
    strcpy(Sp_a4->N, "a4"); // rename that atom a4

    strcpy(filename,"checknamesB.pdb");
    outputAsmblPDB(&Sp, filename);

    superimpose4atoms(&Sp,A);
    
    //Check clash with protein and write out.
    //Check residues below plane if no clash with protein.
    
    find_vdw_clashes_spacers(&Sp,B,namelist[i]->d_name,spaceroutname, &fitterCnt);
    //free(namelist[n]);    

    sprintf(pdboutname, "%s-%dSp.pdb",fileprefix,count); // outputs a file for each change
    printf("outputfile=%s\n",pdboutname);
    outputAsmblPDB_NoRenameResid(&Sp,pdboutname);
    count++;

}
pFile = fopen (spaceroutname,"a+"); // open file and append
fprintf(pFile,"TOTAL fitting spacers=%d\n",fitterCnt);
fclose (pFile);
//free(namelist);
return 0;
}
