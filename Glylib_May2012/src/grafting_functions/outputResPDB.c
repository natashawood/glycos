// hacked by Oliver Grant to only output certain residue numbers
// also checks for an atom named "a9" and does not output it. This is specific for grafting and shouldn't affect anything else.

#include <glylib.h>
#include <molecules.h>

// takes in assembly, name to call the file, list of residue names and number of residues in the list
void outputResPDB(assembly* asmbl, char* file_name, int* residues, int outrn){
//Tells the user what's going on
printf("Writing to file: %s...\n",file_name);
FILE* file;
//char current_line[80] = ""; char temporary[40] = "";
//char* curLine = current_line; char* temp = temporary;
char temp[40]; char curLine[80];
char* resName; char* atmName;
int i,j,k,x,atmTot,resNum,resTot;
atom curAtm; molecule* mol;
file = myfopen(file_name,"w");
for(k = 0; k < (*asmbl).nm; k++){
	mol = ((*asmbl).m+k)[0]; // changed by BLF on 20080622 -- might need revisiting
  	resTot = (*mol).nr;
//	printf("resTot=%d,outrn=%d\n",resTot,outrn);
  	for(i = 0; i < resTot; i++){
		resNum = (*((*mol).r+i)).n;
//		printf("resNum=%d\n",resNum);
		for (x=0; x<outrn;x++){
			if (residues[x]==resNum){ // only output if resNum in the residue list
				resName = (*((*mol).r+i)).N;
				atmTot = (*((*mol).r+i)).na;
				for(j = 0; j < atmTot; j++){
					curAtm = (*((*((*mol).r+i)).a+j));
					atmName = curAtm.N;
					if (strcmp(curAtm.N,"a9")!=0){
						strcpy(curLine,"");
						strcat(curLine,"ATOM"); 
						sprintf(temp,"%d", curAtm.n);
    						//Based on how large the atom # is, depends on the # of spaces
    						strcat(curLine,spacing(strlen(temp),7));
    						//"ATOM atom# atom_name"
						sprintf(temp,"%d  %s", curAtm.n, atmName);
    						strcat(curLine,temp);

			    			//Based on how many char are in the residue name, "   " spaces 
    						strcat(curLine,spacing(strlen(atmName),4));
    						//"ATOM atom# atom_name res_name"
    						strcat(curLine,resName);
   	  			   		sprintf(temp,"%d",resNum);

    				 		//Based on how large the residue # is, "              " spaces
 						strcat(curLine,spacing(strlen(temp),6));
    						//"ATOM atom# atom_name res_name res#"
    						strcat(curLine,temp);
   
   						sprintf(temp,"%.3lf",curAtm.x.i);
  						//Based on how large the x cordinate is, "            " spaces
    						strcat(curLine,spacing(strlen(temp),12));
    						//"ATOM atom# atom_name res_name x"
    						strcat(curLine,temp);
 						sprintf(temp,"%.3lf",curAtm.x.j);
						//Based on how large the y cordinate is, "            " spaces
						strcat(curLine,spacing(strlen(temp),8));
						//"ATOM atom# atom_name res_name x y"
						strcat(curLine,temp);
			
						sprintf(temp,"%.3lf",curAtm.x.k);
						//Based on how large the z cordinate is, "            " spaces
						strcat(curLine,spacing(strlen(temp),8));
	   			 		//"ATOM atom# atom_name res_name x y z 1.00 0.00"
	   			 		sprintf(temp,"%.3lf  1.00  0.00",curAtm.x.k);
	   			 		strcat(curLine,temp);
						//Output the line to the new .pbd file
						fprintf(file,"%s\n",curLine);
//                                                printf("curLine=%s\n",curLine);

					}
  					//fprintf(file,"TER   \n");
				}
			}
		}
	}
}	
fclose(file);
}
char* spacing(int strLen, int totLen)
{
 char returning[15] = ""; char* ret = returning;
 int i;
 for(i = 0; i < totLen - strLen; i++)
  strcat(ret," ");
 return ret;
}
