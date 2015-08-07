#include <glylib.h>


int main (int argc, char *argv[]){
if (argc!=3) {
    printf("Usage: LINK2bond.exe input_filename, output_filename\n");
}
if (argc==3) { // If started correctly

FILE * pFile; 
int count=0;
int i; 

// count the number of lines in the file
pFile = fopen (argv[1],"r");
char dum[200];
if(!pFile) return 1; // bail out if file not found
while(fgets(dum,sizeof(dum),pFile) != NULL){count++;} // count the number of lines in the file
fclose (pFile);

//printf("Count=%d\n",count);

// put all the LINK info into a char array
char LINK_table[count][200];
pFile = fopen (argv[1],"r");

for (i=0;i<count;i++){
	fgets(LINK_table[i],sizeof(LINK_table[i]),pFile);
        }
fclose (pFile);
// Print out the bond info
FILE *file;
file = fopen(argv[2],"w");
for (i=0;i<count;i++){
	//printf("%s\n", LINK_table[i]);
	// Changed as of 2014Mar12 by OG fprintf(file,"bond LIG.%c.%c%c LIG.%c.%c%c\n",LINK_table[i][25],LINK_table[i][13],LINK_table[i][14],LINK_table[i][55],LINK_table[i][43],LINK_table[i][44]);
	fprintf(file,"bond LIG.%c.%c%c LIG.%c.%c%c\n",LINK_table[i][25],LINK_table[i][12],LINK_table[i][13],LINK_table[i][55],LINK_table[i][42],LINK_table[i][43]);
	}
fclose(file); /*done!*/
}
return 0;
}
