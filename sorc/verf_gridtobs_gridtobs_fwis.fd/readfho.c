#include <stdio.h>
#include <stdlib.h>

void readfho(int *numthr, float t[])
{
     FILE *fp;
     int i, num;
     char stat;
     float iii, rnum; 
										
     setbuf(stdout,NULL); /* set output to unbuffered */
     fp = fopen("thresholds", "r"); 


        fscanf(fp, "%f%s%f", &iii, &stat, &rnum);    
/*        printf("%f%s%f\n", iii, stat, rnum);     */
        num = rnum;
        printf("%f\n", iii) ;
/*        printf("%c\n", stat) ;     */
/*        printf("%f\n", rnum);      */   
/*        printf("%d\n", num);       */
    for (i = 0; i <= num-1; i++)
     { 
        fscanf(fp, "%f", &iii);  
        t[i] = iii;
/*        printf("%f%f\n", iii,t[i]);   */   
     }         
        *numthr = i;

    fclose(fp);   

}
