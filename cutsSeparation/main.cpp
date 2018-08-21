
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>

extern "C"
{
	#include "lp.h"
	#include "cut_gpu.h"
    #include "solutionGpu.h"
    #include "prepareGpu.h"
}

/*Pameters
1: Name Instance
2: Precision of decimal
3: Group Size initial
4: Number of max denominator phase 2



*/
int generateVetorSec(Cut_gpu *h_cut, int *vAux, int sz){
    int el, cont = 0,i=0,j=0;
    while(cont<sz){
        for(i = h_cut->ElementsConstraints[j];i<h_cut->ElementsConstraints[j+1];i++){
            el = h_cut->Elements[i];
            if((h_cut->Coefficients[i]!=0)&&(h_cut->xAsterisc[el]!=0)){
               vAux[cont] = j;
               cont++;
               break;
            }

        }
        j++;
        if(j==h_cut->numberConstrains){
            printf("Number Constraints invalided!");
            return -1;
        }
    }
    return 0;
}

int main(int argc, const char *argv[])
{
    if(argc<5){
        printf("Number of parameters invalided\n");
        return 0;
    }
    char name[255] = "../inst/";
    int precision = atoi(argv[2]);
    int szGroup1 = atoi(argv[3]);
    int maxDenomitor =atoi(argv[4]);
    int i;
    printf("%d %d\n",precision,szGroup1);
    getchar();
    strcat(name,argv[1]);
    LinearProgram *lp = lp_create();
    lp_read(lp,name);
    Cut_gpu *h_cut = fillStructPerLP(precision, lp);
    int *vAux = (int*)malloc(sizeof(int)*szGroup1);
    int *idxOriginal;
    int valided;
    valided  = generateVetorSec(h_cut,vAux,szGroup1);
//    for(i=0;i<szGroup1;i++){
//        vAux[i] = i;
//    }

    Cut_gpu *h_cut_group = CreateGroupForVectorNumberConstraints(h_cut,vAux,szGroup1,idxOriginal);
    for(i=0;i<szGroup1;i++){
        show_contraints(h_cut_group,i);
    }
    getchar();
    int cutIni = h_cut_group->numberConstrains;
    int nBlock = 2, nThread = 512;
    h_cut_group = initial_runGPU(h_cut_group,precision,nBlock,nThread,szGroup1);
    h_cut_group = zeroHalf_runGPU(h_cut_group,szGroup1,precision,nThread,nBlock);
    printf("numero inicial: %d, depois : %d\n",cutIni, h_cut_group->numberConstrains) ;
   //s h_cut_group =
    //getchar();
    lp_set_max_seconds( lp, 30 );
    lp_optimize( lp );

//    lp_set_max_seconds(lp,60);
//    lp_optimize(lp);
//    lp_free(&lp);
    free(h_cut);
    lp_free( &lp );
    lp_close_env();

    return 0;
}
