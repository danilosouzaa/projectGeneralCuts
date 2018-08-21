
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

int* calculateViolation(Cut_gpu *h_cut){
    int i, j, el;
    int *vViolation = (int*)malloc(sizeof(int)*h_cut->numberConstrains);
    for(i=0;i<h_cut->numberConstrains;i++){
        vViolation[i] = 0;
        for(j = h_cut->ElementsConstraints[i];j<h_cut->ElementsConstraints[i+1];j++){
            el = h_cut->Elements[j];
            vViolation[i] += h_cut->Coefficients[j]*h_cut->xAsterisc[el];
        }
        vViolation[i] = h_cut->rightSide[i] - vViolation[i];
        printf("Violation: %d\n", vViolation[i]);
    }
    return vViolation;


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
    int *vViolation = calculateViolation(h_cut);
    getchar();
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
    free(idxOriginal);
    free(vAux);
    free(vViolation);
    return 0;
}
