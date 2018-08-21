
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
    for(i=0;i<szGroup1;i++){
        vAux[i] = i;
    }
    Cut_gpu *h_cut_group = CreateGroupForVectorNumberConstraints(h_cut,vAux,szGroup1,idxOriginal);
    for(i=0;i<szGroup1;i++){
        show_contraints(h_cut_group,i);
    }
    getchar();
    h_cut_group = initial_runGPU(h_cut_group,maxDenomitor,precision,1,100,2,szGroup1);

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
