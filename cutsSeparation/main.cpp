
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
    #include "prepareCPU.h"
}

/*Pameters
1: Name Instance
2: Precision of decimal
3: Group Size initial
4: Number of max denominator phase 2



*/

int* calculateViolation(Cut_gpu *h_cut,int precision){
    int i, j, el;
    int *vViolation = (int*)malloc(sizeof(int)*h_cut->numberConstrains);
    for(i=0;i<h_cut->numberConstrains;i++){
        //show_contraints(h_cut,i);
        vViolation[i] = 0;
        for(j = h_cut->ElementsConstraints[i];j<h_cut->ElementsConstraints[i+1];j++){
            el = h_cut->Elements[j];
            vViolation[i] += h_cut->Coefficients[j]*h_cut->xAsterisc[el];
        }
        vViolation[i] = (h_cut->rightSide[i]*precision) - vViolation[i];
        //printf("Violation %d: %d %d\n",i, vViolation[i], h_cut->numberConstrains);
//        if(vViolation[i]<=0)
//            getchar();
    }
    return vViolation;
}

int *contNumberPosible(Cut_gpu *h_cut){
    int *qntPerConstraints = (int*)malloc(sizeof(int)*h_cut->numberConstrains);
    int i = 0, aux;
    for(i=0;i<h_cut->numberConstrains;i++){
        aux = h_cut->ElementsConstraints[i+1] - h_cut->ElementsConstraints[i];
        qntPerConstraints[i] = aux;
    }

    return qntPerConstraints;
}

int* sortPerViolation(int *vViolation, int nConstraints){
    int *pos= (int*)malloc(sizeof(int)*nConstraints);
    int i;
    for(i=0;i<nConstraints;i++){
        pos[i] = i;
    }
    bubble_sort(vViolation,pos,nConstraints);
//    for(i = 0;i<nConstraints;i++){
//        printf("%d - %d\n",pos[i],vViolation[i]);
//    }
    return pos;
}

int main(int argc, const char *argv[])
{
    if(argc<6){
        printf("Number of parameters invalided\n");
        return 0;
    }
    char name[255] = "../inst/";
    int precision = atoi(argv[2]);
    int szGroup1 = atoi(argv[3]);
    int maxDenomitor =atoi(argv[4]);
    double timeMax = atof(argv[5]);
    int i,aux = 0;
    int counterCuts=0;
    strcat(name,argv[1]);
    LinearProgram *lp = lp_create();
    lp_read(lp,name);
    Cut_gpu *h_cut = fillStructPerLP(precision, lp);
    lp_write_lp(lp,"danilo.lp");
    int cutIni = h_cut->numberConstrains;
    time_t timeInitial = clock();
    time_t timeFinal;
    double timeCurrent, obj_Best;
    //lp_optimize_as_continuous(lp);

    obj_Best = lp_obj_value(lp);
    printf("Initial obj: %f\n", obj_Best);
    do{

    int *convertVaribles = (int*)malloc(sizeof(int)*h_cut->cont);
    int numberVariablesInitial = h_cut->numberVariables;
    //printf("Number Variables antes: %d\n", h_cut->numberVariables);
    h_cut = removeNegativeCoefficientsAndSort(h_cut,convertVaribles,precision);
    //printf("Number Variables depois: %d\n", h_cut->numberVariables);
    //printf("Number Constraints: %d\n", h_cut->numberConstrains);


    int nBlock = 2, nThread = 512;

    //h_cut_group = initial_runGPU(h_cut_group,precision,nThread,nBlock,szGroup1);
    //h_cut = generateCutsCover(h_cut,10,10);
    int cutIni2 = h_cut->numberConstrains;
    h_cut = runCPU_Cut_Cover(h_cut,1,cutIni);
    h_cut = returnVariablesOriginals(h_cut,convertVaribles,precision,numberVariablesInitial);
    //printf("numero inicial: %d, depois : %d\n",cutIni, h_cut->numberConstrains);
    insertConstraintsLP(lp,h_cut,cutIni2,&counterCuts);
    lp_write_lp(lp,"danilo2.lp");
    lp_set_max_seconds( lp, 60 );
    lp_set_print_messages(lp,0);
    aux = lp_optimize_as_continuous( lp );
//    aux = lp_optimize( lp );
    double newObj = lp_obj_value(lp);
        free(convertVaribles);
    timeFinal = clock();
    timeCurrent = ((double) (timeFinal - timeInitial)) / CLOCKS_PER_SEC;
    if(newObj!=obj_Best){
        printf("New Objetive: %f\n",newObj);
        printf("Number of Inserted Cuts:  %d \n", h_cut->numberConstrains - cutIni);
        printf("Current Time: %f\n", timeCurrent);
        obj_Best = newObj;

    }

    }while(timeCurrent < timeMax);
//    lp_set_max_seconds(lp,60);
//    lp_optimize(lp);
//    lp_free(&lp);
    free(h_cut);
    //free(convertVaribles);
    //free(convertCoef);
    lp_free( &lp );
    lp_close_env();
   // free(idxOriginal);
   // free(vAux);
   // free(vViolation);
   // free(pos);
    return 0;
}
