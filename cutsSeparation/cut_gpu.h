#ifndef CUT_GPU_H_
#define CUT_GPU_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "gpulib/types.h"
#include "lp.h"

EXTERN_C_BEGIN

typedef int TCoefficients;
typedef int TElements;
typedef int TElementsConstraints;
typedef int TRightSide;
typedef int TXAsterisc;
typedef char TTypeConstraints;
typedef short TInterval;
typedef int TList;
typedef int TPosList;

typedef struct{
    char name[255];
}TNames;

typedef struct {
    int numberVariables;
    int numberConstrains;
    int cont;
    TCoefficients *Coefficients;
    TElements *Elements;
    TElementsConstraints *ElementsConstraints;
    TRightSide *rightSide;
    TXAsterisc *xAsterisc;
    TTypeConstraints *typeConstraints;
}Cut_gpu;

typedef struct{
    int numberConstraints;
    int cont;
    TCoefficients *Coefficients;
    TElementsConstraints *ElementsConstraints;
    TRightSide *rightSide;
}Cover_gpu;

typedef struct{
    int numberVariables;
    int numberConstrains;
    TInterval *intervalMin;
    TInterval *intervalMax;
    TNames *nameElements;
    TNames *nameConstraints;
}Cut_gpu_aux;


typedef struct {
    int nList;
    int nPos;
    TList *list_n;
    TPosList *pos;
}listNeigh;

typedef struct {
                int precision;
                int numberMaxConst;
                int nRuns;
                int maxDenominator;
                int nSizeR1;
}parameters_ccg;

Cut_gpu *AllocationStructCut(int cont, int nConstrains, int nVariables);

Cover_gpu *AllocationStructCover(int cont, int nConstraints);

Cover_gpu *CopyCutToCover(Cut_gpu *h_cut, int nConstraintsIni);

Cut_gpu_aux *AllocationStructCutAux(int nConstrains, int nVariables, int nCont);

listNeigh *AllocationListNeigh(int nConstrains, int nList);

//Cut_gpu* readFile(char *fileName, int precision ,Cut_gpu_aux *cut_aux);

//Cut_gpu *readFile(char *fileName, int precision, Cut_gpu_aux* cut_aux);

Cut_gpu* createGPUcut(const Cut_gpu* h_cut, int nVariables, int nConstrains);

Cover_gpu* createGPUcover(const Cover_gpu* h_cover);

listNeigh *createGPUlist(const listNeigh* list_t, int nConstrains, int nList);

int returnIndVector(TNames *v,char *nome, int sz);

//Cut_gpu* fillStruct(CutCG *ccg_r2,int precision, int numberVariables, int numberConstrains, int cont);

Cut_gpu *CreateGroupForVectorNumberConstraints(Cut_gpu *h_cut, int *vectorConstraints, int szConstraints, int *idxOriginal);

Cut_gpu* fillStructPerLP(int precision, LinearProgram *lp);

void setParameters_ccg(parameters_ccg *parCCG, int mode);

int generateVetorSec(Cut_gpu *h_cut, int *vAux, int sz);

//Cut_gpu *removeNegativeCoefficientsAndSort(Cut_gpu *h_cut, int *convertVector, int precision);

Cut_gpu *removeNegativeCoefficientsAndSort(Cut_gpu *h_cut, int *convertVector, int precision);

Cut_gpu *returnVariablesOriginals(Cut_gpu *h_cut, int *convertVector, int precision, int nVariablesInitial);

int insertConstraintsLP(LinearProgramPtr lp, Cut_gpu *h_cut, int nConstrainsInitial, int *counterCuts);

EXTERN_C_END

#endif // CUT_GPU_H_
