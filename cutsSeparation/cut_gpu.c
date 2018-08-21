#include "cut_gpu.h"

Cut_gpu *AllocationStructCut(int cont, int nConstrains, int nVariables)
{
    size_t size_cut = sizeof(Cut_gpu) +
                      sizeof(TCoefficients)*(cont) +
                      sizeof(TElements)*(cont) +
                      sizeof(TElementsConstraints)*(nConstrains+1) +
                      sizeof(TRightSide)*(nConstrains) +
                      sizeof(TXAsterisc)*(nVariables)+
                      sizeof(TTypeConstraints)*(nConstrains);


    Cut_gpu *cut = (Cut_gpu*)malloc(size_cut);
    assert(cut!=NULL);
    memset(cut,0,size_cut);
    cut->Coefficients = (TCoefficients*)(cut+1);
    cut->Elements = (TElements*)(cut->Coefficients + cont);
    cut->ElementsConstraints = (TElementsConstraints*)(cut->Elements + cont);
    cut->rightSide = (TRightSide*)(cut->ElementsConstraints + (nConstrains+1));
    cut->xAsterisc = (TXAsterisc*)(cut->rightSide + (nConstrains));
    cut->typeConstraints = (TTypeConstraints*)(cut->xAsterisc + (nVariables));
    cut->numberVariables = nVariables;
    cut->numberConstrains = nConstrains;
    cut->cont = cont;
    return cut;
}


Cut_gpu_aux *AllocationStructCutAux(int nConstrains, int nVariables, int nCont)
{
    size_t size_cut_aux = sizeof(Cut_gpu_aux) +
                          sizeof(TInterval)*(nConstrains)+
                          sizeof(TInterval)*(nConstrains)+
                          sizeof(TNames)*(nCont)+
                          sizeof(TNames)*(nConstrains);


    Cut_gpu_aux *cut_aux = (Cut_gpu_aux*)malloc(size_cut_aux);
    assert(cut_aux!=NULL);
    memset(cut_aux,0,size_cut_aux);
    cut_aux->intervalMin = (TInterval*)(cut_aux + 1);
    cut_aux->intervalMax = (TInterval*)(cut_aux->intervalMin + (nConstrains));
    cut_aux->nameElements = (TNames*)(cut_aux->intervalMax + (nConstrains));
    cut_aux->nameConstraints = (TNames*)(cut_aux->nameElements+(nCont));
    cut_aux->numberVariables = nVariables;
    cut_aux->numberConstrains = nConstrains;
    return cut_aux;

}

listNeigh *AllocationListNeigh(int nConstrains, int nList)
{
    size_t size_list = sizeof(listNeigh) +
                       sizeof(TList)*(nList)+
                       sizeof(TPosList)*(nConstrains+1);
    listNeigh *list_t = (listNeigh*)malloc(size_list);
    assert(list_t!=NULL);
    memset(list_t, 0,size_list);
    list_t->list_n = (TList*)(list_t + 1);
    list_t->pos = (TPosList*)(list_t->list_n + nList);
    list_t->nList = nList;
    list_t->nPos = nConstrains + 1;
    return list_t;
}

int returnIndVector(TNames *v,char *nome, int sz)
{
    int i;
    for(i=0; i<sz; i++)
    {
        if(strcmp(v[i].name,nome)==0)
            return i;
    }
    return -1;
}

/*Cut_gpu* fillStruct(CutCG *ccg_r2,int precision, int numberVariables, int numberConstrains, int cont)
{
   Cut_gpu *ccg =  AllocationStructCut(cont, numberConstrains, numberVariables);
   int i,pos = 0, el, elem;
   ccg->ElementsConstraints[0] = 0;
   for(i = 0; i < VStr_size(ccg_r2->rname); i++)
   {
       ccg->rightSide[i] = VDbl_get(ccg_r2->rowrhs,i);
       ccg->typeConstraints[i] = VInt_get(ccg_r2->rowtype,i);
       if(i<numberConstrains-1)
           ccg->ElementsConstraints[i+1] = ccg->ElementsConstraints[i] + VDbl_size(ccg_r2->rowCoef[i]);
       else
           ccg->ElementsConstraints[i+1] = ccg->cont;
       for( el = 0 ; el < VDbl_size(ccg_r2->rowCoef[i]) ; el++)
       {
           elem = VInt_get(ccg_r2->rowElem[i],el);
           ccg->Coefficients[pos] = VDbl_get(ccg_r2->rowCoef[i],el);
           ccg->Elements[pos] = elem;
           pos++;
       }

   }
   for(i=0; i<ccg->numberVariables; i++)
   {
       ccg->xAsterisc[i] = precision * VDbl_get(ccg_r2->xfElemPP,i);
   }
   return ccg;
}
*/

Cut_gpu* fillStructPerLP(int precision, LinearProgram *lp)
{
    int numberVariables, numberConstrains, numberNonZero;
    int i, j, k;
    numberVariables = lp_cols(lp);
    numberConstrains = lp_rows(lp);
    numberNonZero = lp_nz(lp);
    int *idx = (int*)malloc(sizeof(int)*numberVariables);
    double *coef = (double*)malloc(sizeof(double)*numberVariables);
    double rhs;
    double *xTemp;
    //printf("testes: %d %d %d", numberConstrains, numberVariables, numberNonZero);
    lp_optimize_as_continuous(lp);
    xTemp = lp_x(lp);
    Cut_gpu *h_cut;
    h_cut = AllocationStructCut(numberNonZero,numberConstrains,numberVariables);
    for(i=0; i<numberVariables; i++)
    {
        coef[i]=0.0;
        idx[i] = 0;
        h_cut->xAsterisc[i] = xTemp[i] * precision;
        //printf("x = %d", h_cut->xAsterisc[i]);
        //getchar();
    }
    //getchar();
   // lp_row(lp,0,idx,coef);
    int aux = 0;
    h_cut->ElementsConstraints[0] = 0;
    for(i=0; i<numberConstrains; i++)
    {
        h_cut->typeConstraints[i] = 0;
        lp_row(lp,i,idx,coef);
        for(j=0; j<numberVariables; j++)
        {
            if(coef[j]!=0.0)
            {
              //  printf("Coef: %f idx: %d \t ", coef[j],idx[j]);
                h_cut->Coefficients[aux] = coef[j];
                h_cut->Elements[aux] = idx[j];
                aux++;
            }
            coef[j] = 0.0;
            idx[j] = 0;
        }
        h_cut->ElementsConstraints[i+1] = aux;
        rhs = lp_rhs(lp,i);
        h_cut->rightSide[i] = rhs;
        //printf("rhs: %f\n", rhs);
    }
    //printf("AUx: %d total: %d\n", aux, numberNonZero);
    //getchar();
    free(idx);
    free(coef);
    //free(xTemp);
    return h_cut;
}

Cut_gpu *CreateGroupForVectorNumberConstraints(Cut_gpu *h_cut, int *vectorConstraints, int szConstraints, int *idxOriginal){
        int cont = 0, aux = 0;
        int i, j, el, constraint;
        int *vAux = (int*)malloc(sizeof(int)*h_cut->numberVariables);
        memset(vAux,0,sizeof(int)*h_cut->numberVariables);
        for(i=0;i<szConstraints;i++){
            constraint = vectorConstraints[i];
            for(j=h_cut->ElementsConstraints[constraint];j<h_cut->ElementsConstraints[constraint+1];j++){
                el = h_cut->Elements[j];
                if(h_cut->xAsterisc[el]!=0){
                    cont++;
                    vAux[el] = 1;
                }
            }
        }
        for(i=0;i<h_cut->numberVariables;i++){
            if(vAux[i]==1){
                aux++;
            }
        }
        idxOriginal = (int*)malloc(sizeof(int)*aux);
        Cut_gpu *h_cut_group = AllocationStructCut(cont,szConstraints,aux);
        aux = 0;
        for(i=0;i<h_cut->numberConstrains;i++){
            if(vAux[i]==1){
                idxOriginal[aux] = i;
                h_cut_group->xAsterisc[aux] = h_cut->xAsterisc[i];
                aux++;
            }
        }
        cont = 0;
        h_cut_group->ElementsConstraints[0] = 0;
        for(i=0;i<szConstraints;i++){
            constraint = vectorConstraints[i];
            for(j=h_cut->ElementsConstraints[constraint];j<h_cut->ElementsConstraints[constraint+1];j++){
                el = h_cut->Elements[j];
                if(vAux[el]==1){
                    h_cut_group->Coefficients[cont] = h_cut->Coefficients[j];
                    h_cut_group->Elements[cont] = h_cut->Elements[j];
                    cont++;
                }
            }
            h_cut_group->rightSide[i] = h_cut->rightSide[constraint];
            h_cut_group->ElementsConstraints[i+1] = cont;
            h_cut_group->typeConstraints[i] = RES_RR; //Só para testes.
        }
        return h_cut_group;
}

void setParameters_ccg(parameters_ccg *parCCG, int mode)
{
    parCCG->precision=1000;
    if(mode==0)
    {
        parCCG->numberMaxConst = 8; //m'
        parCCG->nRuns  = 7000; //it
        parCCG->maxDenominator = 100; //max denom
        parCCG->nSizeR1 = 6; //n restricao de cover res
    }
    else
    {
        parCCG->numberMaxConst = 16;
        parCCG->nRuns  = 50000;
        parCCG->maxDenominator = 100;
        parCCG->nSizeR1 = 12;
    }
//    printf("parameters: \n precision: %d, numberMax: %d, nRuns: %d, MaxDen: %d, sizeR1: %d\n",parCCG->precision,parCCG->numberMaxConst,parCCG->nRuns,parCCG->maxDenominator,parCCG->nSizeR1);
//    getchar();

}

