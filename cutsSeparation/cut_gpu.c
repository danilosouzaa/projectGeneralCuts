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

Cover_gpu *AllocationStructCover(int cont, int nConstraints)
{
    size_t size_cover = sizeof(Cover_gpu) +
                        sizeof(TCoefficients)*(cont) +
                        sizeof(TElementsConstraints)*(nConstraints+1) +
                        sizeof(TRightSide)*(nConstraints);
    Cover_gpu *h_cover = (Cover_gpu*)malloc(size_cover);
    assert(h_cover!=NULL);
    memset(h_cover,0,size_cover);
    h_cover->Coefficients = (TCoefficients*)(h_cover+1);
    h_cover->ElementsConstraints = (TElementsConstraints*)(h_cover->Coefficients + cont);
    h_cover->rightSide = (TRightSide*) (h_cover->ElementsConstraints + (nConstraints+1) );
    h_cover->numberConstraints = nConstraints;
    h_cover->cont = cont;
    return h_cover;
}


Cover_gpu *CopyCutToCover(Cut_gpu *h_cut, int nConstraintsIni)
{
    //Cover_gpu *cover = AllocationStructCover(h_cut->cont,h_cut->numberConstrains);
    Cover_gpu *cover = AllocationStructCover(h_cut->cont,nConstraintsIni);
    int i;
    for(i=0; i< h_cut->cont; i++)
    {
        cover->Coefficients[i] = h_cut->Coefficients[i];
    }
    for(i=0; i<nConstraintsIni; i++)
    {
        cover->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
        cover->rightSide[i] = h_cut->rightSide[i];
    }
    cover->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
    return cover;
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
    int i, j, k, numberConstrainsLeft = 0, numberNonZeroNew = 0;
    numberVariables = lp_cols(lp);
    numberConstrains = lp_rows(lp);
    numberNonZero = lp_nz(lp);
    int flag;
    int *idx = (int*)malloc(sizeof(int)*numberVariables);
    double *coef = (double*)malloc(sizeof(double)*numberVariables);
    for(i=0; i<numberConstrains; i++)
    {
        if(lp_sense(lp,i)=='L')
        {
            numberConstrainsLeft++;
            for(j=0; j<numberVariables; j++)
            {
                coef[j]=0.0;
                idx[j] = 0;
            }
            lp_row(lp,i,idx,coef);
            flag = 0;
            for(j=0; j<numberVariables; j++)
            {
                if(coef[j]!=0)
                {
                    flag++;
                    numberNonZeroNew ++;
                }
            }
            if(flag<2){
                numberConstrainsLeft--;
                numberNonZeroNew -= flag;
            }
        }
    }
//   printf("Novos valores: %d %d", numberConstrainsLeft,numberNonZeroNew);
//    getchar();
    double rhs;
    double *xTemp;
    //printf("testes: %d %d %d", numberConstrains, numberVariables, numberNonZero);
    lp_optimize_as_continuous(lp);
    xTemp = lp_x(lp);
    Cut_gpu *h_cut;
    h_cut = AllocationStructCut(numberNonZeroNew,numberConstrainsLeft,numberVariables);
    int *v_aux = (int*)malloc(sizeof(int)*(numberVariables + 1));
    int tam;
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
    int contador = 0;
    for(i=0; i<numberConstrains; i++)
    {
        //printf("Sense: %c\n",lp_sense(lp,i));
        //getchar();

        if(lp_sense(lp,i)=='L')
        {
            lp_row(lp,i,idx,coef);
            rhs = lp_rhs(lp,i);

            tam = 0;
            for(j=0; j<numberVariables; j++)
            {
                if(coef[j]!=0.0)
                {
                    v_aux[tam] = coef[j];
                    //printf("%d \t",v_aux[tam]);
                    tam++;
                }
            }
            if(tam<2){
                continue;
            }
            v_aux[tam] = rhs;
            //tam++;
            int mdc = CutP_maxDivisorCommonVector(v_aux, tam);
            //printf("Mdc: %d\t",mdc);
            for(j=0; j<numberVariables; j++)
            {
                if(coef[j]!=0.0)
                {

                    //      printf("Coef: %f idx: %d \t ", coef[j],idx[j]);
                    h_cut->Coefficients[aux] = coef[j]/mdc;
                    h_cut->Elements[aux] = idx[j];
                    aux++;
                }
                coef[j] = 0.0;
                idx[j] = 0;
            }

            // printf("\n rhs: %f\n",rhs);

            h_cut->ElementsConstraints[contador+1] = aux;
            h_cut->rightSide[contador] = rhs/mdc;
            h_cut->typeConstraints[contador] = 0;
            contador++;
        }

    }
    //printf("contador: %d aux: %d const: %d cont: %d\n",contador,aux,h_cut->numberConstrains, h_cut->cont);

    //printf("AUx: %d total: %d\n", aux, numberNonZero);
    //getchar();
    free(idx);
    free(coef);
    free(v_aux);

    //free(xTemp);
    //printf("\n");
    //for(i=0;i<=h_cut->numberConstrains;i++){
    //    printf("%d \t", h_cut->ElementsConstraints[i]);
    //}
    //getchar();
    //free(xTemp);
    return h_cut;
}


void quicksortCof(int *values, int *idc, int began, int end)
{
    int i, j, pivo, aux;
    i = began;
    j = end-1;
    pivo = values[(began + end) / 2];
    while(i <= j)
    {
        while(values[i] > pivo && i < end)
        {
            i++;
        }
        while(values[j] < pivo && j > began)
        {
            j--;
        }
        if(i <= j)
        {
            aux = values[i];
            values[i] = values[j];
            values[j] = aux;

            aux = idc[i];
            idc[i] = idc[j];
            idc[j] = aux;

            i++;
            j--;
        }
    }
    if(j > began)
        quicksortCof(values, idc,  began, j+1);
    if(i < end)
        quicksortCof(values, idc, i, end);
}


void SortByCoefficients(Cut_gpu *h_cut)
{
    int i = 0, j = 0;
    int el;
    for(i = 0; i < h_cut->numberConstrains; i++)
    {
        quicksortCof(h_cut->Coefficients, h_cut->Elements, h_cut->ElementsConstraints[i], h_cut->ElementsConstraints[i+1]);
    }
}

int insertConstraintsLP(LinearProgramPtr lp, Cut_gpu *h_cut, int nConstrainsInitial, int *counterCuts)
{
    int nRows, i, j, w, k = 0 ;
    int *idx;
    double *Coef;
    char name[255];
    for(i = nConstrainsInitial; i < h_cut->numberConstrains; i++)
    {
        int sz = h_cut->ElementsConstraints[i+1] - h_cut->ElementsConstraints[i];
        idx = (int*)malloc(sizeof(int)*sz);
        Coef = (double*)malloc(sizeof(double)*sz);
        for(j = h_cut->ElementsConstraints[i]; j<h_cut->ElementsConstraints[i+1]; j++)
        {
            idx[w] = h_cut->Elements[j];
            Coef[w] = h_cut->Coefficients[j];
            w++;
        }
        double rhs = h_cut->rightSide[i];

        //if(h_cut->typeConstraints[i]==LPC_CCOVER){
        int v = *counterCuts;
        sprintf(name, "CCOVER(%d)",v);
        //printf("%s\n",name);
        (*counterCuts)++;
        lp_add_row(lp,sz,idx,Coef,name,'L',rhs);
        //}
        w = 0;
        free(idx);
        free(Coef);
    }
}

Cut_gpu *removeNegativeCoefficientsAndSort(Cut_gpu *h_cut, int *convertVector, int precision)
{
    int i,j;
    //convertVector = (int*)(malloc(sizeof(int)*h_cut->cont));
    //  convertCoef = (int*)(malloc(sizeof(int)*h_cut->cont));
    int qntX = h_cut->numberVariables;
    int qntNegative = 0, rhs = 0;
    int el = 0;
    for(i=0; i< h_cut->cont; i++)
    {
        if(h_cut->Coefficients[i]<0)
        {
            qntNegative++;
        }
    }
    qntNegative += h_cut->numberVariables;
    Cut_gpu *Cut_new = AllocationStructCut(h_cut->cont,h_cut->numberConstrains, qntNegative);

    for(i=0; i<h_cut->numberVariables; i++)
    {
        Cut_new->xAsterisc[i] = h_cut->xAsterisc[i];
    }


    for(i = 0; i<h_cut->numberConstrains; i++)
    {
        rhs = h_cut->rightSide[i];
        for(j= h_cut->ElementsConstraints[i]; j<h_cut->ElementsConstraints[i+1]; j++)
        {
            if(h_cut->Coefficients[j]<0)
            {
                Cut_new->Coefficients[j] = h_cut->Coefficients[j]*(-1);
//                printf("%d %d \t", Cut_new->Coefficients[j], j);

                rhs += Cut_new->Coefficients[j];
                Cut_new->Elements[j] = qntX;
                el = h_cut->Elements[j];
                Cut_new->xAsterisc[qntX] = precision - h_cut->xAsterisc[el];
                convertVector[qntX-h_cut->numberVariables] = h_cut->Elements[j];
                //convertCoef[qntX-h_cut->numberVariables] = Cut_new->Coefficients[j];
                qntX++;
            }
            else
            {
                Cut_new->Coefficients[j] = h_cut->Coefficients[j];
                Cut_new->Elements[j] = h_cut->Elements[j];
                //convertVector[j] = 0 ;
                //convertCoef[j] = 0;
            }
        }
        Cut_new->rightSide[i] = rhs;
        Cut_new->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
        Cut_new->typeConstraints[i] = h_cut->typeConstraints[i];
    }
    Cut_new->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
    SortByCoefficients(Cut_new);
    free(h_cut);
    return Cut_new;
}


Cut_gpu *returnVariablesOriginals(Cut_gpu *h_cut, int *convertVector, int precision, int nVariablesInitial )
{
    Cut_gpu *Cut_new = AllocationStructCut(h_cut->cont, h_cut->numberConstrains, nVariablesInitial);
    int i, j, rhs = 0, el;

    for(i=0; i<nVariablesInitial; i++)
    {
        Cut_new->xAsterisc[i] = h_cut->xAsterisc[i];
    }
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        rhs = h_cut->rightSide[i];
//        printf("\nrhs %d",rhs);
//        getchar();
        for(j= h_cut->ElementsConstraints[i]; j<h_cut->ElementsConstraints[i+1]; j++)
        {
            el = h_cut->Elements[j];
            if(el>=nVariablesInitial)
            {
                Cut_new->Coefficients[j] = h_cut->Coefficients[j]*(-1);
//                printf("%d %d %d \n",h_cut->Coefficients[j], Cut_new->Coefficients[j], j);
                rhs -= h_cut->Coefficients[j];
                Cut_new->Elements[j] = convertVector[el - nVariablesInitial];
                Cut_new->xAsterisc[ Cut_new->Elements[j] ] = precision - h_cut->xAsterisc[el];
            }
            else
            {
                Cut_new->Coefficients[j] = h_cut->Coefficients[j];
                Cut_new->Elements[j] = h_cut->Elements[j];

            }
        }
//        printf("rhs depois: %d\n",rhs);
        Cut_new->rightSide[i] = rhs;
        Cut_new->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
        Cut_new->typeConstraints[i] = h_cut->typeConstraints[i];

    }
    Cut_new->ElementsConstraints[i] = h_cut->ElementsConstraints[i];
    free(h_cut);
//    free(h_cut);
    return Cut_new;


}

int returnPos(int *v,int n, int x)
{
    int e, m, d;                              // 1
    e = 0;
    d = n-1;                           // 2
    while (e <= d)                            // 3
    {
        m = (e + d)/2;                         // 4
        if (v[m] == x) return m;               // 5
        if (v[m] < x) e = m + 1;               // 6
        else d = m - 1;                        // 7
    }                                         // 8
    return -1;                                // 9
}

Cut_gpu *CreateGroupForVectorNumberConstraints(Cut_gpu *h_cut, int *vectorConstraints, int szConstraints, int *idxOriginal)
{
    int cont = 0, aux = 0;
    int i, j, el, constraint;
    int *vAux = (int*)malloc(sizeof(int)*h_cut->numberVariables);
    memset(vAux,0,sizeof(int)*h_cut->numberVariables);
    for(i=0; i<szConstraints; i++)
    {
        constraint = vectorConstraints[i];
        for(j=h_cut->ElementsConstraints[constraint]; j<h_cut->ElementsConstraints[constraint+1]; j++)
        {
            el = h_cut->Elements[j];
            if(h_cut->xAsterisc[el]!=0)
            {
                cont++;
                vAux[el] = 1;
            }
        }
    }
    for(i=0; i<h_cut->numberVariables; i++)
    {
        if(vAux[i]==1)
        {
            aux++;
        }
    }
    idxOriginal = (int*)malloc(sizeof(int)*aux);
    Cut_gpu *h_cut_group = AllocationStructCut(cont,szConstraints,aux);
    aux = 0;
    for(i=0; i<h_cut->numberVariables; i++)
    {
        if(vAux[i]==1)
        {
            idxOriginal[aux] = i;
            h_cut_group->xAsterisc[aux] = h_cut->xAsterisc[i];
            aux++;
        }
    }
    cont = 0;
    h_cut_group->ElementsConstraints[0] = 0;
    for(i=0; i<szConstraints; i++)
    {
        constraint = vectorConstraints[i];
        for(j=h_cut->ElementsConstraints[constraint]; j<h_cut->ElementsConstraints[constraint+1]; j++)
        {
            el = h_cut->Elements[j];
            if(vAux[el]==1)
            {
                h_cut_group->Coefficients[cont] = h_cut->Coefficients[j];
                h_cut_group->Elements[cont] = returnPos(idxOriginal,aux,el);//h_cut->Elements[j];
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

int generateVetorSec(Cut_gpu *h_cut, int *vAux, int sz)
{
    int el, cont = 0,i=0,j=0;
    while(cont<sz)
    {
        for(i = h_cut->ElementsConstraints[j]; i<h_cut->ElementsConstraints[j+1]; i++)
        {
            el = h_cut->Elements[i];
            if((h_cut->Coefficients[i]!=0)&&(h_cut->xAsterisc[el]!=0))
            {
                vAux[cont] = j;
                cont++;
                break;
            }

        }
        j++;
        if(j==h_cut->numberConstrains)
        {
            printf("Number Constraints invalided!");
            return -1;
        }
    }
    return 0;
}

