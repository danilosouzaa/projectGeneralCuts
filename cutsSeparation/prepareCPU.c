#include "prepareCPU.h"

void runCPUR1(Cut_gpu *h_cut, solutionGpu *h_solution, int precision, double timeconst)
{

    int *constraints;
    int pos;
    int violation = 0, i, j;
    pos = 0;
    constraints = (int*)malloc(sizeof(int)*h_cut->numberConstrains);
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        //if((h_cut->typeConstraints[i] == RES_RR) || (h_cut->typeConstraints[i] == LPC_RR) || (h_cut->typeConstraints[i] == RES_NR) || (h_cut->typeConstraints[i] == RES_MODE) || (h_cut->typeConstraints[i] == RES_CONFCL)  )
        // {
        if((h_cut->typeConstraints[i] == RES_RR) )
        {
            constraints[pos] = i;
            pos++;

        }
    }
    int k;
    //  #pragma omp parallel
    // {
    //   #pragma omp for
    for(k = 0; k<pos; k++)
    {
        int res = constraints[k];
        int sz_u = h_cut->ElementsConstraints[res+1] - h_cut->ElementsConstraints[res] + 1;
        int *u = (int*)malloc(sizeof(int)* sz_u );
        int n1=-1, d1=-1,el, rhs, aux,value_tes;
        aux = 0;

        for(j = h_cut->ElementsConstraints[ res ] ; j < h_cut->ElementsConstraints[ res +1 ]; j++)
        {
            u[aux] = h_cut->Coefficients[j];
            aux++;
        }
        u[sz_u - 1] = h_cut->rightSide[res];

        value_tes = sz_u;
        for( i = 0; i < sz_u; i++ )
        {
            for( j = 0; j < sz_u; j++ ) // sempre 1 elemento à frente
            {
                // se o (x > (x+1)) então o x passa pra frente (ordem crescente)


                if((u[i]!=0)&&(u[j]!=0)&&(i!=j))
                {
                    if(u[i]%u[j]==0)
                    {
                        u[j]=0;
                        value_tes--;
                    }
                    if(u[j]%u[i]==0)
                    {
                        u[i]=0;
                        value_tes--;
                    }
                }
                if (( u[i] < u[j] )&&(j>i))
                {
                    aux = u[i];
                    u[i] = u[j];
                    u[j] = aux;
                }
            }
        }
        sz_u = value_tes;

        int nBest=-1, dBest=-1, violation_best=0;
        for(j = 0 ; j < sz_u; j++)
        {
            d1 = u[j];
            n1 = 1;
            while(n1<d1)
            {
                rhs = 0;
                violation = 0;
                value_tes = 0;
                rhs = h_cut->rightSide[ res ]* n1;
                if(rhs%d1 !=0)
                {

                    for(i = h_cut->ElementsConstraints[ res ]; i<h_cut->ElementsConstraints[ res + 1 ]; i++)
                    {
                        el = h_cut->Elements[i];
                        aux = h_cut->Coefficients[i] * n1;
                        if( ((aux>0&&d1<0)||(aux<0&&d1>0))&&(aux%d1!=0))
                        {
                            aux = (aux/d1) -1;
                        }
                        else
                        {
                            aux = aux/d1;
                        }
                        //aux = aux< 0 ? (aux/d1) - 1 : aux/d1;
                        value_tes += aux*h_cut->xAsterisc[el];
                    }
                    //rhs = h_cut->rightSide[ res ]* n1;
                    if( ((rhs>0&&d1<0)||(rhs<0&&d1>0))&&(rhs%d1!=0))
                    {
                        rhs = (rhs/d1) -1;
                    }
                    else
                    {
                        rhs = rhs/d1;
                    }

                    if(value_tes>rhs*precision)
                    {
                        violation = value_tes - (rhs*precision);
                        if(violation>violation_best)
                        {
                            violation_best = violation;
                            nBest=n1;
                            dBest=d1;
                        }
                    }

                }
                n1++;
            }
        }

        if(violation_best!=0)
        {
            h_solution->SConst[k] = res;
            h_solution->SMult[k] = nBest;
            h_solution->SPAux[k] = dBest;
        }
        else
        {
            h_solution->SConst[k] = -1;
            h_solution->SMult[k] = -1;
            h_solution->SPAux[k] = -1;
        }
        free(u);
    }
    //  }
    free(constraints);
}


Cut_gpu* initial_runCPU(Cut_gpu *h_cut, int maxDenominator, int precision, double timeconst)
{


    Cut_gpu* out_h_cut;
    int nRuns;
    int i, numberC = 0 ;
    int nC_initial = 0, cont_ini;
    for(i = 0; i<h_cut->numberConstrains; i++)
    {

        // if((h_cut->typeConstraints[i] == RES_RR) || (h_cut->typeConstraints[i] == LPC_RR) || (h_cut->typeConstraints[i] == RES_NR) || (h_cut->typeConstraints[i] == RES_MODE) || (h_cut->typeConstraints[i] == RES_CONFCL)  )
        if((h_cut->typeConstraints[i] == RES_RR) )
        {
            numberC++;
        }
    }
    //RES_RR,LPC_RR,RES_NR,RES_MODE
    nC_initial = h_cut->numberConstrains;
    cont_ini = h_cut->cont;

    nRuns =  numberC + (h_cut->numberConstrains - nC_initial);

    solutionGpu *h_solution_r1 = allocationStructSolution1(h_cut,nRuns);
    runCPUR1(h_cut, h_solution_r1,  precision, timeconst);
    int cont=0;
    for(i=0; i<nRuns; i++)
    {
        if(h_solution_r1->SConst[i]!=-1)
        {
            cont++;
        }
    }

    if(cont>0)
    {
        //printf("Number cuts generated in the phase 1: %d ", cont);
        //out_h_cut = createCutsStrongPhaseOne(h_cut, h_solution_r1, cont,precision,nRuns,nC_initial, cont_ini);
        out_h_cut = createCutsOfPhaseOne(h_cut, h_solution_r1, cont,precision,nRuns);

        free(h_solution_r1);
        free(h_cut);
        return out_h_cut;
    }
    else
    {
        // printf("No cuts generate\n");
        free(h_solution_r1);
        return h_cut;
    }



}


void runCPUR2(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int *setConstraint, int precision, int maxDenominator, int nRuns,double timeGPU)
{

    int k;

//    printf("TIME: %lf\n",timeGPU);
//    timeGPU = ((double) timeGPU - (omp_get_wtime() - timeCurrent ));
//    if(timeGPU<1)
//    {
//        for(k=0; k<nRuns; k++)
//        {
//            for(i=0; i<numberMaxConst; i++)
//            {
//                h_solution->SConst[i + k*numberMaxConst] = -1;
//            }
//            h_solution->SPAux[k] = 0;
//            h_solution->SMult[0 + k*4] = -1;
//            h_solution->SMult[1 + k*4] = -1;
//            h_solution->SMult[2 + k*4] = -1;
//            h_solution->SMult[3 + k*4] = -1;
//        }
//        free(Coef);
//        free(Coef2);
//        return ;
//
//    }

    //  #pragma omp parallel
    // {
    //  #pragma omp for
    for(k=0; k<nRuns; k++)
    {

        int mult_1, mult_2, rest_a,rest_b, i, j, el, rhs1, rhs2, value_tes, violation = 0, aux, n1_best = -1, n2_best = -1, d1_best = -1, qnt_1, d2_best=-1;//, cont=0;
        // double timeCurrent = omp_get_wtime();
        int Numerator[20];
        int Denominator[20];
        int *Coef = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
        int *Coef2 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));

        violation = 0;
        n1_best = -1;
        n2_best = -1;
        d1_best = -1;
        d2_best = -1;
        for(i=0; i<20; i++)
        {
            Denominator[i]= rand()%maxDenominator + 2;
            Numerator[i] = rand()%(Denominator[i]-1);
        }
        for(mult_1=0; mult_1<20; mult_1++)
        {
            memset(Coef,0,sizeof(int)*h_cut->numberVariables);
            rhs1 = 0;
            for(rest_a = 0; rest_a< numberMaxConst; rest_a++)
            {
                for(i=h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_a] ]; i<h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_a] + 1]; i++)
                {

                    el = h_cut->Elements[i];
                    Coef[el] += h_cut->Coefficients[i] * Numerator[mult_1];
                }
                rhs1 += h_cut->rightSide[ setConstraint[k*numberMaxConst+rest_a] ] * Numerator[mult_1];
                for(mult_2 = 0; mult_2<20; mult_2++)
                {
                    memset(Coef2,0,sizeof(int)*h_cut->numberVariables);
                    value_tes = 0;
                    rhs2 = 0;
                    for(rest_b = rest_a + 1; rest_b<numberMaxConst; rest_b++)
                    {
                        for(j=h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_b] ]; j<h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_b] + 1]; j++)
                        {
                            el = h_cut->Elements[j];
                            Coef2[el] += h_cut->Coefficients[j] * Numerator[mult_2];
                        }
                        rhs2 += h_cut->rightSide[ setConstraint[k*numberMaxConst + rest_b] ]* Numerator[mult_2];
                    }
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if( (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) < 0 )
                        {
                            aux = (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) - 1;

                        }
                        else
                        {
                            aux = (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]);
                        }
                        value_tes += aux*h_cut->xAsterisc[j];


//                        aux = Coef[j]<0 ? Coef[j]/Denominator[mult_1] - 1 : Coef[j]/Denominator[mult_1];
//                        value_tes += aux*h_cut->xAsterisc[j];
//                        aux = Coef2[j]<0 ? Coef2[j]/Denominator[mult_2] - 1 : Coef2[j]/Denominator[mult_2];
//                        value_tes += aux*h_cut->xAsterisc[j];
                    }

                    if( (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) < 0)
                    {
                        aux = (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) - 1;
                    }
                    else
                    {
                        aux = (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]);
                    }

//                    aux = rhs1<0 ? rhs1/Denominator[mult_1]-1 : rhs1/Denominator[mult_1];
//                    aux +=  rhs2<0 ? rhs2/Denominator[mult_2]-1 : rhs2/Denominator[mult_2];


                    if((value_tes>aux*precision)&&(value_tes-(aux*precision)>violation))
                    {
                        violation = value_tes-(aux*precision);
                        //printf("violation %d\n",violation);
                        n1_best = Numerator[mult_1];
                        d1_best = Denominator[mult_1];
                        n2_best = Numerator[mult_2];
                        d2_best = Denominator[mult_2];
                        qnt_1 = rest_a;
                    }


                }
            }

        }

        if(violation>0)
        {
            for(i=0; i<numberMaxConst; i++)
            {
                h_solution->SConst[i + k*numberMaxConst] = setConstraint[k*numberMaxConst + i];//CPU ja vai ter
            }

            h_solution->SPAux[k] = qnt_1;
            h_solution->SMult[0 + k*4] = n1_best;
            h_solution->SMult[1 + k*4] = d1_best;
            h_solution->SMult[2 + k*4] = n2_best;
            h_solution->SMult[3 + k*4] = d2_best;

        }
        else
        {
            for(i=0; i<numberMaxConst; i++)
            {
                h_solution->SConst[i + k*numberMaxConst] = -1;
            }
            h_solution->SPAux[k] = 0;
            h_solution->SMult[0 + k*4] = -1;
            h_solution->SMult[1 + k*4] = -1;
            h_solution->SMult[2 + k*4] = -1;
            h_solution->SMult[3 + k*4] = -1;
        }

        free(Coef);
        free(Coef2);
    }
    // }
}

int createSolutionsInitial(int *h_solution, int sz)
{
    struct timeval time;
    gettimeofday(&time,NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int i;
    for(i = 0; i<sz; i++)
    {
        h_solution[i] = rand()%2;
    }


}

int* calcCover(Cover_gpu *h_cover, int *h_solution, int qnt_Cover_per_Thread, int tolerance)
{
    int *fillBag = (int*)malloc(sizeof(int)*qnt_Cover_per_Thread*h_cover->numberConstraints);
    int i,j,k, counter = 0;
    for(i= 0; i < h_cover->numberConstraints; i++)
    {
        for(j = 0; j<qnt_Cover_per_Thread; j++)
        {
            counter = 0;
            for(k = h_cover->ElementsConstraints[i]; k< h_cover->ElementsConstraints[i+1]; k++)
            {
                counter += h_cover->Coefficients[k] * h_solution[k + j*h_cover->cont];
                //printf("%d %d %d\t", h_solution[k+j*h_cover->cont], h_cover->Coefficients[k], counter);
            }
            if(counter + tolerance <= h_cover->rightSide[i])
            {
                for(k = h_cover->ElementsConstraints[i]; k<h_cover->ElementsConstraints[i+1]; k++)
                {
                    if(h_solution[k+j*h_cover->cont]==0)
                    {
                        counter += h_cover->Coefficients[k];
                        h_solution[k+j*h_cover->cont] = 1;
                        if(counter + tolerance > h_cover->rightSide[i])
                        {
                            break;
                        }
                    }
                }
            }
            fillBag[j + i*qnt_Cover_per_Thread] = counter;

            //printf("Fill bag: %d Capacity bag: %d\n", fillBag[j + i*qnt_Cover_per_Thread], h_cover->rightSide[i]);
        }
        //getchar();
    }
    //printf("number Constraints: %d qnt %d \n", h_cover->numberConstraints, qnt_Cover_per_Thread);
    return fillBag;


}

Cut_gpu *runCPU_Cut_Cover(Cut_gpu *h_cut, int qnt_Cover_per_Thread, int nConstraintsInitial)
{

    //Cover_gpu *h_cover_new = AllocationStructCover(h_cover->cont*qnt_Cover_per_Thread,h_cover->numberConstraints*qnt_Cover_per_Thread);

    int i, j, k, w;
    double b = 0.0, bdc = 0.0, a_barra = 0.0;
    int counter, qnt, aux = 0;
    int *h_solution = (int*)malloc(sizeof(int)*h_cut->cont*qnt_Cover_per_Thread);
    createSolutionsInitial(h_solution,h_cut->cont*qnt_Cover_per_Thread);
    int *fillBag;
   // int nConstraintsInitial = h_cut->numberConstrains;
    Cover_gpu *h_cover = CopyCutToCover(h_cut,nConstraintsInitial);
    fillBag = calcCover(h_cover,h_solution,qnt_Cover_per_Thread,0);
    int contInitial = h_cover->cont;
    for(j=0; j<qnt_Cover_per_Thread; j++)
    {

        for(i=0; i< nConstraintsInitial; i++)
        {
            if(fillBag[j + i*qnt_Cover_per_Thread]<h_cut->rightSide[i]){
                continue;
            }
            qnt = 0;
            int el;
            for(k = h_cover->ElementsConstraints[i]; k < h_cover->ElementsConstraints[i+1]; k++)
            {
                    qnt+=h_solution[k + j*contInitial];
                //printf("%d ", h_solution[k + j*contInitial]);
            }
           //printf("\n");
            //printf("%d\n",qnt);
            //getchar();
            int *n_coef = (int*)malloc(sizeof(int)*qnt);
            int *n_el = (int*)malloc(sizeof(int)*qnt);
            qnt = 0;
            for(k = h_cover->ElementsConstraints[i]; k < h_cover->ElementsConstraints[i+1]; k++)
            {
                if(h_solution[k + j*contInitial]==1)
                {
                    n_coef[qnt] = h_cover->Coefficients[k];
                    n_el[qnt] = k;
                    qnt++;
                }
            }



            b = h_cover->rightSide[i];
            bdc = (float)h_cover->rightSide[i] / (float)qnt;
            int sz =  qnt;
            //int el = h_cover->ElementsConstraints[i];
            float delta = 0;
            float phi = (float)fillBag[j + i*qnt_Cover_per_Thread] - h_cover->rightSide[i];
            k = 1;
            //printf("fill bag %f \n", phi);
            a_barra = (double)n_coef[0];
            for( w = 1 ; w < qnt; w++ )
            {

                delta = a_barra - n_coef[w];
                if((float)k*delta<phi)
                {
                    a_barra = n_coef[w];
                    phi = phi - (float)k*delta;
                }
                else
                {
                    a_barra = a_barra - (phi/k);
                    phi = 0;
                    break;
                }
                k++;

            }
            if(phi>0)
            {
                a_barra = b/qnt;
            }
            //printf("a_barra = %f\n", a_barra);
            int *c_menus = (int*)malloc(sizeof(int)*qnt);
            int *c_mais = (int*)malloc(sizeof(int)*qnt);
            float *S_barra = (float*)malloc(sizeof(float)*(qnt+1) );
            int id1 = 0,id2 = 0, id3 = 0;
            for(w = 0; w < qnt; w++ )
            {   //printf("%d ", n_coef[w]);
                if((float)n_coef[w] <= a_barra)
                {
                    c_menus[id1] = w;
                    id1++;
                }
                else
                {
                    c_mais[id2] = w;
                    id2++;
                }
            }
            //printf("\n");
            S_barra[id3] = 0;
            id3++;
            for(w = 0; w<id2; w++)
            {
                S_barra[id3] = S_barra[id3-1] + a_barra;
                id3++;
            }
            for(w = 0; w<id1; w++)
            {
                S_barra[id3] = S_barra[id3-1] + (float)n_coef[ c_menus[w] ];
                id3++;
            }
            int ini = 0, fim = 0, meio = 0;
            //printf("S_barra: ");
           //for(w=0;w<=qnt;w++){
            //    printf("%f ", S_barra[w]);
            //}
            //printf("= %d",h_cover->rightSide[i]);
           // printf("\n");
            for(w = h_cover->ElementsConstraints[i]; w<h_cover->ElementsConstraints[i+1]; w++)
            {
                ini  = 0;
                fim  = id3 - 1;
                while(ini<=fim)
                {
                    meio = (ini + fim)/2;
                    if( (h_cover->Coefficients[w] - 1e-6 <= S_barra[meio])&&(h_cover->Coefficients[w] - 1e-6>S_barra[meio-1]) )
                    {
                        h_cover->Coefficients[w] = meio-1;
                        break;
                    }
                    else
                    {
                        if(h_cover->Coefficients[w]<S_barra[meio])
                        {
                            fim = meio - 1;
                        }
                        if(h_cover->Coefficients[w]>S_barra[meio])
                        {
                            ini = meio + 1;
                        }
                    }
                }
            }
            for(w=0; w<id1; w++)
            {
                el = n_el[ c_menus[w] ];
                h_cover->Coefficients[ el ] = 1;
            }
            //d_cover->rightSide[j] = qnt - 1;


            h_cover->rightSide[i] = qnt - 1;
            free(c_menus);
            free(c_mais);
            free(S_barra);
            free(n_coef);
            free(n_el);

        }
        int qnt_cuts_cover = 0;
        int *idc_cover = (int*)malloc(sizeof(int)*h_cover->numberConstraints);
        //int j;
        for(i=0; i<h_cover->numberConstraints; i++)
        {
            idc_cover[i] = 0;
            if(h_cover->rightSide[i] != h_cut->rightSide[i])
            {
                idc_cover[i] = 1;
                qnt_cuts_cover++;
            }
        }


        h_cut = createCutsCover(h_cut,h_cover,idc_cover,qnt_cuts_cover);
        free(h_cover);
        h_cover = CopyCutToCover(h_cut,nConstraintsInitial);
        free(idc_cover);

    }
    free(h_solution);
    free(fillBag);
    free(h_cover);
    return h_cut;

}



Cut_gpu* second_phase_runCPU(Cut_gpu *h_cut, int numberMaxConst, int nRuns, int maxDenominator, int precision, int szR, double timeLeft)
{

    //printf("FASE 2 in CPU\n"); fflush(stdout);
    int *consR1;
    int *consNR1;
    int *nElemR1;

    double startT = omp_get_wtime();
    double _time = 0;

    Cut_gpu* out_cut_gpu;

    int n_r = 0, n_nr = 0, i;//j;
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        return h_cut;
    }
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if((h_cut->typeConstraints[i]==RES_RR)||(h_cut->typeConstraints[i]==LPC_CGGPU))
        {
            n_r++;
        }
        else
        {
            n_nr++;
        }
    }
    consR1 = (int*)malloc(sizeof(int)*n_r);
    nElemR1 = (int*)malloc(sizeof(int)*n_r);
    consNR1 = (int*)malloc(sizeof(int)*n_nr);

    n_r = 0;
    n_nr = 0;
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consR1);
        free(nElemR1);
        free(consNR1);
        return h_cut;

    }

    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if((h_cut->typeConstraints[i]==RES_RR)||(h_cut->typeConstraints[i]==LPC_CGGPU))
        {
            consR1[n_r] = i;
            nElemR1[n_r] = h_cut->ElementsConstraints[i+1] - h_cut->ElementsConstraints[i];
            n_r++;
        }
        else
        {
            if(h_cut->typeConstraints[i]!=LPC_CGCPU)
            {
                consNR1[n_nr]=i;
                n_nr++;
            }
        }
    }
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consNR1);
        free(nElemR1);
        free(consR1);
        return h_cut;

    }

    bubble_sort(nElemR1,consR1,n_r);
    int *Similar = returnOrdConstrainsNR(h_cut);
    float *folga = returnFolga(h_cut);
    solutionGpu *h_solution_r2 = allocationStructSolution2(h_cut,numberMaxConst,nRuns);
    int pos_R1 = 0;
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consNR1);
        free(nElemR1);
        free(consR1);
        free(Similar);
        free(folga);
        free(h_solution_r2);
        return h_cut;

    }

    int *setConstraint = (int*)malloc(sizeof(int)*numberMaxConst*nRuns);
    calcSetConstraint(setConstraint, &pos_R1,numberMaxConst, h_cut->numberConstrains, consR1, consNR1, n_r, n_nr, Similar, folga,  nRuns, szR);
    shuffle_Set(setConstraint, numberMaxConst, numberMaxConst*nRuns);
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consNR1);
        free(nElemR1);
        free(consR1);
        free(Similar);
        free(folga);
        free(h_solution_r2);
        free(setConstraint);
        return h_cut;

    }



    runCPUR2(h_cut, h_solution_r2, numberMaxConst,setConstraint,precision,maxDenominator,nRuns,_time);
    int cont=0;
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {

        free(consNR1);
        free(nElemR1);
        free(consR1);
        free(Similar);
        free(folga);
        free(h_solution_r2);
        free(setConstraint);
        return h_cut;
    }
    for(i=0; i<nRuns; i++)
    {
        if(h_solution_r2->SConst[0 + i*numberMaxConst]!=-1)
        {
            cont++;
        }
    }
    if(cont>0)
    {
        //  printf("Number of Cuts in the second phase in CPU:%d %lf\n",cont,timeGPU);
        //out_cut_gpu = createCutsOfPhaseTwo(h_cut,cut_aux,h_solution_r2,numberMaxConst,cont,precision,nRuns);
        int nT = nRuns;
        int nB = 1;
        out_cut_gpu = createCutsStrongPhaseTwo(h_cut,h_solution_r2,numberMaxConst,cont,precision,nRuns,nT,nB);
    }
    else
    {
        free(consR1);
        free(consNR1);
        free(Similar);
        free(folga);
        free(nElemR1);
        free(setConstraint);
        free(h_solution_r2);

        return h_cut;
    }

    free(consR1);
    free(consNR1);
    free(Similar);
    free(folga);
    free(nElemR1);
    free(setConstraint);
    free(h_solution_r2);
    free(h_cut);
    return out_cut_gpu;
}

