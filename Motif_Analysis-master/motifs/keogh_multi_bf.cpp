/******************************************************************************
 *******************************************************************************
 ******                                                                  *******
 ******     This code is written by Abdullah Al Mueen at the department  *******
 ******     of Computer Science and Engineering of University of         *******
 ******     California - RIverside.                                      *******
 ******                                                                  *******
 *******************************************************************************
 ******************************************************************************/

/*#############################################################################
 * ######                                                                  #######
 * ######     This code is open to use, reuse and distribute at the user's #######
 * ######     own risk and as long as the above credit is ON. The complete #######
 * ######     description of the algorithm and methods applied can be      #######
 * ######     found in the paper - EXACT DISCOVERY OF TIME SERIES MOTIFS   #######
 * ######     by Abdullah Mueen, Eamonn Keogh and Qiang Zhu.               #######
 * ######                                                                  #######
 * #############################################################################*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <iostream>
#include <string.h>
#include "mex.h"

#define INF 999999999999.0                                                      //Pseudo Infitinte number for this code

FILE *fp;                                                                       //The input file pointer

double **data;
long long loc1 , loc2;

long long TIMESERIES;
int LENGTH;
int dim;
double bsf;
double total = 0.0;
int abandon = 1;

/* Calculates the distance between two time series x and y. If the distance is
 * larger than the best so far (bsf) it stops computing and returns the approximate
 * distance. To get exact distance the bsf argument should be omitted.*/

double distance(double *x, double *y, int length , double best_so_far = INF ) {
    int i;
    double sum = 0;
    
    if( abandon == 1 ) {
        double bsf2 = best_so_far*best_so_far;
        for ( i = 0 ; i < length && sum < bsf2  ; i++ )
            sum += (x[i]-y[i])*(x[i]-y[i]);
        return sqrt(sum);
    }
    else {
        for ( i = 0 ; i < length  ; i++ )
            sum += (x[i]-y[i])*(x[i]-y[i]);
        return sqrt(sum);
    }
    
}

/*Comparison function for qsort function. Compares two time series by using their
 * distances from the reference time series. */

void error(int id) {
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
        printf("ERROR : Invalid Number of Arguments!!!\n\n");
    
    exit(1);
    
}



void stop_exec(int sig) {
    printf("Current Motif is (%lld", loc1);
    printf(",%lld)", loc2);
    printf(" and the distance is %lf\n", bsf);
    exit(1);
}

/*arguments passed in plhs should be as follows:
 * 1 - the raw time series
 * 2 - length of the time series
 * 3 - number of dimensions in time series
 * 4 - length of the motif to search for
 * 5 - window to exclude to prevent motif overlap.  Should normally be set to motif length.
 * 6 - whether to print out all output
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double d;
    long long i , j , n;
    int  r = 0;
    double ex , ex2 , mean, std;
    long long length;
    int clear = 1;
    double t1, t2;
    int verbos = 0;
    
    bsf = INF;
    i = 0;
    j = 0;
    ex = ex2 = 0;
    
    /* Check for proper number of input and output arguments */
    if (nrhs != 6) {
        mexErrMsgTxt("Six input arguments required.");
    }
    if (nlhs != 3){
        mexErrMsgTxt("Three output arguments required.");
    }
    
    t1 = clock();
    signal(SIGINT, stop_exec);
    
    
    //fp = fopen(filename,"r");
    double* temp;
    double* rawTS = mxGetPr(prhs[0]);
    temp = mxGetPr(prhs[1]);
    length = (int)temp[0];
    temp = mxGetPr(prhs[2]);
    dim = (int)temp[0];
    temp = mxGetPr(prhs[3]);
    LENGTH = temp[0];
    temp = mxGetPr(prhs[4]);
    clear = temp[0];
    temp = mxGetPr(prhs[5]);
    verbos = temp[0];
    
    if( verbos == 1 )
        printf("\nLength of the Time Series : %lld\nLength of Subsequences : %d\n\n", length, LENGTH);
    
    TIMESERIES = length-LENGTH+1;
    
    data = (double **)mxMalloc(sizeof(double *)*TIMESERIES);
    if( data == NULL )
        error(1);
    //no longer need to read this in...but we do have to fill data
    for(i=0;i<length;i++) {
        if( i < TIMESERIES ) {
            data[i] = (double *)mxMalloc(sizeof(double)*LENGTH*dim);
            if( data[i] == NULL )
                error(1);
        }
        
        //add each dimension after the other
        for(n = 0; n < dim; n++) {
            
            d = rawTS[i + n*length];
            
            ex += d;
            ex2 += d*d;
            
            //make an array of all subsequences of length of motif
            int k = LENGTH-1;
            for(j=((i>=TIMESERIES)?TIMESERIES-1:i);j>=((i>=LENGTH)?(i-k):0);j--)
                data[j][i-j + n*LENGTH] = d;
            /*
             * if( i >= k )
             * {
             * mean = ex/LENGTH;
             * std = ex2/LENGTH;
             * std = sqrt(std-mean*mean);
             * ex -= data[i-k][0];
             * ex2 -= data[i-k][0]*data[i-k][0];
             * for( int u = 0 ; u < LENGTH ; u++ )
             * data[i-k][u] = (data[i-k][u]-mean)/std;
             * }
             */
        }
    }
    
    
    
    if(verbos == 1)
        printf("Data Have been Read but not Normalized\n\n");
    
    for( i = 0 ; i < TIMESERIES; i++ )
        for( j = i + clear; j < TIMESERIES ; j++ ) {
        double x = distance(data[i], data[j], LENGTH*dim, bsf);
        if( abandon == 0 )
            total += x;
        //if( x == 0 ) continue;
        if( bsf > x ) {
            bsf = x;
            loc1 = i;
            loc2 = j;
            if(verbos == 1)
                printf("New best-so-far is %lf and (%lld , %lld) are the new motif pair\n", bsf, loc1, loc2);
        }
        }
    
    //Convert to matlab indexing
    loc1 = loc1 + 1;
    loc2 = loc2 + 1;
    
    double* curOut = mxGetPr(plhs[0]);
    curOut[0] = loc1;
    curOut = mxGetPr(plhs[1]);
    curOut[0] = loc2;
    curOut = mxGetPr(plhs[2]);
    curOut[0] = bsf;
    if(verbos == 1) {
        printf("\n\nFinal Motif is the pair ( %lld", loc1);
        printf(", %lld ) and the Motif distance is %lf\n", loc2, bsf);
    }
    t2 = clock();
    if(verbos == 1)
        printf("\nExecution Time was : %lf seconds\n", (t2-t1)/CLOCKS_PER_SEC);
    
    //free memory
    mxFree(data);
    
}

