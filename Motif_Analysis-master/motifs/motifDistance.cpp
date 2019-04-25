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


/*arguments passed in plhs should be as follows:
 * 1 - the raw time series
 * 2 - length of the time series
 * 3 - number of dimensions in time series
 * 4 - the raw motif to be compared to the time series
 * 5 - length of the motif to search for
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double d;
    long long i , j , n;
    long long length;
    int verbos = 0;
    
    bsf = INF;
    i = 0;
    j = 0;
    
    /* Check for proper number of input and output arguments */
    if (nrhs != 5) {
        mexErrMsgTxt("Five input arguments required.");
    }
    if (nlhs != 1){
        mexErrMsgTxt("One output argument required.");
    }
    
    
    //fp = fopen(filename,"r");
    double* temp;
    double* rawTS = mxGetPr(prhs[0]);
    temp = mxGetPr(prhs[1]);
    length = (int)temp[0];
    temp = mxGetPr(prhs[2]);
    dim = (int)temp[0];
    double* motifData = mxGetPr(prhs[3]);
    temp = mxGetPr(prhs[4]);
    LENGTH = temp[0];
    
    TIMESERIES = length-LENGTH+1;
    
    data = (double **)mxMalloc(sizeof(double *)*TIMESERIES);
    if( data == NULL )
        error(1);
    //fill in time series data
    for(i=0;i<length;i++) {
        if( i < TIMESERIES ) {
            data[i] = (double *)mxMalloc(sizeof(double)*LENGTH*dim);
            if( data[i] == NULL )
                error(1);
        }
        
        //add each dimension after the other
        for(n = 0; n < dim; n++) {
            
            d = rawTS[i + n*length];
            
            //make an array of all subsequences of length of motif
            int k = LENGTH-1;
            for(j=((i>=TIMESERIES)?TIMESERIES-1:i);j>=((i>=LENGTH)?(i-k):0);j--)
                data[j][i-j + n*LENGTH] = d;
        }
    }
    
    for( i = 0 ; i < TIMESERIES; i++ ) {
        double x = distance(data[i], motifData, LENGTH*dim, bsf);
        if( abandon == 0 )
            total += x;
        if( bsf > x )
            bsf = x;
    }
    
    double* curOut = mxGetPr(plhs[0]);
    curOut[0] = bsf;
    
    //free memory
    mxFree(data);
    
}

