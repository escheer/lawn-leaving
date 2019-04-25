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
double **dref;
double **indices;
long long *ind;                                                                 //The sorted indices (I-ordering) to the time series in order of their distances to ref
int *rInd;                                                                      //The sorted indices (Z-ordering) to the reference points in order of their std of distances

int ref;                                                                        //reference point that has maximum std
double *stdRef;                                                                 //Standard Deviations of all the reference points

long long TIMESERIES;                                                           //The length of the Long Time Series.
int LENGTH;                                                                     //Length of each time series
int MAXREF;                                                                     //Number of reference points to be used
int K;
double X;

long long loc1 , loc2;                                                          //location of the current best pair
double bsf;                                                                     //best-so-far




/* Calculates the distance between two time series x and y. If the distance is
 * larger than the best so far (bsf) it stops computing and returns the approximate
 * distance. To get exact distance the bsf argument should be omitted.*/

double distance(double *x, double *y, int length , double best_so_far = INF ) {
    int i;
    double sum = 0;
    double bsf2 = best_so_far*best_so_far;
    for ( i = 0 ; i < length && sum < bsf2 ; i++ )
        sum += (x[i]-y[i])*(x[i]-y[i]);
    return sqrt(sum);
}

/*Comparison function for qsort function. Compares two time series by using their
 * distances from the reference time series. */

int comp1(const void *a, const void *b) {
    long long *x=(long long *)a;
    long long *y=(long long *)b;
    
    if (indices[ref][*x]>indices[ref][*y])
        return 1;
    else if (indices[ref][*x]<indices[ref][*y])
        return -1;
    else
        return 0;
}


int comp2(const void *a, const void *b) {
    int *x=(int *)a;
    int *y=(int *)b;
    
    if (stdRef[*x]<stdRef[*y])
        return 1;
    else if (stdRef[*x]>stdRef[*y])
        return -1;
    else
        return 0;
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



void stop_exec(int sig) {
    printf("Current Motif is (%lld", loc1);
    printf(",%lld)", loc2);
    printf(" and the distance is %lf\n", bsf);
    exit(1);
}

/*arguments passed in plhs should be as follows:
 * 1 - the raw time series
 * 2 - a matched length string of where allowable starting positions are
 * 3 - length of the time series
 * 4 - dimension of data
 * 5 - length of the subsequences
 * 6 - # allowable reference points
 * 7 - window to exclude to prevent motif overlap.  Should normally be set to motif length.
 * 8 - whether to print out all output
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double d;
    long long i , j , n;
    long long offset = 0;
    int abandon = 0 ,  r = 0;
    double ex , ex2 , mean, std;
    long long length;
    int clear = 0;
    int dim = 1;
    double *cntr;
    double t1, t2;
    char fname[500];
    int verbos = 0;
    double count = 0;
    
    bsf = INF;
    i = 0;
    j = 0;
    ex = ex2 = 0;
    
    /* Check for proper number of input and output arguments */
    if (nrhs != 8) {
        mexErrMsgTxt("Eight input arguments required.");
    }
    if (nlhs != 3){
        mexErrMsgTxt("Three output arguments required.");
    }
    
    t1 = clock();
    signal(SIGINT, stop_exec);
    
    
    //fp = fopen(filename,"r");
    double* temp;
    double* rawTS = mxGetPr(prhs[0]);
    //indicates which are allowed starting points
    double* rawAllowed = mxGetPr(prhs[1]);
    temp = mxGetPr(prhs[2]);
    length = temp[0];
    temp = mxGetPr(prhs[3]);
    dim = (int)temp[0];
    temp = mxGetPr(prhs[4]);
    LENGTH = temp[0];
    temp = mxGetPr(prhs[5]);
    MAXREF = temp[0];
    temp = mxGetPr(prhs[6]);
    clear = temp[0];
    temp = mxGetPr(prhs[7]);
    verbos = temp[0];
    
    if( verbos == 1 )
        printf("\nLength of the Time Series : %lld\nLength of Subsequences : %d\n\n", length, LENGTH);
    
    TIMESERIES = length-LENGTH+1;
    
    data = (double **)mxMalloc(sizeof(double *)*TIMESERIES);
    ind = (long long *)mxMalloc(sizeof(long long)*TIMESERIES);
    if( data == NULL || ind == NULL )
        error(1);
    //no longer need to read this in...but we do have to fill data
    for(i=0;i<length;i++) {
        if( i < TIMESERIES ) {
            data[i] = (double *)mxMalloc(sizeof(double)*LENGTH*dim);
            if( data[i] == NULL )
                error(1);
            ind[i] = i;
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
    
    dref = (double **)mxMalloc(MAXREF*sizeof(double *));
    indices = (double **)mxMalloc(MAXREF*sizeof(double *));
    stdRef = (double *)mxMalloc(MAXREF*sizeof(double));
    cntr = (double *)mxMalloc(MAXREF*sizeof(double));
    rInd = (int *)mxMalloc(MAXREF*sizeof(int));
    
    if ( dref == NULL || indices == NULL || stdRef == NULL || cntr == NULL || rInd == NULL )
        error(1);
    //////////////////////////////////////////////////////////////////////////
    
    
    /*Generating the reference time series. Here it is a random time series*/
    
    srand( time(NULL) );
    
    for( r = 0 ; r < MAXREF ; r++ ) {
        
        dref[r] = (double *)mxMalloc(sizeof(double)*LENGTH*dim);
        indices[r] = (double *)mxMalloc(sizeof(double)*TIMESERIES);
        if( dref[r] == NULL || indices[r] == NULL )
            error(1);
        
        long long random_pick = rand() % TIMESERIES;
        for( i = 0 ; i < LENGTH*dim ; i++ )
            dref[r][i] = data[random_pick][i];
        
        if( verbos == 1 )
            printf("\nThe %lldth subsequence is chosen as %dth reference\n", random_pick, r);
        
        
        /*Creating the indices array which is a 2 by TIMESERIES
         * sized vector having the indices (to the data array) and distances (from the
         * reference time series) in each row.*/
        
        
        ex = 0;
        ex2 = 0;
        
        
        for( i = 0 ; i < TIMESERIES ; i++ ) {
            if( i == random_pick || rawAllowed[i]!=1 )
            { indices[r][i] = INF; continue; }
            d = indices[r][i] = distance(data[i], dref[r], LENGTH*dim);
            count = count + 1;
            ex += d;
            ex2 += d*d;
            if ( abs((double)(i - random_pick) ) < clear )  continue;
            if ( d < bsf ) {
                bsf = d; loc1 = i; loc2 = random_pick;
                if(verbos == 1)
                    printf("New best-so-far is %lf and (%lld , %lld) are the new motif pair\n", bsf, loc1, loc2);
            }
            
        }
        
        ex = ex/(TIMESERIES-1);
        ex2 = ex2/(TIMESERIES-1);
        std = sqrt(ex2-ex*ex);
        
        
        
        rInd[r] = r;
        stdRef[r] = std;
        cntr[r] = 0;
        //  printf("%d mean is %lf std is %lf \n",r,ex,std);
        ////////////////////////////////////////////////////////////////////
    }
    
    
    if(verbos == 1)
        printf("\nReferences are picked and Dist has been Computed\n\n");
    
    /*Sort the standard Deviations*/
    qsort(rInd, MAXREF, sizeof(int), comp2);
    
    ref = rInd[0];
    
    long long remaining = TIMESERIES;
    
    
    /*Sort indices using the distances*/
    qsort(ind, TIMESERIES, sizeof(long long), comp1);
    
    
    if(verbos == 1)
        printf("Orderings have been Computed and Search has begun\n\n");
    
    /*Motif Search loop of the algorithm that finds the motif. The algorithm
     * computes the distances between a pair of time series only when it thinks
     * them as a potential motif'*/
    
    
    //  printf("Loop starts\n");
    offset = 0;
    abandon = 0;
    
    while (!abandon && offset < remaining) {
        abandon = 1;
        offset++;
        
        for(i = 0 ; i < remaining - offset ; i++ ) {
            long long left = ind[i];
            long long right = ind[i + offset];
            if( abs((double)(left-right)) < clear )
                continue;
            
            //According to triangular inequality distances between left and right
            //is obviously greater than lower_bound.
            double lower_bound = 0;
            r = 0;
            do {
                lower_bound = fabs(indices[rInd[r]][right] - indices[rInd[r]][left]);
                r++;
            }while( r < MAXREF && lower_bound < bsf );
            
            
            if (r >= MAXREF && lower_bound < bsf && rawAllowed[left]==1 &&rawAllowed[right]==1) {
                
                abandon = 0;
                count =  count + 1;
                d = distance( data[left] , data[right] , LENGTH*dim, bsf );
                signal(SIGINT, SIG_IGN);
                if (d < bsf ) {
                    
                    bsf = d;
                    loc1 = left;
                    loc2 = right;
                    
                    if(verbos == 1)
                        printf("New best-so-far is %lf and (%lld , %lld) are the new motif pair\n", bsf, loc1, loc2);
                    
                }
                signal(SIGINT, stop_exec);
            }
        }
    }
    
    //convert to Matlab indexing
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
    mxFree(ind);
    mxFree(dref);
    mxFree(indices);
    mxFree(stdRef);
    mxFree(cntr);
    mxFree(rInd);
    
}

