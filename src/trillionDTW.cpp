// TrillionDTW codebase. See README for attribution.
#include "trillionDTW.h"
#include "deque.h"
#include "util.h"

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) ((x-y)*(x-y))

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(double *t, int len, int r, double *l, double *u)
{
    struct deque du(2*r + 2);
    struct deque dl(2*r + 2);

    du.push_back(0);
    dl.push_back(0);

    for (int i = 1; i < len; i++)
    {
        if (i > r)
        {
            u[i-r-1] = t[du.front()];
            l[i-r-1] = t[dl.front()];
        }
        if (t[i] > t[i-1])
        {
            du.pop_back();
            while (!du.empty() && t[i] > t[du.back()])
                du.pop_back();
        }
        else
        {
            dl.pop_back();
            while (!dl.empty() && t[i] < t[dl.back()])
                dl.pop_back();
        }
        du.push_back(i);
        dl.push_back(i);
        if (i == 2*r + 1 + du.front())
            du.pop_front();
        else if (i == 2*r + 1 + dl.front())
            dl.pop_front();
    }
    for (int i = len; i < len+r+1; i++)
    {
        u[i-r-1] = t[du.front()];
        l[i-r-1] = t[dl.front()];
        if (i-du.front() >= 2*r + 1)
            du.pop_front();
        if (i-dl.front() >= 2*r + 1)
            dl.pop_front();
    }
}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give siginifant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
double lb_kim_hierarchy(double *t, double *q, int j, int len, double mean, double std, double bsf)
{
    /// 1 point at front and back
    double d, lb;
    double x0 = (t[j] - mean) / std;
    double y0 = (t[(len-1+j)] - mean) / std;
    lb = dist(x0,q[0]) + dist(y0,q[len-1]);
    if (lb >= bsf)   return lb;

    /// 2 points at front
    double x1 = (t[(j+1)] - mean) / std;
    d = min(dist(x1,q[0]), dist(x0,q[1]));
    d = min(d, dist(x1,q[1]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 2 points at back
    double y1 = (t[(len-2+j)] - mean) / std;
    d = min(dist(y1,q[len-1]), dist(y0, q[len-2]) );
    d = min(d, dist(y1,q[len-2]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 3 points at front
    double x2 = (t[(j+2)] - mean) / std;
    d = min(dist(x0,q[2]), dist(x1, q[2]));
    d = min(d, dist(x2,q[2]));
    d = min(d, dist(x2,q[1]));
    d = min(d, dist(x2,q[0]));
    lb += d;
    if (lb >= bsf)   return lb;

    /// 3 points at back
    double y2 = (t[(len-3+j)] - mean) / std;
    d = min(dist(y0,q[len-3]), dist(y1, q[len-3]));
    d = min(d, dist(y2,q[len-3]));
    d = min(d, dist(y2,q[len-2]));
    d = min(d, dist(y2,q[len-1]));
    lb += d;

    return lb;
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(int* order, double *t, double *uo, double *lo, double *cb, int j, int len, double mean, double std, double best_so_far)
{
    double lb = 0;
    double x, d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        x = (t[(order[i]+j)] - mean) / std;
        d = 0;
        if (x > uo[i])
            d = dist(x,uo[i]);
        else if(x < lo[i])
            d = dist(x,lo[i]);
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelop of the current data
double lb_keogh_data_cumulative(int* order, double *tz, double *qo, double *cb, double *l, double *u, int len, double mean, double std, double best_so_far)
{
    double lb = 0;
    double uu,ll,d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        uu = (u[order[i]]-mean)/std;
        ll = (l[order[i]]-mean)/std;
        d = 0;
        if (qo[i] > uu)
            d = dist(qo[i], uu);
        else
        {   if(qo[i] < ll)
            d = dist(qo[i], ll);
        }
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* A, double* B, double *cb, int m, int r, double bsf)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=dist(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(A[i],B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i+r < m-1 && min_cost + cb[i+r+1] >= bsf)
        {   free(cost);
            free(cost_prev);
            return min_cost + cb[i+r+1];
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}

/// Main Calculation Function
int calculate(const char *dataPath, const char *queryPath, int queryLength, int r)
{
    FILE *fp;            /// data file pointer
    FILE *qp;            /// query file pointer
    double bsf;          /// best-so-far
    double *t, *q;       /// data array and query array
    int *order;          ///new order of the query
    double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;


    double d;
    long long i , j;
    double ex , ex2 , mean, std;
    int m=-1;
    long long loc = 0;
    double t1,t2;
    int kim = 0,keogh = 0, keogh2 = 0;
    double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
    double *buffer, *u_buff, *l_buff;
    Index *Q_tmp;

    int EPOCH = 100000;
    m=queryLength;

    t1 = clock();
    fp = fopen(dataPath,"r");
    if( fp == NULL )
        error(2);

    qp = fopen(queryPath,"r");
    if( qp == NULL )
        error(2);

    /// malloc everything here
    q = (double *)malloc(sizeof(double)*m);
    if( q == NULL )
        error(1);
    qo = (double *)malloc(sizeof(double)*m);
    if( qo == NULL )
        error(1);
    uo = (double *)malloc(sizeof(double)*m);
    if( uo == NULL )
        error(1);
    lo = (double *)malloc(sizeof(double)*m);
    if( lo == NULL )
        error(1);

    order = (int *)malloc(sizeof(int)*m);
    if( order == NULL )
        error(1);

    Q_tmp = (Index *)malloc(sizeof(Index)*m);
    if( Q_tmp == NULL )
        error(1);

    u = (double *)malloc(sizeof(double)*m);
    if( u == NULL )
        error(1);

    l = (double *)malloc(sizeof(double)*m);
    if( l == NULL )
        error(1);

    cb = (double *)malloc(sizeof(double)*m);
    if( cb == NULL )
        error(1);

    cb1 = (double *)malloc(sizeof(double)*m);
    if( cb1 == NULL )
        error(1);

    cb2 = (double *)malloc(sizeof(double)*m);
    if( cb2 == NULL )
        error(1);

    u_d = (double *)malloc(sizeof(double)*m);
    if( u == NULL )
        error(1);

    l_d = (double *)malloc(sizeof(double)*m);
    if( l == NULL )
        error(1);

    t = (double *)malloc(sizeof(double)*m*2);
    if( t == NULL )
        error(1);

    tz = (double *)malloc(sizeof(double)*m);
    if( tz == NULL )
        error(1);

    buffer = (double *)malloc(sizeof(double)*EPOCH);
    if( buffer == NULL )
        error(1);

    u_buff = (double *)malloc(sizeof(double)*EPOCH);
    if( u_buff == NULL )
        error(1);

    l_buff = (double *)malloc(sizeof(double)*EPOCH);
    if( l_buff == NULL )
        error(1);

    /// Read query file
    bsf = INF;
    i = 0;
    j = 0;
    ex = ex2 = 0;

    while(fscanf(qp,"%lf",&d) != EOF && i < m)
    {
        ex += d;
        ex2 += d*d;
        q[i] = d;
        i++;
    }
    fclose(qp);

    /// Do z-normalize the query, keep in same array, q
    mean = ex/m;
    std = ex2/m;
    std = sqrt(std-mean*mean);
    for( i = 0 ; i < m ; i++ )
         q[i] = (q[i] - mean)/std;

    /// Create envelop of the query: lower envelop, l, and upper envelop, u
    cout<<"Envelop created"<<endl;
    lower_upper_lemire(q, m, r, l, u);

    /// Sort the query one time by abs(z-norm(q[i]))
    for( i = 0; i<m; i++)
    {
        Q_tmp[i].value = q[i];
        Q_tmp[i].index = i;
    }
    qsort(Q_tmp, m, sizeof(Index), Index::comp);

    /// also create another arrays for keeping sorted envelop
    for( i=0; i<m; i++)
    {   int o = Q_tmp[i].index;
        order[i] = o;
        qo[i] = q[o];
        uo[i] = u[o];
        lo[i] = l[o];
    }
    free(Q_tmp);

    /// Initial the cummulative lower bound
    for( i=0; i<m; i++)
    {   cb[i]=0;
        cb1[i]=0;
        cb2[i]=0;
    }

    i = 0;          /// current index of the data in current chunk of size EPOCH
    j = 0;          /// the starting index of the data in the circular array, t
    ex = ex2 = 0;
    bool done = false;
    int it=0, ep=0, k=0;
    long long I;    /// the starting index of the data in current chunk of size EPOCH

    while(!done)
    {
        /// Read first m-1 points
        ep=0;
        if (it==0)
        {   for(k=0; k<m-1; k++)
                if (fscanf(fp,"%lf",&d) != EOF)
                    buffer[k] = d;
        }
        else
        {   for(k=0; k<m-1; k++)
                buffer[k] = buffer[EPOCH-m+1+k];
        }

        /// Read buffer of size EPOCH or when all data has been read.
        ep=m-1;
        while(ep<EPOCH)
        {   if (fscanf(fp,"%lf",&d) == EOF)
                break;
            buffer[ep] = d;
            ep++;
        }

        /// Data are read in chunk of size EPOCH.
        /// When there is nothing to read, the loop is end.
        if (ep<=m-1)
        {   done = true;
        } else
        {   lower_upper_lemire(buffer, ep, r, l_buff, u_buff);

            /// Just for printing a dot for approximate a million point. Not much accurate.
            if (it%(1000000/(EPOCH-m+1))==0)
                fprintf(stderr,".");

            /// Do main task here..
            ex=0; ex2=0;
            for(i=0; i<ep; i++)
            {
                /// A bunch of data has been read and pick one of them at a time to use
                d = buffer[i];

                /// Calcualte sum and sum square
                ex += d;
                ex2 += d*d;

                /// t is a circular array for keeping current data
                t[i%m] = d;

                /// Double the size for avoiding using modulo "%" operator
                t[(i%m)+m] = d;

                /// Start the task when there are more than m-1 points in the current chunk
                if( i >= m-1 )
                {
                    mean = ex/m;
                    std = ex2/m;
                    std = sqrt(std-mean*mean);

                    /// compute the start location of the data in the current circular array, t
                    j = (i+1)%m;
                    /// the start location of the data in the current chunk
                    I = i-(m-1);

                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_kim = lb_kim_hierarchy(t, q, j, m, mean, std, bsf);

                 //   cout<<"BSF "<<bsf<<endl;
                    if (lb_kim < bsf)
                    {
                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                        /// uo, lo are envelop of the query.
                        lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, bsf);
                        if (lb_k < bsf)
                        {
                            /// Take another linear time to compute z_normalization of t.
                            /// Note that for better optimization, this can merge to the previous function.
                            for(k=0;k<m;k++)
                            {   tz[k] = (t[(k+j)] - mean)/std;
                            }

                            /// Use another lb_keogh to prune
                            /// qo is the sorted query. tz is unsorted z_normalized data.
                            /// l_buff, u_buff are big envelop for all data in this chunk
                            lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff+I, u_buff+I, m, mean, std, bsf);
                            if (lb_k2 < bsf)
                            {
                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                /// Note that cb and cb2 will be cumulative summed here.
                                if (lb_k > lb_k2)
                                {
                                    cb[m-1]=cb1[m-1];
                                    for(k=m-2; k>=0; k--)
                                        cb[k] = cb[k+1]+cb1[k];
                                }
                                else
                                {
                                    cb[m-1]=cb2[m-1];
                                    for(k=m-2; k>=0; k--)
                                        cb[k] = cb[k+1]+cb2[k];
                                }

                                /// Compute DTW and early abandoning if possible
                                dist = dtw(tz, q, cb, m, r, bsf);

                                if( dist < bsf )
                                {   /// Update bsf
                                    /// loc is the real starting location of the nearest neighbor in the file
                                    bsf = dist;
                                    loc = (it)*(EPOCH-m+1) + i-m+1;
                                }
                            } else
                                keogh2++;
                        } else
                            keogh++;
                    } else
                        kim++;

                    /// Reduce obsolute points from sum and sum square
                    ex -= t[j];
                    ex2 -= t[j]*t[j];
                }
            }

            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            if (ep<EPOCH)
                done=true;
            else
                it++;
        }
    }

    i = (it)*(EPOCH-m+1) + ep;
    fclose(fp);

    free(q);
    free(u);
    free(l);
    free(uo);
    free(lo);
    free(qo);
    free(cb);
    free(cb1);
    free(cb2);
    free(tz);
    free(t);
    free(l_d);
    free(u_d);
    free(l_buff);
    free(u_buff);

    t2 = clock();
    printf("\n");

    /// Note that loc and i are long long.
    cout << "Location : " << loc << endl;
    cout << "Distance : " << sqrt(bsf) << endl;
    cout << "Data Scanned : " << i << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;

    /// printf is just easier for formating ;)
    printf("\n");
    printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
    printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
    printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
    printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/i*100));
    return 0;
}
