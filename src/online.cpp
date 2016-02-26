#include "online.h"
#include "deque.h"
#include "util.h"
#include "trillionDTW.h"
#include "TimeSeries.h"
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <map>
#include <set>
#include <algorithm>

#define INF 1e20       //Pseudo Infitinte number for this code

//#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code

using namespace std;

//constructor passing the query data length
//ST, Ndata number of time series in data file, L is length of time series
onlineOperation::onlineOperation(double STnew, int Ndata, int Ldata)
{
   ST=STnew;        //new ST`
   N=Ndata;
   L=Ldata;
}

//type =true means to normalize, for ONEX, type =false mean run Spring, do not normalize in spring
int onlineOperation::readTimeSeries(const char* Filename, bool type)
{
    inputFile.open(Filename); //Opening the file
    if(!inputFile)
    {
        cout << "there's something wrong with the file!\n"; //Throw exception
        return -1;
    }
    else
    {
        //read the input file into time series structure
        //Read numbers from file into the array
        for(int i = 0; i < N; i++)
        {
            for(int j=0; j<L;j++)
            {
                //normalize at the same time
                inputFile>>timeSeries[i][j];
                if(type==true)
                    timeSeries[i][j]=(timeSeries[i][j]-min)/(max-min);

            }
        }
        inputFile.close();

        return 1;
    }

}

//reads the query file
//type =true means normalize and false means running for Spring don't normalize
void onlineOperation::readQueryFile(int m, const char * Filename, bool type)
{
    // Array for keeping the query data
    Q = (double *)malloc(sizeof(double)*m);
    queryLength=m;
    queryFile.open(Filename);
    //read the query data and normalize at the same time
    for(int i=0;i<m;i++)
    {
        queryFile>>Q[i];
        if(type==true)
            Q[i]=(Q[i]-min)/(max-min);
        tempQ.push_back(Q[i]);
    }
    queryFile.close();
}
bool sortByDist(const kBest &lhs, const kBest &rhs) { return lhs.dist < rhs.dist; }
//finds k similar time series, query Flag=true query sequence is given
//timeFlag=true time interval Exact, false means Any

void onlineOperation::kSimilar(bool queryFlag, int k,bool timeFlag,int startTime, int endTime)
{

    double bsf=INF;

    int bsfIndex;       //index of bsf centroid
    int bsfTSIndex;
    double currentDist=0;
  /*  if(queryFlag==true)
    {
        //read the query file
        readQueryFile("data/query.txt",queryLength);
    }
*/
    if(timeFlag==false)    //ANY time interval
    {
        //call that embedded code

        //compare the query with all the centroids
        for(unsigned int i=0;i<groupArray.size();i++)
        {
            currentDist=simpleDTW(groupArray[i].centroid,tempQ);
            if(currentDist<bsf)
            {
                bsf=currentDist;
                bsfIndex=i;
            }
        }
    }
    else
    {
        //EXACT match
        //iterate through the groups in the time interval provided
        int numGroups=Time[startTime][endTime].count;
        for(int i=0;i<numGroups;i++)
        {
            int groupID=Time[startTime][endTime].groupid[i];
            currentDist=simpleDTW(tempQ,groupArray[groupID].centroid);
            if(currentDist<bsf)
            {
                bsf=currentDist;
                bsfIndex=groupID;
            }
        }

    }
    //now we have the closest centroid
    bsf=INF;
    vector<double> tempSeries;
    bsfTSIndex=0;
    int kbestCount=0;
    kBest tempBest;

    if(k==1)
    {
        //compare the DTW between all the time series in this group
        for(int i=0;i<N;i++)
        {
            if(groupArray[bsfIndex].bitVector[i]==1)    //this time series is present in the group
            {
                //copy this time series in tempSeries
                for(int j=0;j<groupArray[bsfIndex].endT-groupArray[bsfIndex].startT+1;j++)
                {
                    tempSeries.push_back(timeSeries[i][groupArray[bsfIndex].startT+j]);
                }
                currentDist=simpleDTW(tempSeries,tempQ);
                if(currentDist<bsf)
                {
                    //update the best index and bsf
                    bsf=currentDist;
                    bsfTSIndex=i;
                }
                tempSeries.clear();
            }
        }
        //have the most similar time series
        cout<<"Most similar TS is "<<bsfTSIndex<<" Interval "<<groupArray[bsfIndex].startT<<" "<<groupArray[bsfIndex].endT<<" "<<bsf<<endl;
    }
    else
    {
        //need to find k most similar from this group
        if(groupArray[bsfIndex].count==k)   //
        {
            cout<<"K most similar TS are";
            for(int i=0;i<N;i++)
            {
                if(groupArray[bsfIndex].bitVector[i]==1)    //this time series is present in the group
                {
                    cout<<i<<" ";
                }
            }
        }
        else if(groupArray[bsfIndex].count>k)       //
        {
            for(int i=0;i<N;i++)
            {
                if(groupArray[bsfIndex].bitVector[i]==1)    //this time series is present in the group
                {
                    //copy this time series in tempSeries
                    for(int j=0;j<groupArray[bsfIndex].endT-groupArray[bsfIndex].startT+1;j++)
                    {
                        tempSeries.push_back(timeSeries[i][groupArray[bsfIndex].startT+j]);
                    }
                    currentDist=simpleDTW(tempQ,tempSeries);
                    if(kbestCount<k)
                    {
                        //add this TS to k best
                        kbestCount++;
                        tempBest.dist=currentDist;
                        tempBest.id=i;
                        kbestArray.push_back(tempBest);

                    }
                    else
                    {
                        sort(kbestArray.begin(),kbestArray.end(),sortByDist);
                        double tempD=kbestArray[kbestCount].dist;    //getting the last distance
                        if(tempD>currentDist)
                        {
                            tempBest.dist=currentDist;
                            tempBest.id=i;
                            kbestArray[kbestCount]=tempBest;

                        }
                    }
                    tempSeries.clear();
                }
            }
            cout<<"K most similar are";
            for(int i=0;i<k;i++)
                cout<<kbestArray[i].id<<" ";
            cout<<endl;

        }
        else
        {
            //count of TS in group is < k
            //finds second best group
        }
    }



}

//finds dominant or outlier TS
//For outlier the time interval is specified
int onlineOperation::dominantOutlier(int startTime, int endTime)
{
    int numGroups=Time[startTime][endTime].count;   //get the number of groups for this time interval
    //check if any of the groups has only 1 time series
    double msf;     //max so far
    vector<double>::iterator pos;
    //getting the position and max distance
    pos = max_element (Time[startTime][endTime].DcDistance.begin(), Time[startTime][endTime].DcDistance.end());
    msf=*pos;
    for(int i=0;i<numGroups;i++)
    {
        int groupID=Time[startTime][endTime].groupid[i];
        if(groupArray[groupID].count==1)
        {
            //check if this groups has maximum distance with others
            if(Time[startTime][endTime].DcDistance[i]==msf)
            {
                //return the TS in this group
                cout<<"Outlier is ";
                for(int j=0;j<N;j++)
                {
                    if(groupArray[groupID].bitVector[j]==1)
                    {
                        cout<<j<<endl;
                        return j;
                    }
                }
            }
        }
    }

    cout<<"No Outlier found"<<endl;
    return -1;
}

//reads time information from the file containing the count of groups in each time interval
void onlineOperation::readTime()
{

    inputTimeFile.open("data/Time.txt");
    //read till the end of file
    string temp;
    int startT;
    int endT;
    int count;
    double dist;
    while(inputTimeFile)
    {

        inputTimeFile>>temp;    //read T
        inputTimeFile>>startT;  //read start interval
        inputTimeFile>>endT;  //read start interval
        inputTimeFile>>count;  //count of groups in this interval
        Time[startT][endT].count=count;
        for(int i=0;i<count;i++)
        {
            inputTimeFile>>temp;    //read G
            int id;
            inputTimeFile>>id;
            Time[startT][endT].groupid.push_back(id);
        }
        inputTimeFile>>temp;    //read DcDist
        for(int i=0;i<count;i++)
        {
            inputTimeFile>>dist;
            Time[startT][endT].DcDistance.push_back(dist);
        }

    }
    inputFile.close();


}
//reads the groups attributes from a file and saves in group structure
void onlineOperation::readGroups()
{

    groupFile.open("data/Groups.txt");
    string temp;
    //int startT;
    //int endT;
    int count;
    onlineGroup *newGroup;//=new Group();
    groupFile>>temp;        //reading total
    groupFile>>count;
    groupFile>>temp;        //reading min
    groupFile>>min;
    groupFile>>temp;        //reading max
    groupFile>>max;
    while(groupFile.good())
    {
        newGroup=new onlineGroup();
        groupFile>>temp;        //reading G
        groupFile>>newGroup->id;    //group id
        groupFile>>newGroup->count; //count of TS
        groupFile>>newGroup->startT;
        groupFile>>newGroup->endT;
        //read bit vector
        for(int i=0;i<N;i++)
        {
            int bit;
            
            groupFile>>bit;
            newGroup->bitVector[i]=bit;
        }
        //read centroid
        groupFile>>count;   //centroid length
        newGroup->centroidLen=count;
        for(int i=0;i<count;i++)
        {
            double value;            
            groupFile>>value;
            newGroup->centroid.push_back(value);
        }
        groupArray.push_back(*newGroup);
        delete newGroup;
    }
    groupFile.close();
}

double distFunc(double x, double y)
{
        return sqrt(pow((x - y), 2));

}
//taken from http://bytefish.de/blog/dynamic_time_warping/
double onlineOperation::simpleDTW(const vector<double>& t1, const vector<double>& t2)
{
    //double temp;
    int m = t1.size();
            int n = t2.size();

            // create cost matrix
            double cost[m][n];

            cost[0][0] = distFunc(t1[0], t2[0]);
       //     temp=cost[0][0];
            // calculate first row
            for(int i = 1; i < m; i++)
            {

                cost[i][0] = cost[i-1][0] + distFunc(t1[i], t2[0]);
         //       temp=cost[i][0];
            }
            // calculate first column
            for(int j = 1; j < n; j++)
            {

                cost[0][j] = cost[0][j-1] + distFunc(t1[0], t2[j]);
           //     temp=cost[0][j];
            }
            // fill matrix
            for(int i = 1; i < m; i++)
            {
                for(int j = 1; j < n; j++)
                {
                    cost[i][j] = distFunc(t1[i],t2[j])+ std::min(cost[i-1][j],std::min(cost[i][j-1], cost[i-1][j-1]));
                }
            }

       //     temp=cost[m-1][n-1];
            return cost[m-1][n-1];
}
//implementation of spring outlier
int onlineOperation::springOutlier(int startT,int endT)
{
    //for each subsequence in this interval
    //get each subsequence
    double msf=0; //max distance so far
    int maxTS;    //Sequence
    double tempDist=0;
    vector<double> temp1;
    vector<double> temp2;
     for (int i=0;i<N;i++)
     {
         temp1.clear();
         temp2.clear();
         //copy  the subsequence in temp
         for(int p=startT;p<=endT;p++)
             temp1.push_back(timeSeries[i][p]);
         //compute the distance of this subsequence with all sequences
         for(int j=0;j<N;j++)
         {
             if(j!=i)
             {
                 for(int p=startT;p<=endT;p++)
                     temp2.push_back(timeSeries[j][p]);
                 tempDist+=simpleDTW(temp2,temp1);
                 temp2.clear();

             }
         }
         if(tempDist>msf)
         {
             msf=tempDist;
             maxTS=i;
         }
         tempDist=0;

     }
     cout<<"Outlier found by SPRING "<<maxTS<<" "<<msf<<endl;
     return maxTS;


}

//implementation of Spring
//query Type 1 mean k similar
void onlineOperation::springCalculation(int queryType, int k)
{
    double bsf=INF;
    int bsfIndex;       //index of bsf time series
    double currentDist=INF;
    int bestIntervalS;
    int bestIntervalE;
    vector<double> temp;    //temporary subsequence
    int kbestCount=0;
    kBest tempBest;

    if(queryType==1)
    {
        if(k==1)        //most similar
        {
            //get each subsequence
            for(int j=1;j<L;j++)
            {
                for (int l=0, m=l+j; m<L; l++, m++)
                {

                     for (int i=0;i<N;i++)
                     {
                         //copy  the subsequence in temp
                         for(int p=0;p<=m-l;p++)
                             temp.push_back(timeSeries[i][l+p]);
                         currentDist=simpleDTW(temp,tempQ);
                         if(currentDist<bsf)
                         {
                             bsf=currentDist;
                             bsfIndex=i;        //best time series index
                             bestIntervalS=l;   //record the best interval
                             bestIntervalE=m;
                         }
                         //clear the temp array
                         temp.clear();
                     }

                 }
             }
            cout<<"TS "<<bsfIndex<<" Interval "<<bestIntervalS<<" "<<bestIntervalE<<" Dist "<<bsf<<endl;
         }
        else
        {
            //have to find more than 1 similar time series
            //get each subsequence
            for(int j=1;j<L;j++)
            {
                for (int l=0, m=l+j; m<L-1; l++, m++)
                {

                     for (int i=0;i<N;i++)
                     {
                         //copy  the subsequence in temp
                         for(int p=0;p<m-l;p++)
                             temp.push_back(timeSeries[i][l+p]);
                         currentDist=simpleDTW(tempQ,temp);
                         if(kbestCount<k)
                         {
                             //add this TS to k best
                             kbestCount++;
                             tempBest.dist=currentDist;
                             tempBest.id=i;
                             kbestArray.push_back(tempBest);

                         }
                         else
                         {
                             sort(kbestArray.begin(),kbestArray.end(),sortByDist);
                             double tempD=kbestArray[kbestCount].dist;    //getting the last distance
                             if(tempD>currentDist)
                             {
                                 tempBest.dist=currentDist;
                                 tempBest.id=i;
                                 kbestArray[kbestCount]=tempBest;

                             }
                         }
                         temp.clear();
                     }

                 }
             }

        }
    }

}


/// Main Calculation Function
int onlineOperation::calculateTrillionDTW()
{
 //   FILE *fp;            /// data file pointer
 //   FILE *qp;            /// query file pointer
    double bsf;          /// best-so-far
    double *t, *q;       /// data array and query array
    int *order;          ///new order of the query
    double *u, *l, *qo, *uo, *lo,*tz,*cb, *cb1, *cb2,*u_d, *l_d;


    double d;
    long long i , j;
    double ex , ex2 , mean, std;
    int m=-1, r=-1;
    long long loc = 0;
    double t1,t2;
    int kim = 0,keogh = 0, keogh2 = 0;
    double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
    double *buffer, *u_buff, *l_buff;
    Index *Q_tmp;

    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    /// If not enough input, display an error.

    /// read size of the query
    m=queryLength;
    /// read warping windows
    r=floor(0.05);

    /// start the clock
    t1 = clock();

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

    for(int i=0;i<queryLength;i++)
    {
        d=tempQ[i];
        ex += d;
        ex2 += d*d;
        q[i] = d;

    }


    /// Do z-normalize the query, keep in same array, q
    mean = ex/m;
    std = ex2/m;
    std = sqrt(std-mean*mean);
    for( i = 0 ; i < m ; i++ )
         q[i] = (q[i] - mean)/std;

    /// Create envelop of the query: lower envelop, l, and upper envelop, u
    lower_upper_lemire(Q, m, r, l, u);

    /// Sort the query one time by abs(z-norm(q[i]))
    for( i = 0; i<m; i++)
    {
        Q_tmp[i].value = Q[i];
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
    int groupIndex=0;
    while(!done)
    {
        /// Read first m-1 points
        ep=0;
        int numPoints=groupArray[groupIndex].centroid.size();
        for(k=0; k<numPoints; k++)
                buffer[k] = d;

        /// Read buffer of size EPOCH or when all data has been read.
        ep=numPoints-1;
/*        while(ep<EPOCH)
        {   if (fscanf(fp,"%lf",&d) == EOF)
                break;
            buffer[ep] = d;
            ep++;
        }
*/
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
    //fclose(fp);

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
