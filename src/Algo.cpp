#include "Algo.h"
#include "TimeSeries.h"
#include "Grouping.h"

#include <iostream>
#include <algorithm>

double _distFunc(double x, double y)
{
    double t = x - y;
    return (t < 0)? -t : t;
}

//taken from http://bytefish.de/blog/dynamic_time_warping/
double simpleDTW(const vector<double>& t1, const vector<double>& t2)
{
    //double temp;
    int m = t1.size();
    int n = t2.size();

    // create cost matrix
    double cost[m][n];

    cost[0][0] = _distFunc(t1[0], t2[0]);
    //     temp=cost[0][0];
    // calculate first row
    for(int i = 1; i < m; i++)
        {

            cost[i][0] = cost[i-1][0] + _distFunc(t1[i], t2[0]);
            //       temp=cost[i][0];
        }
    // calculate first column
    for(int j = 1; j < n; j++)
        {

            cost[0][j] = cost[0][j-1] + _distFunc(t1[0], t2[j]);
            //     temp=cost[0][j];
        }
    // fill matrix
    for(int i = 1; i < m; i++)
        {
            for(int j = 1; j < n; j++)
                {
                    cost[i][j] = _distFunc(t1[i],t2[j])+ std::min(cost[i-1][j],std::min(cost[i][j-1], cost[i-1][j-1]));
                }
        }

    //     temp=cost[m-1][n-1];
    return cost[m-1][n-1];
}

//implementation of spring outlier
int TimeSeriesSetSlice::springOutlier(void)
{
    double *timeSeries = set->getRawData(0, 0);
    int seqLen = set->seqLength;
    int seqCnt = set->seqCount;

    int startT = slice.start;
    int endT = slice.end;
    //for each subsequence in this interval
    //get each subsequence
    double msf = 0; //max distance so far
    int maxTS = -1;    //Sequence
    double tempDist = 0;
    vector<double> temp1;
    vector<double> temp2;
     for (int i = 0; i < seqCnt; i++)
     {
         temp1.clear();
         temp2.clear();
         //copy  the subsequence in temp
         for (int p = startT; p <= endT; p++)
             temp1.push_back(timeSeries[i * seqLen + p]);
         //compute the distance of this subsequence with all sequences
         for (int j = 0; j < seqCnt; j++)
         {
             if (j!=i)
             {
                 for (int p = startT; p <= endT; p++)
                     temp2.push_back(timeSeries[j * seqLen + p]);
                 tempDist+=simpleDTW(temp2,temp1);
                 temp2.clear();
             }
         }
         if (tempDist > msf)
         {
             msf=tempDist;
             maxTS=i;
         }
         tempDist=0;
     }

     cout << "Outlier found by SPRING " << maxTS << " " << msf << endl;

     return maxTS;
}

bool _sortByDist(const kBest &lhs, const kBest &rhs)
{
    return lhs.dist < rhs.dist;
}

//implementation of Spring
//query Type 1 mean k similar
vector<kBest> TimeSeriesSetSlice::springkSimilar(vector<double> tempQ, int queryType, int k)
{
    vector<kBest> kbestArray;
    double *timeSeries = set->getRawData(0, 0);
    int N = set->seqCount;
    int L = set->seqLength;

    double bsf=INF;
    int bsfIndex = -1;       //index of bsf time series
    double currentDist=INF;
    int bestIntervalS = -1;
    int bestIntervalE = -1;
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
                             temp.push_back(timeSeries[i*L + l+p]);
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
                             temp.push_back(timeSeries[i*L + l+p]);
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
                             sort(kbestArray.begin(),kbestArray.end(), _sortByDist);
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
    return kbestArray;
}
