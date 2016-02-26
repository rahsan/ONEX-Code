#ifndef ONLINE_H
#define ONLINE_H

#include <vector>
#include <set>
#include <string>
#include <map>
#include <fstream>

using namespace std;
//const int SeriesN=3;
//const int SeriesL=5;
const int SeriesN=67;   //Italy Power
const int SeriesL=24;
//#define INF 1e20       //Pseudo Infitinte number for this code
//const int SeriesN=1;   //ECG
//const int SeriesL=10;
//const int SeriesN=400;   //Synthetic
//const int SeriesL=100;


//const int SeriesN=3;
//const int SeriesL=24;

struct onlineintraDistance
{
    double distance;   //distance with centroid
    bool flag;      //flag true means +ve, false means -ve
    bool operator< (const onlineintraDistance  &e) const {return distance > e.distance;}
};
struct onlinetimeEntry
{
    vector<int> groupid;       //an array containing the id of each group in this time
    int count;          //count of groups
    vector<double> DcDistance;  //sum of ED distance with other centroids
};
struct kBest
{
    int id;         //id of best time series
    double dist;    //distance
};

struct onlineGroup
{
    int id;     //group number
    int startT; //start time interval
    int endT;   //end time interval
    int count;  //count of time series in this group
    int bitVector[SeriesN];   //bit vector containing 1 or 0 for each time series
    vector<double> centroid;    //centroid element array
    vector<double> sum;         //maintains the sum of each point of TS in this group to calculate centroid
    int centroidLen;    //number of points in centroid
 //   set<onlineintraDistance> seqCentroidDist;   //distance of sequences to centroid

};

class onlineOperation
{
public:
    onlineOperation(double STnew, int Ndata, int Ldata);
    ~onlineOperation();

    void readTime();
    void readGroups();
    void readQueryFile(int m, const char * Filename, bool type);
    int readTimeSeries(const char* Filename, bool type);

    void kSimilar(bool queryFlag, int k,bool timeFlag,int startTime, int endTime);
    int dominantOutlier(int startTime, int endTime);
    void criticalTime(bool queryFlag, double epsilon);
    void stableTimeInterval(bool queryFlag, double epsilon);
    void evolvingSimilarity(bool queryFlag);
    void simThresholdRec(bool simDegreeFlag, char category);


    int calculate();
    double simpleDTW(const vector<double>& t1, const vector<double>& t2);
    void springCalculation(int queryType, int k);
    int calculateTrillionDTW();
    int springOutlier(int startT,int endT);


public:
    //variables
    ifstream inputTimeFile;         //input Time data file
    ifstream queryFile;         //input query file
    ifstream groupFile;         //input Group file
    ifstream inputFile;         //input time series data file

    vector<onlineGroup> groupArray;
    vector<kBest> kbestArray;

    int N;  //number of time series
    int L;  //length of time series
    double ST;            //Similarity threshold
    double min;         //min value in data
    double max;         //max value in data for normalization
    int centroidDist[SeriesN][SeriesN];     //Dc distance
    onlinetimeEntry Time[SeriesL][SeriesL];       //Time
    double timeSeries[SeriesN][SeriesL];

    double *Q;      //query data
    int queryLength;          //size of query data
    vector<double> tempQ;
};

#endif // ONLINE_H
