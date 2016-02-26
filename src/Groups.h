#ifndef GROUPS_H
#define GROUPS_H

#include <iostream>
#include <stdio.h>
#include<fstream>
#include <vector>
#include <math.h>

using namespace std;

//const int timeSeriesN=24;      //ECG
//const int timeSeriesL=96;
const int timeSeriesN=67;       //Italy Power
const int timeSeriesL=24;
//const int timeSeriesN=560;      //Synthetic
//const int timeSeriesL=132;

//const int timeSeriesN=67;      //Small dataset
//const int timeSeriesL=24;
//distance of centroid with the centroids is the same group
struct intraDistance
{
    double distance;
    double flag;      //flag true means +ve, false means -ve
};

struct Group
{
    int id;     //group number
    int startT; //start time interval
    int endT;   //end time interval
    int count;  //count of time series in this group
    int bitVector[timeSeriesN];   //bit vector containing 1 or 0 for each time series
    vector<double> centroid;    //centroid element array
    vector<double> sum;         //maintains the sum of each point of TS in this group to calculate centroid
    int centroidLen;    //number of points in centroid
//    intraDistance seqCentroidDist[timeSeriesN];   //distance of sequence to centroid

};

struct timeEntry
{
    vector<int> groupid;       //an array containing the id of each group in this time
    int count;          //count of groups
    vector<double> DcDistance;  //sum of ED distance with other centroids
};

class GroupOperation
{
public:
    // Construction, destruction, and file initialization.
    GroupOperation(double STinput);
    ~GroupOperation();

    int readFile(char * Filename, int N, int length);    //filename, number of time series and length of each time series

    // Print operations.
    void printTS();     //prints the time series
    void printTime(int level); //prints all the times and asociated information (group ids, count of groups)
    void printGroups(); //prints alll groups with associated information
    void printOneGroup(int groupId);    //prints information of the group id specified

    //group Operations
    void groupOp();     //groups the time series
    void intersectGroup(int level);  //finds the intersection of groups
    void groupByTime(int startTime,int endTime);    //groups all the time series for the specified time interval
    double distED(int TSStartIndex, int TSEndIndex, double centroidStartValue, double centroidEndValue, int m); //calculate ED between TS and centroid
    void normalize();               //normalize the time series

    //variables
    int groupCounter;       //maintains the the total number of groups so far.
    int N;  //number of time series
    int L;  //length of time series
    double ST;            //Similarity threshold
    double min;         //min value in data
    double max;         //max value in data for normalization
    int centroidDist[timeSeriesN][timeSeriesN];     //Dc distance
    timeEntry Time[timeSeriesL][timeSeriesL];       //Time
    double timeSeries[timeSeriesN][timeSeriesL];
    ifstream inputFile;         //input data file
    ifstream queryFile;         //input query file
    ofstream timeFile;          //output the Global Time index
    ofstream groupFile;         //output file for groups and its attributes
    //ofstream minmaxFile;        //output min and max for dataset
 //   Group *myGroup;
    vector<Group> groupArray;   //

};

#endif // GROUPS_H
