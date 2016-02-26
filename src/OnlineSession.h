#ifndef ONLINESESSION_H
#define ONLINESESSION_H

#include "TimeSeries.h"
#include "Grouping.h"

#include <map>
#include <string>
#include <iostream>
#include <vector>

using namespace std;

class OnlineSession
{
protected:
    vector<TimeSeriesSet*> dataSets;
    vector<TimeSeriesGrouping*> groupings;

public:
    OnlineSession(double defaultST=1.5);
    ~OnlineSession(void);

    int loadDataSet(TimeSeriesSet *data);
    int saveDataSet(int index, const char *path, bool old=false);
    int dropDataSet(int index); // nonzero on error.

    int normalizeDataSet(int index);

    void genGrouping(int index, seqitem_t ST);
    void genNaiveGrouping(int index, seqitem_t ST);
    int loadGrouping(int index, const char *path);
    int saveGrouping(int index, const char *path);
    void dropGrouping(int index);

    double defaultST;
    int debug;

    void printDataSets(ostream *out);
    void checkIndex(int index);

    int run(istream *in, ostream *out, bool interactive=true);

    void kSimilar(int dbindex, vector<double> qdata, TimeInterval interval, int k, int strict);
    void dominantOutlier(int dbindex, TimeInterval interval);
    void criticalTime(int dbindex, double epsilon);
    void evolvingSimilarity(int dbindex);
    void simThresholdRec(bool simDegreeFlag, char category);
};

#endif // ONLINESESSION_H
