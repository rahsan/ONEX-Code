#ifndef GROUPING_H
#define GROUPING_H

#include "TimeSeries.h"

#include <vector>

#include <stdint.h>

using namespace std;

struct kBest
{
    int id;         //id of best time series
    double dist;    //distance
};

// All slices from a dataset in a given TimeInterval.
class TimeSeriesSetSlice
{
protected:
    TimeSeriesSet *set;
    TimeInterval slice;

public:

    TimeSeriesSetSlice(TimeSeriesSet *dataset, TimeInterval slice);

    TimeSeriesInterval operator[](int index);

    int length(void);
    int seqCount(void);


    // Spring calculation methods.
    int springOutlier(void); // Returns farthest series number.
    vector<kBest> springkSimilar(vector<double> tempQ, int queryType, int k);
};


// A centroid for a time series interval. Efficient centroid calculation.
class TimeSeriesIntervalCentroid
{
protected:
    int count, length;
    vector<seqitem_t> sum;
    bool cacheValid;
    vector<seqitem_t> cachedCentroid;
public:
    TimeSeriesIntervalCentroid(int length);
    TimeSeriesIntervalCentroid(const TimeSeriesIntervalCentroid &other);
    TimeSeriesIntervalCentroid(vector<seqitem_t> &centroid, int n);

    int getCount(void);
    int getSliceLength(void);

    void addVector(const vector<seqitem_t> data);
    void addArray(const seqitem_t *data);

    vector<seqitem_t> &getCentroid(void);
};

// Implements trillionDTW style cascading lower bounds comparison. See README for references.
class TimeSeriesIntervalEnvelope
{
protected:
    //bool order_cache_valid;
    //vector<int> index_order;
    //void genIndexOrder(void);

    bool kimfl_cache_valid;
    seqitem_t kimfl_f, kimfl_l, kimfl_min, kimfl_max;
    void genKimFL(void);

    bool keogh_cache_valid;
    vector<seqitem_t> keogh_lower, keogh_upper;
    void genKeoghLU(void);
public:
    TimeSeriesInterval interval;

 TimeSeriesIntervalEnvelope(void):interval(NULL,TimeInterval()){}
    TimeSeriesIntervalEnvelope(TimeSeriesInterval interval);

    void setInterval(TimeSeriesInterval interval);

    //vector<int> &getIndexOrder(void);
    vector<seqitem_t> &getKeoghLower(void);
    vector<seqitem_t> &getKeoghUpper(void);

    seqitem_t kimFLDist(TimeSeriesIntervalEnvelope &other, double dropout=INF);
    seqitem_t keoghDist(TimeSeriesIntervalEnvelope &other, double dropout=INF);
    seqitem_t crossKeoghDist(TimeSeriesIntervalEnvelope &other, double dropout=INF);
    seqitem_t cascadeDist(TimeSeriesIntervalEnvelope &other, double dropout=INF);
};

struct kSim {
    double distance;
    int index;

    bool operator<(kSim &other)
    {
        return distance < other.distance;
    }
};


// vector<bool> optimized for our usage.
struct BoolVector
{
    int size;
    int32_t *data;

    BoolVector(int size);
    BoolVector(const BoolVector &other);
    ~BoolVector(void);

    bool get(int index);
    void set(int index, bool value=1);
    void reset(int index, bool value=0);

    vector<int> intersection(const BoolVector &other);
};


// A subgroup of a TimeSeriesSet slice. Used in grouping. Calculates centroid.
class TimeSeriesSliceGroup
{
protected:
    TimeSeriesIntervalCentroid centroid;
    BoolVector member;

    bool interval_cache_valid;
    TimeSeriesIntervalEnvelope cached_interval;

public:
    TimeSeriesSetSlice *slice;
    int count;

    TimeSeriesSliceGroup(TimeSeriesSetSlice *setslice);
    //TimeSeriesSliceGroup(const TimeSeriesSliceGroup &other);
    ~TimeSeriesSliceGroup(void);

    void addMember(int index, bool update=true);
    bool isMember(int index);

    //vector<int> members(void);
    vector<int> intersect(TimeSeriesSliceGroup *other);

    // Distance of interval from this centroid
    seqitem_t distance(int index, SeriesDistanceMetric *metric=&euclideanDist, seqitem_t dropout=INF);
    seqitem_t distance(vector<seqitem_t> data, SeriesDistanceMetric *metric=&euclideanDist, seqitem_t dropout=INF);

    TimeSeriesIntervalEnvelope &getEnvelope(void);

    vector<kSim> getSortedSimilar(TimeSeriesIntervalEnvelope query, int count, double dropout=INF);

    vector<seqitem_t> getCentroid(void);

    void toFile(ostream *out);
    void fromFile(istream *in);
};


// A complete (disjoint) grouping for a TimeSeriesSet slice.
class TimeSeriesSliceGrouping
{
protected:
    TimeSeriesSetSlice slice;

    void updateDiffs(void);

public:
    vector<TimeSeriesSliceGroup*> groups;
    vector<double> centroidTotalDiffs;
    double maxDiff;

    TimeSeriesSliceGrouping(TimeSeriesSetSlice setslice);

    void hgroup(TimeSeriesSliceGrouping *left, TimeSeriesSliceGrouping *right);
    void naiveGroup(seqitem_t ST);

    // Get the index of the group with the closest centroid. Returns -1 if none or >dropout.
    int getBestGroup(TimeSeriesIntervalEnvelope query, double *dist, double dropout=INF);

    void free(void);

    void toFile(ostream *out);
    void fromFile(istream *in); // Must already have ST and slice ready.
};


// A grouping for every TimeSeriesSet in a time series. (All intervals are grouped.)
class TimeSeriesGrouping
{
protected:
    seqitem_t ST;
    TimeSeriesSet &set;

public:
    TimeSeriesSliceGrouping *groups; // Triangular, [start*seqLength + end].

    TimeSeriesGrouping(TimeSeriesSet &database, seqitem_t ST);
    ~TimeSeriesGrouping(void);

    bool valid(void);
    void naiveGroup(bool debug);
    void hgroup(bool debug);

    int fromFile(const char *path);
    int toFile(const char *path);

    TimeSeriesSliceGrouping &getGroup(TimeInterval slice);

    seqitem_t getST(void);
    void setST(seqitem_t ST);
};

#endif // GROUPING_H
