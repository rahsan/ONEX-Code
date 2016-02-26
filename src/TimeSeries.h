#ifndef TIMESERIES_H
#define TIMESERIES_H


#include <iostream>
#include <stdlib.h>


#define INF 1e20       //Pseudo Infitite number for this code
typedef double seqitem_t;


struct SeriesDistanceMetric;
struct TimeInterval;
class TimeSeriesInterval;
class TimeSeriesSet;



//// Distance metrics.
struct SeriesDistanceMetric
{
public:
    // For pretty printing.
    const char *name;
    const char *desc;

    // Find the distance betwee two intervals.
    double run(TimeSeriesInterval &a, TimeSeriesInterval &b, double dropout=INF);

    // Because initialization is bothersome.
    SeriesDistanceMetric(const char *name,
                         const char *desc,
                         double (*calc)(TimeSeriesInterval &a, TimeSeriesInterval &b, double dropout),
                         double (*each)(double a, double b),
                         double (*final)(double sum, int count));

protected:
    // If not NULL, use instead of linear method.
    double (*calc)(TimeSeriesInterval &a, TimeSeriesInterval &b, double dropout);

    double (*each)(double a, double b);
    double (*final)(double sum, int count);
};

// Available distance metrics.
extern SeriesDistanceMetric naiveDTWDist; // DTW
extern SeriesDistanceMetric taxicabDist; // l1 distance metric.
extern SeriesDistanceMetric euclideanDist; // l2 distance metric.

extern SeriesDistanceMetric *availableDistMetrics[]; // List of available distance methods. NULL terminated.

SeriesDistanceMetric *getDistMetric(const char *name);
void printDistMetrics(void);



//// Time interval [a,b].
struct TimeInterval
{
    int start;
    int end;

    TimeInterval(int start=0, int end=0);

    TimeInterval subinterval(TimeInterval other);
    int length(void);
};



//// Reference to a slice of a timeseries.
class TimeSeriesInterval
{
protected:
    seqitem_t *data;

    int seqNum;

    TimeInterval interval;

public:

    TimeSeriesSet *dataset;

    TimeSeriesInterval(TimeSeriesSet *dataset, int seqNum, TimeInterval interval);
    TimeSeriesInterval(seqitem_t *data, TimeInterval interval);
    ~TimeSeriesInterval(void);
    int length(void);

    seqitem_t &operator[](int index); // Get a reference to a part of the index.
    seqitem_t *getData(int index = 0); // Get an actual pointer to the data.
    TimeSeriesInterval subinterval(TimeInterval interval); // Get a slice of this interval.

    double dist(TimeSeriesInterval &other, SeriesDistanceMetric *dist=&euclideanDist, double dropout=INF);
};



//// Time series database.
class TimeSeriesSet
{
public:

    seqitem_t *data;
    seqitem_t min, max;

    int seqCount;
    int seqLength;

    TimeSeriesSet(int seqCount, int seqLength, seqitem_t *data=NULL);
    TimeSeriesSet(const char *path);
    TimeSeriesSet(const char *path, int seqCount, int seqLength, int drop=0);
    ~TimeSeriesSet(void);

    int toFile(const char *path, bool addSize=true);


    static TimeSeriesSet &randomSet(int seqCount=100, int seqLength=800, int range=100);


    void normalize(void);
    void zero(void);
    void recalcMinMax(void);

    int valid(void); // Did we load properly?
    TimeSeriesInterval getInterval(int seqNum, TimeInterval interval);
    seqitem_t &getData(int seqNum, int seqIndex);
    seqitem_t *getRawData(int seqNum, int seqIndex);

    void printDesc(std::ostream *out); // prints seqcount/length, min/max.
    void printData(std::ostream *out); // prints *everything*

    char *name; // May not be available.
};


#endif // TIMESERIES_H
