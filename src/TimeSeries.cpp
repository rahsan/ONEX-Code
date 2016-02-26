#include "TimeSeries.h"

#include <algorithm>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


using namespace std;

SeriesDistanceMetric::SeriesDistanceMetric(const char *name,
                                           const char *desc,
                                           double (*calc)(TimeSeriesInterval &a,
                                                          TimeSeriesInterval &b, double dropout),
                                           double (*each)(double a, double b),
                                           double (*final)(double sum, int count))
{
    this->name = name;
    this->desc = desc;
    this->calc = calc;
    this->each = each;
    this->final = final;
}

double SeriesDistanceMetric::run(TimeSeriesInterval &a, TimeSeriesInterval &b, double dropout)
{
    if (calc != NULL) {
        return calc(a, b, dropout);
    }

    int len = a.length();

    if (b.length() != len) {
        printf("Attempted to get linear distance from two intervals of unequal size.\n");
        return INF;
    }

    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += each(a[i], b[i]);
        if (sum > dropout)
            return INF;
    }

    return final(sum, len);
}

double _abs(double a)
{
    return (a > 0)? a : -a;
}

double _naiveDTWcalc(TimeSeriesInterval &a, TimeSeriesInterval &b, double dropout=INF)
{
    int m = a.length();
    int n = b.length();

    // Fastpath for baseintervals.
    if (m == 1 && n == 1)
        return (b[0] - a[0]) * (b[0] - a[0]);

    // create cost matrix
    double cost[m][n];

    cost[0][0] = (b[0] - a[0]) * (b[0] - a[0]);

    // calculate first row
    for(int i = 1; i < m; i++) {

        cost[i][0] = cost[i-1][0] + (b[0] - a[i]) * (b[0] - a[i]);
    }

    // calculate first column
    for(int j = 1; j < n; j++) {

        cost[0][j] = cost[0][j-1] + (b[j] - a[0]) * (b[j] - a[0]);
    }

    // fill matrix. If using dropout, keep tabs on min cost per row.
    if (dropout != INF) {

        for(int i = 1; i < m; i++) {

            double min = cost[i][0];

            for(int j = 1; j < n; j++) {
                cost[i][j] = (b[j] - a[i]) * (b[j] - a[i])
                    + std::min(cost[i-1][j],
                               std::min(cost[i][j-1],
                                        cost[i-1][j-1]));

                min = std::min(min, cost[i][j]);
            }

            if (min > dropout) // Short circuit calculation.
                return INF;
        }

    } else {

        for(int i = 1; i < m; i++) {
            for(int j = 1; j < n; j++) {
                cost[i][j] = (b[j] - a[i]) * (b[j] - a[i])
                    + std::min(cost[i-1][j],
                               std::min(cost[i][j-1],
                                        cost[i-1][j-1]));
            }
        }
    }

    return cost[m-1][n-1];
}

SeriesDistanceMetric naiveDTWDist = SeriesDistanceMetric(
    "naivedtw",
    "DTW distance between vectors.",
    _naiveDTWcalc,
    NULL,
    NULL
);      

double _l2each(double a, double b)
{
    double diff = _abs(b - a);
    return diff * diff;
}

double _l2final(double sum, int count)
{
    return sqrt(sum/count);
}

SeriesDistanceMetric euclideanDist = SeriesDistanceMetric(
    "euclidean",
    "Euclidean (l2) distance between vectors.",
    NULL,
    _l2each,
    _l2final
);

SeriesDistanceMetric *availableDistMetrics[3] = {
    &naiveDTWDist,
    &euclideanDist,
    NULL
};

SeriesDistanceMetric *getDistMetric(const char *name)
{
    SeriesDistanceMetric **t;
    for (t = availableDistMetrics; *t != NULL; t++)
        if (!strcmp((*t)->name, name))
            return *t;

    return NULL;
}

void printDistMetrics(void)
{
    printf("Available distance metrics:\n");

    for (int i = 0; availableDistMetrics[i] != NULL; i++) {
        SeriesDistanceMetric *t = availableDistMetrics[i];
        printf("[%d] '%10s': %s\n", i, t->name, t->desc);
    }
}

TimeInterval::TimeInterval(int start, int end)
{
    this->start = start;
    this->end = end;
}

TimeInterval TimeInterval::subinterval(TimeInterval other)
{
    int s = start + other.start;
    int e = start + other.end;

    return TimeInterval(s, e);
}

int TimeInterval::length(void)
{
    return end - start + 1; // We're going inclusive, apparently. Ouch.
}


TimeSeriesInterval::TimeSeriesInterval(TimeSeriesSet *dataset, int seqNum, TimeInterval interval)
{
    this->dataset = dataset;

    this->interval = interval;

    this->seqNum = seqNum;

    data = dataset->getRawData(seqNum, interval.start);
}

TimeSeriesInterval::TimeSeriesInterval(seqitem_t *data , TimeInterval interval)
{
    this->dataset = NULL;
    this->interval = interval;
    this->data = data;
    this->seqNum = 0;

}

TimeSeriesInterval::~TimeSeriesInterval(void) {}

seqitem_t &TimeSeriesInterval::operator[](int index)
{
    return data[index];
}

seqitem_t *TimeSeriesInterval::getData(int index)
{
    return &data[index];
}

int TimeSeriesInterval::length(void)
{
    return interval.length();
}

TimeSeriesInterval TimeSeriesInterval::subinterval(TimeInterval interval)
{
    TimeInterval t = this->interval.subinterval(interval);
    return TimeSeriesInterval(data, t);
}


double TimeSeriesInterval::dist(TimeSeriesInterval &other, SeriesDistanceMetric *dist, double dropout)
{
    return dist->run(*this, other, dropout);
}

TimeSeriesSet::TimeSeriesSet(int seqCount, int seqLength, seqitem_t *data)
{
    this->seqCount = seqCount;
    this->seqLength = seqLength;

    int dataSize = seqCount * seqLength;

    if (data == NULL)
        this->data = (seqitem_t*) calloc(dataSize, sizeof(seqitem_t));
    else
        this->data = data;

    this->name = strdup("<memory>");
}

TimeSeriesSet::TimeSeriesSet(const char *path)
{
    ifstream in(path);

    if (!in.is_open()) {
        printf("Failed to open file for reading: %s.\n", path);
        seqCount = seqLength = 0;
        data = NULL;
        return;
    }

    in >> seqCount >> seqLength;
    if (seqCount <= 0 || seqLength <= 0) {
        printf("Failed to read file %s: Given invalid sequence count or length.\n", path);
        seqCount = seqLength = 0;
        data = NULL;
        return;
    }

    data = (seqitem_t*) calloc(seqLength * seqCount, sizeof(seqitem_t));

    seqitem_t minItem = INF, maxItem = -INF;
    seqitem_t t;
    for (int s = 0; s < seqCount; s++) {
        for (int i = 0; i < seqLength; i++) {
            in >> t;
            maxItem = std::max(maxItem, t);
            minItem = std::min(minItem, t);
            data[s*seqLength + i] = t;
        }
    }
    min = minItem;
    max = maxItem;

    name = strdup(path);
}

TimeSeriesSet::TimeSeriesSet(const char *path, int seqCount, int seqLength, int drop)
{
    ifstream in(path);

    if (!in.is_open()) {
        printf("Failed to open file for reading: %s.\n", path);
        seqCount = seqLength = 0;
        data = NULL;
        return;
    }

    this->seqLength = seqLength;
    this->seqCount = seqCount;

    data = (seqitem_t*) calloc(seqLength * seqCount, sizeof(seqitem_t));

    seqitem_t minItem = INF, maxItem = -INF;
    seqitem_t t;
    for (int s = 0; s < seqCount; s++) {
        for (int i = 0; i < drop; i++) in >> t;
        for (int i = 0; i < seqLength; i++) {
            in >> t;
            maxItem = std::max(maxItem, t);
            minItem = std::min(minItem, t);
            data[s*seqLength + i] = t;
        }
    }
    min = minItem;
    max = maxItem;

    name = strdup(path);
}

int TimeSeriesSet::toFile(const char *path, bool addSize)
{
    ofstream out(path);

    if (!out.is_open()) {
        printf("Failed to open file for writing: %s.\n", path);
        return -1;
    }

    if (addSize) {
        out << seqCount << " ";
        out << seqLength << " ";
        out << std::endl;
    }

    for (int s = 0; s < seqCount; s++) {
        for (int i = 0; i < seqLength; i++) {
            out << data[s*seqLength + i] << " ";
        }
        out << std::endl;
    }

    out.close();

    return 0;
}

TimeSeriesSet::~TimeSeriesSet(void)
{
    if (data)
        free(data);
    if (name)
        free(name);
}

TimeSeriesSet &TimeSeriesSet::randomSet(int seqCount, int seqLength, int range)
{
    int size = seqCount * seqLength;

    seqitem_t *data = (seqitem_t*) calloc(size, sizeof(seqitem_t));
    for (int i = 0; i < size; i++)
        data[i] = (rand() % range) + 1;

    TimeSeriesSet *t = new TimeSeriesSet(seqCount, seqLength, data);
    t->recalcMinMax();
    t->name = strdup("<random>");

    return *t;
}


int TimeSeriesSet::valid(void)
{
    return data != NULL;
}

TimeSeriesInterval TimeSeriesSet::getInterval(int seqNum, TimeInterval interval)
{
    return TimeSeriesInterval(this, seqNum, interval);
}

void TimeSeriesSet::zero(void)
{
    for (int i = 0; i < seqCount * seqLength; i++)
        data[i] = 0.0;
    min = max = 0;
}

void TimeSeriesSet::normalize(void)
{
    double diff = max - min;

    if (diff == 0.0)
        return;

    for (int i = 0; i < seqCount * seqLength; i++)
        data[i] = (data[i] - min)/diff;

    min = 0.0;
    max = 1.0;
}

void TimeSeriesSet::recalcMinMax(void)
{
    seqitem_t minItem = INF;
    seqitem_t maxItem = -INF;

    for (int i = 0; i < seqCount * seqLength; i++) {
        minItem = std::min(minItem, data[i]);
        maxItem = std::max(maxItem, data[i]);
    }

    this->min = minItem;
    this->max = maxItem;
}

seqitem_t *TimeSeriesSet::getRawData(int seqNum, int seqIndex)
{
    return &data[(seqNum * seqLength) + seqIndex];
}

seqitem_t &TimeSeriesSet::getData(int seqNum, int seqIndex)
{
    return data[(seqNum * seqLength) + seqIndex];
}


void TimeSeriesSet::printDesc(ostream *out)
{
    *out << "Time series set:" << endl;
    *out << "Name: '" << name << "'" << endl;
    *out << "Sequences: " << seqCount << endl;
    *out << "Sequence length: " << seqLength << endl;
    *out << "Min / max values: " << min << ", " << max << endl;
}

void TimeSeriesSet::printData(ostream *out)
{
    int w, p;
    w = out->width();
    p = out->precision();
    *out << "Printing time series set '" << name;
    *out << "' with " << seqCount <<" sequences of length " << seqLength << ":";

    for (int s = 0; s < seqCount; s++) {
        *out << endl << "Sequence[" << s << "]:";
        for (int i = 0; i < seqLength; i++) {
            *out << " ";
            out->width(6);
            out->precision(4);
            *out << data[(s * seqLength) + i];
            out->width(w);
            out->precision(p);
        }
    }
    *out << endl << "End of series set." << endl;
}
