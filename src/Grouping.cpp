#include "Grouping.h"

#include "trillionDTW.h"

#include <iostream>
#include <fstream>
#include <algorithm>

#include <stdlib.h>

TimeSeriesSetSlice::TimeSeriesSetSlice(TimeSeriesSet *dataset, TimeInterval slice)
{
    this->set = dataset;
    this->slice = slice;
}

TimeSeriesInterval TimeSeriesSetSlice::operator[](int index)
{
    return set->getInterval(index, slice);
}

int TimeSeriesSetSlice::length(void)
{
    return slice.length();
}

int TimeSeriesSetSlice::seqCount(void)
{
    return set->seqCount;
}

TimeSeriesIntervalCentroid::TimeSeriesIntervalCentroid(int length)
{
    count = 0;
    this->length = length;
    sum = vector<seqitem_t>(length, 0.0);
    cachedCentroid = vector<seqitem_t>(length, 0.0);
    cacheValid = true;
}

TimeSeriesIntervalCentroid::TimeSeriesIntervalCentroid(const TimeSeriesIntervalCentroid &other):sum(other.sum),
                                                                                                cachedCentroid(other.cachedCentroid)
{
    count = other.count;
    length = other.length;
    cacheValid = other.cacheValid;
}

TimeSeriesIntervalCentroid::TimeSeriesIntervalCentroid(vector<seqitem_t> &centroid, int n)
{
    count = n;
    length = centroid.size();
    cachedCentroid = vector<seqitem_t>(centroid);
    sum = vector<seqitem_t>(length, 0);

    for (int i = 0; i < length; i++)
        sum[i] = centroid[i] * count;

    cacheValid = true;
}

int TimeSeriesIntervalCentroid::getCount(void)
{
    return count;
}

int TimeSeriesIntervalCentroid::getSliceLength(void)
{
    return length;
}

void TimeSeriesIntervalCentroid::addVector(const vector<seqitem_t> data)
{
    addArray(data.data());
}

void TimeSeriesIntervalCentroid::addArray(const seqitem_t *data)
{
    for (unsigned int i = 0; i < sum.size(); i++)
        sum[i] += data[i];

    count++;
    cacheValid = false;
}

vector<seqitem_t> &TimeSeriesIntervalCentroid::getCentroid(void)
{
    if (cacheValid == true) {
        return cachedCentroid;
    }

    for (unsigned int i = 0; i < sum.size(); i++) {
        cachedCentroid[i] = sum[i] / count;
    }

    cacheValid = true;
    return cachedCentroid;
}

TimeSeriesIntervalEnvelope::TimeSeriesIntervalEnvelope(TimeSeriesInterval interval): interval(interval)
{
    /*order_cache_valid = */
    kimfl_cache_valid = keogh_cache_valid = false;
}

void TimeSeriesIntervalEnvelope::setInterval(TimeSeriesInterval interval)
{
    this->interval = interval;
    /*order_cache_valid = */
    kimfl_cache_valid = keogh_cache_valid = false;
}

/*
vector<int> &TimeSeriesIntervalEnvelope::getIndexOrder(void)
{
    if (!order_cache_valid) {
        genIndexOrder();
    }

    return index_order;
}
*/

vector<seqitem_t> &TimeSeriesIntervalEnvelope::getKeoghLower(void)
{
    if (!keogh_cache_valid) {
        genKeoghLU();
    }

    return keogh_lower;
}

vector<seqitem_t> &TimeSeriesIntervalEnvelope::getKeoghUpper(void)
{
    if (!keogh_cache_valid) {
        genKeoghLU();
    }

    return keogh_upper;
}

//void TimeSeriesIntervalEnvelope::genIndexOrder(void) {}
// TODO: Remove, we don't use this.
void TimeSeriesIntervalEnvelope::genKimFL(void)
{
    kimfl_f = interval[0];
    kimfl_l = interval[interval.length() - 1];

    double min = INF;
    double max = -INF;
    for (int i = 0; i < interval.length(); i++) {
        double t = interval[i];
        max = (t > max)? t : max;
        min = (t < min)? t : min;
    }

    kimfl_min = min;
    kimfl_max = max;

    kimfl_cache_valid = true;
}

void TimeSeriesIntervalEnvelope::genKeoghLU(void)
{
    keogh_lower.resize(interval.length(), 0);
    keogh_upper.resize(interval.length(), 0);

    int R = 0; // What value to use?

    // Function provided by trillionDTW codebase. See README.
    lower_upper_lemire(interval.getData(), interval.length(), R,
                       keogh_lower.data(), keogh_upper.data());

    keogh_cache_valid = true;
}

inline double _min(double a, double b)
{
    return (a < b)? a : b;
}

inline double _dist(double a, double b)
{
    return (b - a) * (b - a);
}

// Methods based on trillionDTW codebase. See README.
seqitem_t TimeSeriesIntervalEnvelope::kimFLDist(TimeSeriesIntervalEnvelope &other, double dropout)
{
    int al = interval.length();
    int bl = other.interval.length();
    int l = _min(al,bl);

    if (l == 0)
        return 0;

    if (l == 1)
        return _dist(interval[0], other.interval[0]);

    double lb = 0;

    lb += _dist(interval[0], other.interval[0]);
    lb += _dist(interval[al-1], other.interval[bl-1]);
    if (lb >= dropout)
        return INF;

    lb += _min(_min(_dist(interval[0], other.interval[1]),
                    _dist(interval[1], other.interval[1])),
               _dist(interval[1], other.interval[0]));

    if (lb >= dropout)
        return INF;

    lb += _min(_min(_dist(interval[al-1], other.interval[bl-2]),
                    _dist(interval[al-2], other.interval[bl-2])),
               _dist(interval[al-2], other.interval[bl-1]));

    if (lb >= dropout)
        return INF;

    if (l == 4)
        return lb;

    lb += _min(_min(_min(_dist(interval[0], other.interval[2]),
                         _dist(interval[1], other.interval[2])),
                    _min(_dist(interval[2], other.interval[2]),
                         _dist(interval[2], other.interval[1]))),
               _dist(interval[2], other.interval[0]));

    if (lb >= dropout)
        return INF;

    lb += _min(_min(_min(_dist(interval[al-1], other.interval[bl-3]),
                         _dist(interval[al-2], other.interval[bl-3])),
                    _min(_dist(interval[al-3], other.interval[bl-3]),
                         _dist(interval[al-3], other.interval[bl-2]))),
               _dist(interval[al-3], other.interval[bl-1]));

    return lb;
}

// Methods based on trillionDTW implementation. See README.
seqitem_t TimeSeriesIntervalEnvelope::keoghDist(TimeSeriesIntervalEnvelope &other, double dropout)
{
    if (!other.keogh_cache_valid)
        other.genKeoghLU();

    int al = interval.length();
    int bl = other.interval.length();
    int len = _min(al, bl);

    double lb = 0;
    double x, d, l, u;

    for (int i = 0; i < len && lb < dropout; i++)
    {
        x = interval[i];
        l = other.keogh_lower[i];
        u = other.keogh_upper[i];
        d = 0;
        if (x > u)
            d = _dist(x,u);
        else if(x < l)
            d = _dist(x,l);
        lb += d;
    }

    return lb;
}

// Based on trillionDTW methods. See README.
seqitem_t TimeSeriesIntervalEnvelope::crossKeoghDist(TimeSeriesIntervalEnvelope &other, double dropout)
{
    double lb = keoghDist(other, dropout);

    if (lb >= dropout)
        return INF;
    else
        return lb + other.keoghDist(*this, dropout - lb);
}

// Based on trillionDTW methods. See README.
seqitem_t TimeSeriesIntervalEnvelope::cascadeDist(TimeSeriesIntervalEnvelope &other, double dropout)
{
    double lb = kimFLDist(other, dropout);
    if (lb >= dropout)
        return INF;

    lb += crossKeoghDist(other, dropout - lb);
    if (lb >= dropout)
        return INF;

    return interval.dist(other.interval, &naiveDTWDist, dropout - lb);
}



BoolVector::BoolVector(int size)
{
    this->size = size;
    data = (int32_t*) calloc(size/32 + 1, sizeof(int32_t));
}

BoolVector::BoolVector(const BoolVector &other)
{
    size = other.size;
    data = (int32_t*) calloc(size/32 + 1, sizeof(int32_t));
    for (int i = 0; i < size/32 + 1; i++)
        data[i] = other.data[i];
}
BoolVector::~BoolVector(void)
{
    size = 0;
    if (data)
        free(data);
}

bool BoolVector::get(int index)
{
    int i = index / 32;
    int32_t d = data[i];

    return (d & (1 << index & 0x1F))? 1 : 0;
}

void BoolVector::set(int index, bool value)
{
    int i = index / 32;
    if (value)
        data[i] |= 1 << (index & 0x1F);
    else
        data[i] &= ~(1 << (index & 0x1F));
}

void BoolVector::reset(int index, bool value)
{
    int i = index / 32;
    if (!value)
        data[i] |= 1 << (index & 0x1F);
    else
        data[i] &= ~(1 << (index & 0x1F));
}

vector<int> BoolVector::intersection(const BoolVector &other)
{
    int32_t *da = data;
    int32_t *db = other.data;

    vector<int> ret;

    for (int i = 0; i < (size / 32) + 1; i++) {
        int32_t d = da[i] & db[i];
        if (d != 0) {
            int32_t o = 0x0001;
            for (int j = 0; j < 32; j++) {
                if (o & d) {
                    ret.push_back(i*32 + j);
                }
                o = o << 1;
            }
        }
    }

    return ret;
}

TimeSeriesSliceGroup::TimeSeriesSliceGroup(TimeSeriesSetSlice *slice):centroid(slice->length()),
                                                                      member(slice->seqCount())
{
    interval_cache_valid = false;
    this->slice = slice;
    count = 0;
}

/*
TimeSeriesSliceGroup::TimeSeriesSliceGroup(const TimeSeriesSliceGroup &other):centroid(other.centroid),
                                                                              member(slice->seqLength())
{
    centroid = other.centroid;
    this->slice = other.slice;
    count = other.count;
}
*/

TimeSeriesSliceGroup::~TimeSeriesSliceGroup(void)
{
}

bool TimeSeriesSliceGroup::isMember(int index)
{
    return member.get(index);
}

void TimeSeriesSliceGroup::addMember(int index, bool update)
{
    member.set(index, 1);
    count++;
    if (update)
        centroid.addArray((*slice)[index].getData());
    interval_cache_valid = false;
}

/*vector<int> TimeSeriesSliceGroup::members(void)
{
    vector<int> t;
    for (int i = 0; i < (slice->seqLength())/8+1; i++) {
        int8_t m = ((int8_t*)member)[i];
        int8_t o = 0x1;
        for(int j = 0; j < 8; j++, o = o << 1) {
            if (m & o) {
                t.push_back(i * 8 + j);
            }
        }
    }

    return t;
    }*/

// Interesting bit array optimizations. Ignore.
vector<int> TimeSeriesSliceGroup::intersect(TimeSeriesSliceGroup *other)
{
    return member.intersection(other->member);
}

vector<kSim> TimeSeriesSliceGroup::getSortedSimilar(TimeSeriesIntervalEnvelope query, int count, double dropout)
{
    vector<kSim> distances(slice->seqCount());
    kSim bsf;
    bsf.distance = dropout;
    bsf.index = -1;

    for (int i = 0; i < slice->seqCount(); i++) {
        if (!member.get(i))
            distances[i].distance = INF;
        else
            distances[i].distance = TimeSeriesIntervalEnvelope((*slice)[i])
                .cascadeDist(query, dropout);
        distances[i].index = i;
        if (count == 1) {
            dropout = _min(dropout, distances[i].distance);
            if (distances[i].distance < bsf.distance) {
                bsf.distance = distances[i].distance;
                bsf.index = i;
            }
        }
    }

    if (count == 1)
        return vector<kSim>(1, bsf);

    sort(distances.begin(), distances.end());

    distances.resize(count);

    return distances;
}

TimeSeriesIntervalEnvelope &TimeSeriesSliceGroup::getEnvelope(void)
{
    if (!interval_cache_valid) {
        vector<seqitem_t> &cent = centroid.getCentroid();
        cached_interval.setInterval(TimeSeriesInterval(cent.data(), TimeInterval(0, cent.size() - 1)));
        interval_cache_valid = true;
    }

    return cached_interval;
}

vector<seqitem_t> TimeSeriesSliceGroup::getCentroid(void)
{
    return centroid.getCentroid();
}

void TimeSeriesSliceGroup::toFile(ostream *out)
{
    (*out) << centroid.getCount();
    for (int i = 0; i < slice->seqCount(); i++)
        (*out) << ((member.get(i))? "1 " : "0 ");

    vector<seqitem_t> centroid = getCentroid();
    for (unsigned int i = 0; i < centroid.size(); i++)
        (*out) << centroid[i] << " ";

    (*out) << endl;
}

void TimeSeriesSliceGroup::fromFile(istream *in)
{
    int count;
    (*in) >> count;
    bool t;

    for (int i = 0; i < slice->seqCount(); i++) {
        (*in) >> t;
        if (t)
            member.set(i);
    }

    vector<seqitem_t> centroid(count, 0.0);
    for (unsigned int i = 0; i < centroid.size(); i++)
        (*in) >> centroid[i];

    this->centroid = TimeSeriesIntervalCentroid(centroid, count);
}

TimeSeriesSliceGrouping::TimeSeriesSliceGrouping(TimeSeriesSetSlice setslice):slice(setslice)
{
}
seqitem_t TimeSeriesSliceGroup::distance(int index, SeriesDistanceMetric *metric, seqitem_t dropout)
{
    vector<seqitem_t> &cent = centroid.getCentroid();
    TimeSeriesInterval a = (*slice)[index];
    TimeSeriesInterval b = TimeSeriesInterval((seqitem_t*) cent.data(),
                                              TimeInterval(0, slice->length() - 1));

    return metric->run(a,b,dropout);
}

seqitem_t TimeSeriesSliceGroup::distance(vector<seqitem_t> data, SeriesDistanceMetric *metric, seqitem_t dropout)
{
    vector<seqitem_t> &cent = centroid.getCentroid();
    TimeSeriesInterval a = TimeSeriesInterval((seqitem_t*) data.data(), TimeInterval(0, data.size() - 1));
    TimeSeriesInterval b = TimeSeriesInterval((seqitem_t*) cent.data(),
                                              TimeInterval(0, slice->length() - 1));

    seqitem_t t = metric->run(a,b,dropout);
    return t;
}

void TimeSeriesSliceGrouping::hgroup(TimeSeriesSliceGrouping *left,
                                     TimeSeriesSliceGrouping *right)
{
    // Clear previous groupings, if any.
    groups.clear();
    centroidTotalDiffs.clear();

    // Iterate over all group pairs in left, right.
    for (unsigned int li = 0; li < left->groups.size(); li++) {
        for (unsigned int ri = 0; ri < right->groups.size(); ri++) {

            vector<int> intersect = left->groups[li]->intersect(right->groups[ri]);

            if (intersect.size() < 2)
                continue;

            TimeSeriesSliceGroup *group = new TimeSeriesSliceGroup(&slice);

            for (unsigned int j = 0; j < intersect.size(); j++)
                group->addMember(intersect[j], true);

            groups.push_back(group);
        }
    }

    updateDiffs();
}

void TimeSeriesSliceGrouping::naiveGroup(seqitem_t ST)
{
    // Clear previous grouping, if any.
    groups.clear();
    centroidTotalDiffs.clear();

    // For each interval to group, insert into first >ST/2 group, or make a new one.
    for (int si = 0; si < slice.seqCount(); si++) {
        unsigned int gi;
        for (gi = 0; gi < groups.size(); gi++) {
            if (groups[gi]->distance(si) <= ST/2)
                break;
        }
        if (gi == groups.size()) {// If no matches, add new.
            groups.push_back(new TimeSeriesSliceGroup(&slice));
        }
        groups[gi]->addMember(si);
    }

    updateDiffs();
}

void TimeSeriesSliceGrouping::updateDiffs(void)
{
    maxDiff = -INF;
    for (unsigned int i = 0; i < groups.size(); i++) {
        double sumDiff = 0;
        for (unsigned int j = 0; j < groups.size(); j++) {
            vector<double> t = groups[j]->getCentroid();
            sumDiff += groups[i]->distance(t);
        }
        centroidTotalDiffs.push_back(sumDiff);
        maxDiff = max(maxDiff, sumDiff);
    }
}

int TimeSeriesSliceGrouping::getBestGroup(TimeSeriesIntervalEnvelope query, double *dist, double dropout)
{
    double bsf = dropout;
    int bsfIndex = -1;
    for (unsigned int i = 0; i < groups.size(); i++) {
        double dl = groups[i]->getEnvelope().cascadeDist(query, bsf);
        if (dl < bsf) {
            bsf = dl;
            bsfIndex = i;
        }
    }

    *dist = bsf;

    return bsfIndex;
}

void TimeSeriesSliceGrouping::free(void) {
    groups.clear();
}


void TimeSeriesSliceGrouping::toFile(ostream *out)
{
    (*out) << groups.size() << " " << slice.length() << endl;
    for (unsigned int i = 0; i < groups.size(); i++)
        groups[i]->toFile(out);
}

void TimeSeriesSliceGrouping::fromFile(istream *in)
{
    groups.clear();
    int gcount, slength;
    (*in) >> gcount >> slength;
    for (int i = 0; i < gcount; i++) {
        groups.push_back(new TimeSeriesSliceGroup(&slice));
        groups[i]->fromFile(in);
    }
}

TimeSeriesGrouping::TimeSeriesGrouping(TimeSeriesSet &database, seqitem_t ST):set(database)
{
    this->ST = ST;

    int slength = set.seqLength;
    int intervals = slength * slength;
    groups = (TimeSeriesSliceGrouping*) calloc(intervals, sizeof(TimeSeriesSliceGrouping));

    if (groups != NULL) {
        for (int s_t = 0; s_t < database.seqLength; s_t++) {
            for (int e_t = s_t; e_t < database.seqLength; e_t++) {
                TimeSeriesSetSlice dbslice = TimeSeriesSetSlice(&set, TimeInterval(s_t, e_t));
                groups[s_t * slength + e_t] = TimeSeriesSliceGrouping(dbslice);
            }
        }
    }
}

TimeSeriesGrouping::~TimeSeriesGrouping(void)
{
    int slength = set.seqLength;

    if (groups != NULL) {
        for (int s_t = 0; s_t < set.seqLength; s_t ++) {
            for (int e_t = s_t; e_t < set.seqLength; e_t++) {
                groups[s_t * slength + e_t].free();
            }
        }
        free(groups);
    }
}


TimeSeriesSliceGrouping &TimeSeriesGrouping::getGroup(TimeInterval slice)
{
    return groups[slice.start * set.seqLength + slice.end];
}

void TimeSeriesGrouping::naiveGroup(bool debug)
{
    int slength = set.seqLength;
    int count = 0;

    if (groups != NULL) {
        for (int s_t = 0; s_t < slength; s_t ++) {
            if (debug)
                cout << "Grouping intervals from start=" << s_t << endl;
            for (int e_t = s_t; e_t < slength; e_t++) {
                groups[s_t * slength + e_t].naiveGroup(ST);
                count += groups[s_t * slength + e_t].groups.size();
            }
        }
        cout << "Finished grouping with a total of " << count << " groups over " << (slength*(slength-1))/2
             << " discrete intervals." << endl;
    }
}

void TimeSeriesGrouping::hgroup(bool debug)
{
    int slength = set.seqLength;
    int count = 0;

    if (groups != NULL) {
        cout << "Grouping base intervals." << endl;
        for (int i = 0; i < slength; i++) {
            groups[i * slength + i].naiveGroup(ST);
            count += groups[i * slength + i].groups.size();
        }

        cout << "Grouping intersected intervals." << endl;

        for (int s = 0; s < slength - 1; s++) {
            if (debug)
                cout << "Grouping intervals from start=" << s << endl;
            for (int e = s + 1; e < slength; e++) {
                groups[s * slength + e].hgroup(&groups[s * slength + e - 1],
                                                   &groups[e * slength + e]);
                count += groups[s * slength + e].groups.size();
            }
        }

        cout << "Finished grouping with a total of " << count << " groups over " << (slength*(slength+1))/2
             << " discrete intervals." << endl;
        /*
        for (int e = 1; e < slength; e++)
            groups[e].hgroup(&groups[e-1],
                             &groups[e*slength + e]);
        */
    }
}



bool TimeSeriesGrouping::valid(void)
{
    return groups != NULL;
}

seqitem_t TimeSeriesGrouping::getST(void)
{
    return ST;
}

void TimeSeriesGrouping::setST(seqitem_t ST)
{
    this->ST = ST;
}

int TimeSeriesGrouping::fromFile(const char *path)
{
    ifstream in;
    in.open(path);

    if (!in.is_open())
        return -1;

    int slen;
    in >> slen >> ST;
    int s,e;
    for (int s_t = 0; s_t < slen; s_t++) {
        for (int e_t = s_t; e_t < slen; e_t++) {
            in >> s >> e;
            groups[s * slen + e].fromFile(&in);
        }
    }

    in.close();

    return 0;
}

int TimeSeriesGrouping::toFile(const char *path)
{
    ofstream out;

    out.open(path);
    if (!out.is_open())
        return -1;

    int slen = set.seqLength;
    out << slen << " " << ST << endl;
    for (int s_t = 0; s_t < slen; s_t++) {
        for (int e_t = s_t; e_t < slen; e_t++) {
            out << s_t << " " << e_t << endl;
            groups[s_t * slen + e_t].toFile(&out);
        }
    }

    out.close();

    return 0;
}
