#ifndef ALGO_H
#define ALGO_H

#include "TimeSeries.h"

#include <vector>

using namespace std;

void kSimilar(bool queryFlag, int k,bool timeFlag,int startTime, int endTime);
int dominantOutlier(int startTime, int endTime);
void criticalTime(bool queryFlag, double epsilon);
void stableTimeInterval(bool queryFlag, double epsilon);
void evolvingSimilarity(bool queryFlag);
void simThresholdRec(bool simDegreeFlag, char category);

#endif // ALGO_H
