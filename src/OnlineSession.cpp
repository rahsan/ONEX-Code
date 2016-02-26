#include "OnlineSession.h"

#include <map>
#include <string>
#include <exception>
#include <stdexcept>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

OnlineSession::OnlineSession(double defaultST)
{
    this->defaultST = defaultST;
}

OnlineSession::~OnlineSession(void)
{
    for (unsigned int i = 0; i < dataSets.size(); i++) {
        if (groupings[i])
            delete groupings[i];
        if (dataSets[i])
            delete dataSets[i];
    }
}

int OnlineSession::loadDataSet(TimeSeriesSet *data)
{
    if (data == NULL)
        return -1;

    if (!data->valid())
        return -1;

    int index = dataSets.size();

    dataSets.push_back(data);
    groupings.push_back(NULL);

    return index;
}

void OnlineSession::checkIndex(int index)
{
    if ((unsigned) index >= dataSets.size())
        throw out_of_range("No dataset with that index.");
}

int OnlineSession::saveDataSet(int index, const char *path, bool old)
{
    checkIndex(index);

    return dataSets[index]->toFile(path, !old);
}

int OnlineSession::dropDataSet(int index)
{
    checkIndex(index);

    dataSets.erase(dataSets.begin() + index);
    groupings.erase(groupings.begin() + index);
    return 0;
}

void OnlineSession::genGrouping(int index, double ST)
{
    checkIndex(index);

    if (groupings[index] != NULL)
        delete groupings[index];

    cout << "Initializing grouping datastructure." << endl;
    groupings[index] = new TimeSeriesGrouping(*dataSets[index], ST);
    groupings[index]->setST(ST);
    cout << "Grouping intervals." << endl;
    groupings[index]->hgroup(debug);
}

void OnlineSession::genNaiveGrouping(int index, double ST)
{
    checkIndex(index);

    if (groupings[index] != NULL)
        delete groupings[index];

    groupings[index] = new TimeSeriesGrouping(*dataSets[index], ST);
    groupings[index]->setST(ST);
    groupings[index]->naiveGroup(debug);
}

int OnlineSession::loadGrouping(int index, const char *path)
{
    checkIndex(index);

    if (groupings[index] != NULL)
        delete groupings[index];
        
    groupings[index] = new TimeSeriesGrouping(*dataSets[index], defaultST);
    return groupings[index]->fromFile(path);
}

int OnlineSession::saveGrouping(int index, const char *path)
{
    checkIndex(index);

    if (groupings[index] == NULL) {
        cerr << "Dataset has not been grouped." << endl;
        return -1;
    }

    return groupings[index]->toFile(path);
}

void OnlineSession::dropGrouping(int index)
{
    checkIndex(index);

    if (groupings[index] != NULL)
        delete groupings[index];

    groupings[index] = NULL;
}

int OnlineSession::normalizeDataSet(int index)
{
    checkIndex(index);

    dataSets[index]->normalize();

    return 0;
}

void OnlineSession::printDataSets(ostream *out)
{
    *out << "Available data sets:" << endl;
    for (unsigned int i = 0; i < dataSets.size(); i++) {
        *out << "[" << i << "] " << ((groupings[i] == NULL)? "UNGROUPED " : "  GROUPED ");
        *out << dataSets[i]->name << "." << endl;
    }
}

enum _commands {

    _LOAD_TSS = 1,   // Load a time series set.
    _OLOAD_TSS,  // Load a time series set from old file format.
    _SAVE_TSS,   // Save a time series set.
    _OSAVE_TSS,  // Save a time series set to old file format.
    _DROP_TSS,   // Drop a time series set.
    _RAND_TSS,   // Generate a random time series set.

    _LIST_TSS,   // List available time series set.

    _NORM_TSS,      // Normalize a time series set.

    _KSIM_TSS,      // Get the k similar neighbors.
    _OUTLIER_TSS,   // Find the dominant outlier.

    _GROUP_TSS,      // Group a time series set with default ST.
    _GROUP_ST_TSS,   // Group a time series set with given ST.
    _NGROUP_TSS,     // Group a time series set with default ST.
    _NGROUP_ST_TSS,  // Group a time series set with default ST.
    _LOAD_GROUP_TSS, // Load a grouping file.
    _SAVE_GROUP_TSS, // Save a grouping file.
    _DROP_GROUP_TSS, // Drop a grouping.

    _DESC_TSS,   // Describe a time series set.
    _PRINT_TSS,  // Print a time series set.

    _SET_DEF_ST, // Set the default ST.
    _GET_DEF_ST, // Get the default ST.

    _DEBUG,      // Enter debug mode (prints more.)

    _HELP,       // Print help.
    _EXIT        // Exit the program.
};

void _print_help(ostream *out)
{
    *out << "ONEX interactive command line interface." << endl;
    *out << "Available commands:" << endl;
    *out << "load <path>          Load a dataset from the given file. New format only." << endl;
    *out << "save <id> <path>     Save a loaded dataset to the given file. New format." << endl;
    *out << "drop | rm <id>       Unload a loaded database.                           " << endl;
    *out << "random <N> <L> <rng> Generate a random NxL dataset with range rng.       " << endl;
    *out << endl;
    *out << "normalize <id>       Normalize the given dataset.                        " << endl;
    *out << endl;
    *out << "group <id>           (Re)Group the given dataset using current ST.           " << endl;
    *out << "ngroup <id>          (Re)Group the given dataset using current ST (naive).   " << endl;
    *out << "groupST <id> <ST>    (Re)Group the given dataset using the given ST.         " << endl;
    *out << "ngroupST <id> <ST>   (Re)Group the given dataset using the given ST (naive). " << endl;
    *out << "gload <id> <path>    Load grouping data for the given dataset from file.     " << endl;
    *out << "gsave <id> <path>    Save grouping data for the given dataset to file.       " << endl;
    *out << "gdrop <id> |grm <id> Delete grouping data for the given dataset.             " << endl;
    *out << endl;
    *out << "ksim <dbid> <qid> <qN> <start> <end> <k> <strict>                        " << endl;
    *out << "                     Get the kSimilar sequences from dataset dbid to the " << endl;
    *out << "                     query string in dataset qid at qN from start to end." << endl;
    *out << "                     Print k results and if strict, only check the same  " << endl;
    *out << "                     intervals.                                          " << endl;
    *out << "outlier <id> <s> <e> Find a dominant outlier in the dataset.             " << endl;
    *out << endl;
    *out << "describe | desc <id> Describe a dataset.                                 " << endl;
    *out << "print <id>           Print a dataset.                                    " << endl;
    *out << "list | ls            List all loaded / available datasets.               " << endl;
    *out << endl;
    *out << "setST <value>        Set the default Similarity Threshold.               " << endl;
    *out << "getST                View the current default Similarity Threshold.      " << endl;
    *out << endl;
    *out << "debug <value>        Set the verbosity. 0=quiet, 1=verbose.              " << endl;
    *out << endl;
    *out << "osave <id> <path>         Save a loaded dataset to the given file. Old format.            " << endl;
    *out << "oload <path> <N> <L> <D>  Load a dataset from a raw file, skip first D in each line.      " << endl;
    *out << endl;
    *out << "help | h | ?         Print this message.                                 " << endl;
    *out << "exit|quit|q          Exit from ONEX CLI.                                 " << endl;
}

void prepareCommands(map<string, int> &commands) {
    
    commands["exit"] = _EXIT;
    commands["quit"] = _EXIT;
    commands["q"] = _EXIT;

    commands["load"] = _LOAD_TSS;
    commands["oload"] = _OLOAD_TSS;
    commands["save"] = _SAVE_TSS;
    commands["osave"] = _OSAVE_TSS;
    commands["rm"] = _DROP_TSS;
    commands["drop"] = _DROP_TSS;
    commands["random"] = _RAND_TSS;

    commands["list"] = _LIST_TSS;
    commands["ls"] = _LIST_TSS;
    commands["help"] = _HELP;
    commands["h"] = _HELP;
    commands["?"] = _HELP;

    commands["desc"] = _DESC_TSS;
    commands["describe"] = _DESC_TSS;
    commands["print"] = _PRINT_TSS;

    commands["ksim"] = _KSIM_TSS;
    commands["outlier"] = _OUTLIER_TSS;
    
    commands["normalize"] = _NORM_TSS;

    commands["group"] = _GROUP_TSS;
    commands["groupST"] = _GROUP_ST_TSS;
    commands["ngroup"] = _NGROUP_TSS;
    commands["ngroupST"] = _NGROUP_ST_TSS;

    commands["gload"] = _LOAD_GROUP_TSS;
    commands["gsave"] = _SAVE_GROUP_TSS;
    commands["gdrop"] = _DROP_GROUP_TSS;
    commands["grm"] = _DROP_GROUP_TSS;

    commands["setST"] = _SET_DEF_ST;
    commands["getST"] = _GET_DEF_ST;

    commands["debug"] = _DEBUG;
}

int OnlineSession::run(istream *in, ostream *out, bool interactive)
{
    map<string, int> commands;
    prepareCommands(commands);

    int res = 0;
    int count = 0;

    int cmd;
    string command;

    string sarg1, sarg2;
    int iarg1, iarg2, iarg3, iarg4, iarg5, iarg6, iarg7;
    double darg1;
    double *dparg1;
    TimeSeriesSet *t;

    clock_t time;

    res = 0;

    while (!in->eof()) {

        if (res != 0)
            *out << "Command returned with status " << res << "." << endl;

        if (interactive) {
            int width = out->width();
            *out << "[";

            out->width(3);
            *out << count;
            out->width(width);

            *out << "] > ";
            out->flush();
        }

        *in >> command;
        time = clock();

        if (command.size() > 0) {
            cmd = commands[command];
        }

        if (command.size() <= 0 || in->eof()) {
            cmd = _EXIT;
            *out << endl;
        }

        if (cmd == 0) {
            *out << "Unknown command '" << command << "'. Type 'help' for help." << endl;
            res = 1;
            count++;
            continue;
        }

        try {
            res = 0;
            switch (cmd) {
            case _EXIT:
                *out << "Quitting..." << endl;

                return iarg1;

            case _HELP:
                _print_help(out);

                break;

            case _LOAD_TSS:
                *in >> sarg1;

                *out << "Loading Time Series Set from file '" << sarg1 << "'." << endl;

                res = loadDataSet(new TimeSeriesSet(sarg1.c_str()));
                if (res >= 0) {
                    *out << "Dataset successfully loaded. Index: " << res << endl;
                    res = 0;
                } else {
                    *out << "Failed to load dataset." << endl;
                }

                break;

            case _SAVE_TSS:
                *in >> iarg1;
                *in >> sarg2;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << "Saving Time Series Set " << iarg1 << ":" << t->name << " to file '" << sarg2 << "'." << endl;

                res = saveDataSet(iarg1, sarg2.c_str());
                if (res == 0)
                    *out << "Dataset successfully saved." << endl;

                break;

            case _OLOAD_TSS:
                *in >> sarg1;
                *in >> iarg1 >> iarg2 >> iarg3;

                *out << "Loading Time Series Set from file '" << sarg1 << "' with N=" << iarg1 << ", L=" << iarg2 << ", and D=" << iarg3 << "." << endl;

                res = loadDataSet(new TimeSeriesSet(sarg1.c_str(), iarg1, iarg2, iarg3));
                if (res >= 0) {
                    *out << "Dataset successfully loaded. Index: " << res << endl;
                    res = 0;
                } else {
                    *out << "Failed to load dataset." << endl;
                }

                break;

            case _OSAVE_TSS:
                *in >> iarg1;
                *in >> sarg2;

                checkIndex(iarg1);
                t = dataSets[iarg1];

                *out << "Saving Time Series Set " << iarg1 << ":" << t->name << " to file '" << sarg2 << "'." << endl;
                res = saveDataSet(iarg1, sarg2.c_str(), true);
                if (res == 0)
                    *out << "Dataset successfully saved." << endl;
                break;

            case _DROP_TSS:
                *in >> iarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << "Dropping Time Series Set " << iarg1 << ":" << t->name << "." << endl;

                res = dropDataSet(iarg1);

                break;

            case _RAND_TSS:
                *in >> iarg1 >> iarg2 >> iarg3;

                *out << "Generating random Time Series Set with N=" << iarg1 << ", L=" << iarg2 << ", and range=" << iarg3 << "." << endl;

                res = loadDataSet(&TimeSeriesSet::randomSet(iarg1, iarg2, iarg3));
                if (res >= 0) {
                    *out << "Dataset successfully loaded. Index: " << res << endl;
                    res = 0;
                } else {
                    *out << "Failed to load dataset." << endl;
                }

                break;

            case _LIST_TSS:
                printDataSets(out);
                break;

            case _NORM_TSS:
                *in >> iarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << "Normalizing Time Series Set " << iarg1 << ":" << t->name << "." << endl;

                res = normalizeDataSet(iarg1);

                break;

            case _KSIM_TSS:
                *in >> iarg1 >> iarg2 >> iarg3 >> iarg4 >> iarg5 >> iarg6 >> iarg7;

                checkIndex(iarg1);
                checkIndex(iarg2);
                t = dataSets[iarg1];
                *out << "Getting the kSimilar for Time Series Set " << iarg1 << ":" << t->name;
                t = dataSets[iarg2];
                *out << ", query string at " << iarg3 << " in dataset " << iarg2 << ":" << t->name;
                *out << " in interval [" << iarg4 << "," << iarg5 << "] with k=" << iarg6 << " and strict="
                     << iarg7 << "." << endl;

                dparg1 = t->getInterval(iarg3,TimeInterval(iarg4, iarg5)).getData();
                kSimilar(iarg1, vector<double>(dparg1, dparg1 + (iarg5 - iarg4 + 1)),
                         TimeInterval(iarg4, iarg5), iarg6, iarg7);

                break;

            case _OUTLIER_TSS:
                *in >> iarg1 >> iarg2 >> iarg3;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << "Gettind dominant outlier on Time Series Set " << iarg1 << ":" << t->name
                     << " in the range [" << iarg2 << ", " << iarg3 << "]." << endl;

                dominantOutlier(iarg1, TimeInterval(iarg2, iarg3));

                break;

            case _GROUP_TSS:
                *in >> iarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];

                *out << ((groupings[iarg1] != NULL)?"Regrouping":"Grouping") << " Time Series Set "
                     << iarg1 << ":" << t->name << " with ST " << defaultST << "." << endl;

                genGrouping(iarg1, defaultST);

                break;

            case _NGROUP_TSS:
                *in >> iarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];

                *out << ((groupings[iarg1] != NULL)?"Regrouping":"Grouping") << " Time Series Set "
                     << iarg1 << ":" << t->name << " with ST " << defaultST << " (naive)." << endl;

                genNaiveGrouping(iarg1, defaultST);

                break;

            case _GROUP_ST_TSS:
                *in >> iarg1;
                *in >> darg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << ((groupings[iarg1]  == NULL)?"Grouping":"Regrouping") << " Time Series Set "
                     << iarg1 << ":" << t->name << " with ST " << darg1 << "." << endl;

                genGrouping(iarg1, darg1);

                break;

            case _NGROUP_ST_TSS:
                *in >> iarg1;
                *in >> darg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << ((groupings[iarg1]  == NULL)?"Grouping":"Regrouping") << " Time Series Set "
                     << iarg1 << ":" << t->name << " with ST " << darg1 << " (naive)." << endl;

                genNaiveGrouping(iarg1, darg1);

                break;

            case _LOAD_GROUP_TSS:
                *in >> iarg1;
                *in >> sarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << "Loading grouping data for Time Series Set "
                     << iarg1 << ":" << t->name << " at " << sarg1 << "." << endl;

                res = loadGrouping(iarg1, sarg1.c_str());

                break;

            case _SAVE_GROUP_TSS:
                *in >> iarg1;
                *in >> sarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << "Saving grouping data for Time Series Set "
                     << iarg1 << ":" << t->name << " to " << sarg1 << "." << endl;

                res = saveGrouping(iarg1, sarg1.c_str());

                break;

            case _DROP_GROUP_TSS:
                *in >> iarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                *out << "Dropping grouping data for Time Series Set "
                     << iarg1 << ":" << t->name << "." << endl;

                dropGrouping(iarg1);

                break;

            case _DESC_TSS:
                *in >> iarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                t->printDesc(out);

                break;

            case _PRINT_TSS:
                *in >> iarg1;

                checkIndex(iarg1);
                t = dataSets[iarg1];
                t->printData(out);

                break;

            case _DEBUG:
                *in >> iarg1;

                *out << "Setting debug to " << iarg1 << "." << endl;
                debug = iarg1;

                break;

            case _SET_DEF_ST:
                *in >> darg1;

                *out << "Setting default ST to " << darg1 << "." << endl;

                defaultST = darg1;

                break;

            case _GET_DEF_ST:
                *out << "default ST = " << defaultST << "." << endl;

                break;
            }

        } catch (exception &e) {
            *out << "Caught exception attempting operation:" << endl;
            *out << e.what() << endl;
        }

        time = clock() - time;
        (*out) << "Command used " << ((float)time/CLOCKS_PER_SEC) << " seconds." << endl << endl;

        count++;
    }

    return 0;
}

void OnlineSession::kSimilar(int dbindex, vector<double> qdata, TimeInterval interval, int k, int strict)
{
    if (groupings[dbindex] == NULL)
        genGrouping(dbindex, defaultST);

    TimeSeriesGrouping *t = groupings[dbindex];

    int slen = dataSets[dbindex]->seqLength;
    int ilen = interval.length();

    TimeSeriesIntervalEnvelope eqdata(TimeSeriesInterval(qdata.data(), TimeInterval(0, qdata.size() - 1)));

    double gb;
    int bsfStart = 0, bsfEnd = 0, bsfIndex = -1;
    double bsf = INF;

    if (strict == 0) {
        for (int i = 0; i < slen; i++) {
            if (debug) {
                cout << "Searching groups from start=" << i << ", bsf=" << bsf << endl;
            }
            for (int j = i; j < slen; j++) {
                int k = t->groups[i * slen + j].getBestGroup(eqdata, &gb, bsf);
                if (k < 0)
                    continue;
                if (gb < bsf) {
                    bsf = gb;
                    bsfStart = i;
                    bsfEnd = j;
                    bsfIndex = k;
                }
            }
        }
    } else if (strict == 1) {
        for (int i = 0; i + ilen - 1 < slen; i++) {
            if (debug)
                cout << "Searching groups from start=" << i << endl;
            int j = i + (ilen - 1);
            int k = t->groups[i * slen + j].getBestGroup(eqdata, &gb, bsf);
            if (k < 0)
                continue;
            if (gb < bsf) {
                bsf = gb;
                bsfStart = i;
                bsfEnd = j;
                bsfIndex = k;
            }
        }
    } else {
        bsfStart = interval.start;
        bsfEnd = interval.end;
        bsfIndex = t->groups[bsfStart * slen + bsfEnd].getBestGroup(eqdata, &gb, bsf);
        if (bsfIndex >= 0)
            bsf = gb;
    }

    if (bsf == INF || bsfIndex < 0) {
        cerr << "kSimilar: Failed to find similar objects. No suitable candidate group centroids." << endl;
        return;
    }

    cout << "Found most similar interval and group: " << bsfIndex << "@"
         << "[" << bsfStart << "," << bsfEnd << "]" << endl;

    vector<kSim> sim = t->groups[bsfStart * slen + bsfEnd].groups[bsfIndex]->getSortedSimilar(eqdata, k);

    cout << "Discovered k similar points:" << endl;
    for (unsigned int i = 0; i < sim.size(); i++) {

        cout << "Series " << sim[i].index << ", interval [" << bsfStart << "," << bsfEnd
             << "] is at distance " << sim[i].distance << "." << endl;

        TimeSeriesInterval interval = (*t->groups[bsfStart * slen + bsfEnd].groups[bsfIndex]->slice)[sim[i].index];

        for (int j = 0; j < interval.length(); j++) {
            cout << interval[j] << " ";
        }
        cout << endl;
    }
}

// void OnlineSession::kSimilar(int dbindex, int qindex, TimeInterval interval, bool sameInterval)
// {
//     if (groupings[dbindex] == NULL)
//         genGroupings(dbindex, defaultST);

//     double bsf=INF;

//     int bsfIndex;       //index of bsf centroid
//     int bsfTSIndex;
//     double currentDist=0;

//     if (sameInterval==false)    //ANY time interval
//     {
//         //call that embedded code

//         //compare the query with all the centroids
//         for(unsigned int i=0;i<groupArray.size();i++)
//         {
//             currentDist=simpleDTW(groupArray[i].centroid,tempQ);
//             if(currentDist<bsf)
//             {
//                 bsf=currentDist;
//                 bsfIndex=i;
//             }
//         }
//     }
//     else
//     {
//         //EXACT match
//         //iterate through the groups in the time interval provided
//         int numGroups=Time[startTime][endTime].count;
//         for(int i=0;i<numGroups;i++)
//         {
//             int groupID=Time[startTime][endTime].groupid[i];
//             currentDist=simpleDTW(tempQ,groupArray[groupID].centroid);
//             if(currentDist<bsf)
//             {
//                 bsf=currentDist;
//                 bsfIndex=groupID;
//             }
//         }

//     }
//     //now we have the closest centroid
//     bsf=INF;
//     vector<double> tempSeries;
//     bsfTSIndex=0;
//     int kbestCount=0;
//     kBest tempBest;

//     if(k==1)
//     {
//         //compare the DTW between all the time series in this group
//         for(int i=0;i<N;i++)
//         {
//             if(groupArray[bsfIndex].bitVector[i]==1)    //this time series is present in the group
//             {
//                 //copy this time series in tempSeries
//                 for(int j=0;j<groupArray[bsfIndex].endT-groupArray[bsfIndex].startT+1;j++)
//                 {
//                     tempSeries.push_back(timeSeries[i][groupArray[bsfIndex].startT+j]);
//                 }
//                 currentDist=simpleDTW(tempSeries,tempQ);
//                 if(currentDist<bsf)
//                 {
//                     //update the best index and bsf
//                     bsf=currentDist;
//                     bsfTSIndex=i;
//                 }
//                 tempSeries.clear();
//             }
//         }
//         //have the most similar time series
//         cout<<"Most similar TS is "<<bsfTSIndex<<" Interval "<<groupArray[bsfIndex].startT<<" "<<groupArray[bsfIndex].endT<<" "<<bsf<<endl;
//     }
//     else
//     {
//         //need to find k most similar from this group
//         if(groupArray[bsfIndex].count==k)   //
//         {
//             cout<<"K most similar TS are";
//             for(int i=0;i<N;i++)
//             {
//                 if(groupArray[bsfIndex].bitVector[i]==1)    //this time series is present in the group
//                 {
//                     cout<<i<<" ";
//                 }
//             }
//         }
//         else if(groupArray[bsfIndex].count>k)       //
//         {
//             for(int i=0;i<N;i++)
//             {
//                 if(groupArray[bsfIndex].bitVector[i]==1)    //this time series is present in the group
//                 {
//                     //copy this time series in tempSeries
//                     for(int j=0;j<groupArray[bsfIndex].endT-groupArray[bsfIndex].startT+1;j++)
//                     {
//                         tempSeries.push_back(timeSeries[i][groupArray[bsfIndex].startT+j]);
//                     }
//                     currentDist=simpleDTW(tempQ,tempSeries);
//                     if(kbestCount<k)
//                     {
//                         //add this TS to k best
//                         kbestCount++;
//                         tempBest.dist=currentDist;
//                         tempBest.id=i;
//                         kbestArray.push_back(tempBest);

//                     }
//                     else
//                     {
//                         sort(kbestArray.begin(),kbestArray.end(),sortByDist);
//                         double tempD=kbestArray[kbestCount].dist;    //getting the last distance
//                         if(tempD>currentDist)
//                         {
//                             tempBest.dist=currentDist;
//                             tempBest.id=i;
//                             kbestArray[kbestCount]=tempBest;

//                         }
//                     }
//                     tempSeries.clear();
//                 }
//             }
//             cout<<"K most similar are";
//             for(int i=0;i<k;i++)
//                 cout<<kbestArray[i].id<<" ";
//             cout<<endl;

//         }
//         else
//         {
//             //count of TS in group is < k
//             //finds second best group
//         }
//     }



// }


//finds dominant or outlier TS
//For outlier the time interval is specified
void OnlineSession::dominantOutlier(int dbindex, TimeInterval interval)
{
    if (groupings[dbindex] == NULL)
        genGrouping(dbindex, defaultST);
    TimeSeriesGrouping *gp = groupings[dbindex];
    TimeSeriesSliceGrouping &sgp = gp->getGroup(interval);

    unsigned int i = 0;
    for (i = 0; i < sgp.groups.size(); i++) {
        if (sgp.centroidTotalDiffs[i] == sgp.maxDiff) {
            if (sgp.groups[i]->count == 1) {
                break;
            }
        }
    }

    if (i == sgp.groups.size()) {
        cerr << "Failed to find an outlier: The farthest group has more than one member." << endl;
        return;
    }

    for (unsigned int j = 0; j < sgp.groups.size(); j++) {
        if (sgp.groups[i]->isMember(j))
            cout << "Dominant outlier found: Group " << j << "." << endl;
    }

    return;
}

