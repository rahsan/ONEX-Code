#include "Groups.h"

#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;

GroupOperation::GroupOperation(double STinput)
{
    //myGroup->seqCentroidDist[1]
    groupCounter=0;
    ST=STinput;
    timeFile.open("data/Time.txt");
    groupFile.open("data/Groups.txt");

}

//reads the time series data from input file and store in structure
//n is number of time series and L is length of time series
int GroupOperation::readFile(char * Filename, int number, int length)
{
    inputFile.open(Filename); //Opening the file
    N=number;
    L=length;
    min=10000;
    max=0;
    if(!inputFile)
    {
        cout << "there's something wrong with the file!\n"; //Throw exception
        return -1;
    }
    else
    {
        //read the input file into time series structure
        //Read numbers from file into the array
        for(int i = 0; i < N; i++)
        {
            for(int j=0; j<L;j++)
            {
                inputFile>>timeSeries[i][j];
                if(timeSeries[i][j]<min)
                    min=timeSeries[i][j];
                if(timeSeries[i][j]>max)
                    max=timeSeries[i][j];
             //   cout<<timeSeries[i][j]<<" ";

            }
        }
        inputFile.close();
        normalize();
        return 1;
    }



}
//normalizes the data such that each value is now between 0 and 1
void GroupOperation::normalize()
{
    cout<<"Min "<<min<<" Max "<<max<<endl;
    for(int i = 0; i < N; i++)
    {
        for(int j=0; j<L;j++)
        {
            timeSeries[i][j]=(timeSeries[i][j]-min)/(max-min);

        }
    }

}
 //prints the time series
void GroupOperation::printTS()
{
    //cout<<"Printing Timeseries"<<endl;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<L;j++)
        {
            cout<<timeSeries[i][j]<<" ";
        }
        cout<<endl;
    }

}
//prints all the times and asociated information (group ids, count of groups)
void GroupOperation::printTime(int level)
{

    //cout<<endl<<endl;
    for (int l=0, m=l+level; m<L; l++, m++)
    {

   //     cout<<"Time intervals: "<<l<<" "<<m<<endl;
        timeFile<<"T"<<" "<<l<<" "<<m<<" "<<Time[l][m].count<<endl;
   //     cout<<"Group Count: "<<Time[l][m].count<<" "<<endl;

            for(int k=0;k<Time[l][m].count;k++)
            {
     //           cout<<"Groups ID: "<<Time[l][m].groupid[k];
                timeFile<<"G"<<" "<<Time[l][m].groupid[k]<<" ";
       //         cout<<endl;
            }
            timeFile<<endl;
            timeFile<<"DcDist"<<" ";
            for(int k=0;k<Time[l][m].count;k++)
            {

                timeFile<<Time[l][m].DcDistance[k]<<" ";
       //         cout<<endl;
            }
            timeFile<<endl;

     }


}
//prints alll groups with associated information
void GroupOperation::printGroups()
{
    cout<<"Total Groups "<<groupArray.size()<<endl;
    groupFile<<"Total"<<" "<<groupArray.size()<<endl;
    groupFile<<"Min "<<min<<" "<<"Max "<<max<<endl;
    for(int i=0;i<groupArray.size();i++)
    {
 //       cout<<"Group "<<groupArray[i].id<<endl;
        groupFile<<"G"<<" "<<groupArray[i].id<<" "<<groupArray[i].count<<endl;
        groupFile<<groupArray[i].startT<<" "<<groupArray[i].endT<<endl;
   //     cout<<"TS Count "<<groupArray[i].count<<endl;
     //   cout<<"Time Interval "<<groupArray[i].startT<<" "<<groupArray[i].endT<<endl;

        for(int k=0;k<timeSeriesN;k++)
        {
            //printing the bit vector
       //     cout<<groupArray[i].bitVector[k];
            groupFile<<groupArray[i].bitVector[k]<<" ";
        }
    //    cout<<endl;
        groupFile<<endl;
    //    cout<<"Group Centroid ";
        groupFile<<groupArray[i].centroid.size()<<endl;
        for(int k =0;k<groupArray[i].centroid.size();k++)
        {
      //      cout<<groupArray[i].centroid[k]<<" ";
            groupFile<<groupArray[i].centroid[k]<<" ";
        }
   //     cout<<endl;
        groupFile<<endl;
    }

    groupFile.close();

}

void GroupOperation::printOneGroup(int groupId)
{

    cout<<"Group "<<groupId<<": Attributes"<<endl;
    cout<<"Start T: "<<groupArray[groupId].startT<<endl;
    cout<<"End T: "<<groupArray[groupId].endT<<endl;
    cout<<"TS Count: "<<groupArray[groupId].count<<endl;
}

//groups the time series
void GroupOperation::groupOp()
{
    //group the time series for only base time interval
    for(int i=0,j=1;i<L-1 && j<L;i++,j++)
    {
        groupByTime(i,j);
       // cout<<i<<" "<<j<<endl;
    }
    //do intersections for other time intervals
    for(int j=1;j<L;j++)
        intersectGroup(j);
    //printing all the groups
    printGroups();

    //

    for(int j=1;j<L;j++)
        printTime(j);

    timeFile.close();
}
void GroupOperation::groupByTime(int startTime,int endTime)
{
    Time[startTime][endTime].count=0;
    Group *newGroup;//=new Group();
    double EDDist=0;
    int bestGroupIndex=0;
    double bsf=0;       //best so far distance
    if(Time[startTime][endTime].count==0)   //no groups in this time interval
    {
        newGroup=new Group();

        //there are no groups

        //make the first group

        newGroup->id=groupCounter++;
        //add this group to the the list of groups for this time interval
        Time[startTime][endTime].groupid.push_back(newGroup->id);
        Time[startTime][endTime].count++;

        newGroup->startT=startTime;
        newGroup->endT=endTime;
        newGroup->count=1;
        newGroup->bitVector[0]=1;

        for(int i=startTime;i<=endTime;i++)
        {
            //first new centroid is first time series
            newGroup->centroid.push_back(timeSeries[0][i]);
            newGroup->sum.push_back(timeSeries[0][i]);
        }
        newGroup->centroidLen=newGroup->centroid.size();
        //push this group in groups vector
        groupArray.push_back(*newGroup);
      //  delete newGroup;
    }
    //iterate over all the time series
    for(int i=1;i<N;i++)
    {
        bestGroupIndex=0;
        bsf=10000;       //best so far distance

        //check with all the groups centroid for this time interval
        for(int k=0;k<Time[startTime][endTime].count;k++)
        {
            EDDist=0;
            for(int j=0;j<=endTime-startTime;j++)
            {
                //compare time series with centroid one at a time
                EDDist+=pow((timeSeries[i][j+startTime]-groupArray[Time[startTime][endTime].groupid[k]].centroid[j]),2.0);
            }
            EDDist=sqrt(EDDist);
            if(EDDist<bsf)
            {
                bsf=EDDist;
                bestGroupIndex=Time[startTime][endTime].groupid[k];
            }

        }
   //     cout<<"Best so Far distance "<<bsf<<endl;
   //     cout<<"Best group Index "<<bestGroupIndex<<endl;
        //checking if distance is within ST add in that group
        if(bsf<=ST/2)
        {
            groupArray[bestGroupIndex].count++; //increase the count of TS in this group

            //update the bit vector
            groupArray[bestGroupIndex].bitVector[i]=1;
            //update the centroid

            for(int j=0;j<groupArray[bestGroupIndex].centroidLen;j++)
            {
                //update centroid by calculating mean of each value.
                double tempSum=groupArray[bestGroupIndex].sum[j];
                tempSum+=timeSeries[i][startTime+j];
                //update the sum
                groupArray[bestGroupIndex].sum[j]=tempSum;
                double temp3=tempSum/groupArray[bestGroupIndex].count;
                groupArray[bestGroupIndex].centroid[j]=temp3;//double (groupArray[bestGroupIndex].centroid[j]+timeSeries[i][startTime+j])/groupArray[bestGroupIndex].centroidLen;
                //cout<<"Centroid of length "<<groupArray[bestGroupIndex].centroidLen<<"is "<<groupArray[bestGroupIndex].centroid[j]<<"of "<<groupArray[bestGroupIndex].centroid[j]<<" "<<timeSeries[i][startTime+j]<<endl;
            }

        }
        else
        {
            //make new group
            Group *nGroup=new Group();
            nGroup->id=groupCounter++;
            //add this group id to this time interval groups
            Time[startTime][endTime].groupid.push_back(nGroup->id);
            Time[startTime][endTime].count++;
            nGroup->count=1;
            nGroup->bitVector[i]=1;
            nGroup->startT=startTime;
            nGroup->endT=endTime;
            for(int k=startTime;k<=endTime;k++)
            {
                //first new centroid is first time series
                nGroup->centroid.push_back(timeSeries[i][k]);
                nGroup->sum.push_back(timeSeries[i][k]);
            }
            nGroup->centroidLen=nGroup->centroid.size();
            //push this group in groups vector
            groupArray.push_back(*nGroup);
       //     delete nGroup;

        }

    }
    //for this time interval for each group calculate the sum of ED distances with other centroids
    double DcEDDist=0;
    double tempDcDist=0;
    for(int i=0;i<Time[startTime][endTime].count;i++)
    {
        //pick first group
        int groupId=Time[startTime][endTime].groupid[i];
        //calculate the distance of its centroid with all other
        for(int j=0;j<Time[startTime][endTime].count;j++)
        {

            if(j!=i)
            {
                int groupId2=Time[startTime][endTime].groupid[j];
                for(int k=0;k<groupArray[groupId].centroid.size();k++)
                {
                    tempDcDist+=pow((groupArray[groupId].centroid[k]-groupArray[groupId2].centroid[k]),2.0);

                }
                DcEDDist+=sqrt(tempDcDist);
                tempDcDist=0;
            }
        }
        //save the DCEDDist in the structure
        Time[startTime][endTime].DcDistance.push_back(DcEDDist);
        DcEDDist=0;


    }

}

void GroupOperation::intersectGroup(int level)  //finds the intersection of groups
{

    for (int l=0, m=l+level; m<L-1; l++, m++)
    {

        Time[l][m+1].count=0;
         for (int i=0;i<Time[l][m].count;i++)
             {
                int gri=Time[l][m].groupid[i];
                for (int j=0;j<Time[m][m+1].count;j++)
                {
                  int grj=Time[m][m+1].groupid[j];


                 //make new group
                   Group *nGroup=new Group();

                   nGroup->startT=l;
                   nGroup->endT=m+1;
                   nGroup->count=0;

                   for (int k=0;k<N;k++)
                   {
                        if ((groupArray[gri].bitVector[k]==groupArray[grj].bitVector[k]) &&(groupArray[gri].bitVector[k]==1))
                        {
                           nGroup->bitVector[k]=1;
                           nGroup->count++;
                           for(int z=l, y=0;z<=m+1;z++,y++)
                           {
                            //   if(nGroup->sum.size()==0)
                               if(nGroup->count==1)
                                nGroup->sum.push_back(timeSeries[k][z]);
                               else
                                nGroup->sum[y]=nGroup->sum[y]+timeSeries[k][z];
                           }

                        }
                        else
                            nGroup->bitVector[k]=0;

                    }
                   if(nGroup->count>=1)
                   {
                       nGroup->centroidLen=groupArray[gri].centroidLen+1;
                        for(int z=0;z<nGroup->centroidLen;z++)
                        {
                            nGroup->centroid.push_back(nGroup->sum[z]/nGroup->count);
                        }


                       nGroup->id=groupCounter++;
                      //push this group in groups vector
                       groupArray.push_back(*nGroup);
                       //push this group id in time interval
                       Time[l][m+1].groupid.push_back(nGroup->id);
                       Time[l][m+1].count++;

                   }
                   delete nGroup;
                }


             }
         //for this time interval for each group calculate the sum of ED distances with other centroids
         double DcEDDist=0;
         double tempDcDist=0;
         for(int i=0;i<Time[l][m+1].count;i++)
         {
             //pick first group
             int groupId=Time[l][m+1].groupid[i];
             //calculate the distance of its centroid with all other
             for(int j=0;j<Time[l][m+1].count;j++)
             {

                 if(j!=i)
                 {
                     int groupId2=Time[l][m+1].groupid[j];
                     for(int k=0;k<groupArray[groupId].centroid.size();k++)
                     {
                         tempDcDist+=pow((groupArray[groupId].centroid[k]-groupArray[groupId2].centroid[k]),2.0);

                     }
                     DcEDDist+=sqrt(tempDcDist);
                     tempDcDist=0;
                 }
             }
             //save the DCEDDist in the structure
             Time[l][m+1].DcDistance.push_back(DcEDDist);
             DcEDDist=0;
         }
   }
}
