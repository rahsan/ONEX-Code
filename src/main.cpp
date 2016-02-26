/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected � 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "Groups.h"
//#include "online.h"
#include "trillionDTW.h"
#include "TimeSeries.h"
#include "OnlineSession.h"

#include "trillionDTW.h"

#include <iostream>

#include <sys/time.h>
#include <math.h>
#include <time.h>

using namespace std;

int main()
{
    cout << "Welcome to ONEX." << endl;

    OnlineSession *session = new OnlineSession();
    session->run(&cin, &cout);

    //double d;
    //double t1, t2;  //time variables

    //clock_t start = clock();
    //calculate("data/ShapesAllS.txt", "data/Query.txt", 128);
    //cout << "Trillion " << (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;


    //GroupOperation *groupObj = new GroupOperation(1.5);
    //groupObj->readFile("data/ItalyPower.txt", 67, 24);
    //groupObj->groupOp();
    //groupObj->printGroups();
    
    //groupObj->printTS();
/*
    GroupOperation *groupObj=new GroupOperation(1.5);     //ST passed as input

    groupObj->readFile("data/ItalyPower.txt", 67,24);
    //groupObj->readFile("data/ShapesAllS.txt", 600,513);
 //   groupObj->readFile("data/Synthetic2.txt", 400,100);
    //groupObj->readFile("data/FaceAll.txt", 560,132);
  //  groupObj->printTS();

    //groupObj->readFile("data/ItalyPower.txt", 67,24);
    clock_t start = clock();
    groupObj->groupOp();
    cout << "Offline Construction " << (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;


    //gettimeofday(&tim, NULL);
   // t1=tim.tv_sec+(tim.tv_usec/1000000.0);

 //   generateData("data/Synthetic.txt");

*/

    // ST, number of time series in original file, length of time series
    //onlineOperation *onlineObj=new onlineOperation(3,67,24);
    //onlineOperation *onlineObj=new onlineOperation(1.5,560,132);
    //onlineOperation *onlineObj=new onlineOperation(1,67,24);
    //onlineOperation *onlineObj=new onlineOperation(3,3,24);

    //onlineObj->readTimeSeries("data/ItalyPower.txt",2);
    //onlineObj->readTimeSeries("data/FaceAll.txt",true);
    //onlineObj->readTimeSeries("data/ItalyPower.txt",1);

    //read time, then groups generated in offline step
    //onlineObj->readTime();
    //onlineObj->readGroups();


    //read time series data from file
    //onlineObj->readQueryFile(3,"data/query.txt");
    //onlineObj->readQueryFile(4,"data/query.txt",true);


    //    clock_t start1 = clock();
    //onlineObj->kSimilar(true,1,false,0,2);
    //cout << "ONEX" << (clock()-start1)/double(CLOCKS_PER_SEC) << "s]" << endl;

    /*
    cout<<"Spring Calculation"<<endl;
    onlineObj->springCalculation(1,1);
    cout << "Spring " << (clock()-start1)/double(CLOCKS_PER_SEC) << "s]" << endl;
    */


    /*
    onlineObj->dominantOutlier(2,20);
    cout << "ONEX "<< (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;
    onlineObj->springOutlier(0,90);
    cout << "Spring " << (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;
    */
/*
    vector<double> temp1;
    vector<double> temp2;
    temp1.push_back(4);
    temp1.push_back(5);
    temp2.push_back(4);
    temp2.push_back(5);
    temp2.push_back(2);
    cout<<onlineObj->simpleDTW(temp1,temp2);
    //Close the file


    if( fp == NULL )
    {
        cout<<"File not found";
        exit(2);
    }
    int i=0;
    int j=0;
    int temp=fscanf(fp,"%lf",&d);
    while(temp != EOF)
    {
        cout<<"temp "<<temp<<" ";
        if(temp == EOL)
        {
            i++;
            j=0;
        }

        //cout<<d<<" ";
        timeSeries[i][j]=d;
    //    cout<<timeSeries[i][j]<<" ";
        j++;
        temp=fscanf(fp,"%lf",&d);
    }
    cout<<endl;
    cout<<"i "<<i<<" j "<<j;

    double *Q;             // query array
    double *T;             // array of current data
    int *order;            // ordering of query by |z(q_i)|
    double bsf;            // best-so-far
    int m;                 // length of query
    long long loc = 0;     // answer: location of the best-so-far match

    double d;
    long long i , j ;
    double ex , ex2 , mean, std;

    double t1,t2;

    t1 = clock();

    bsf = INF;
    i = 0;
    j = 0;
    ex = ex2 = 0;

  //  if (argc<=3)      error(4);

    //fp = fopen(argv[1],"r");
    //input file
    fp = fopen("data/input.txt","r");
    if( fp == NULL )
    {
        cout<<"File not found";
        exit(2);
    }

    //query file
    //qp = fopen(argv[2],"r");
    qp = fopen("data/temp.txt","r");
    if( qp == NULL )
        exit(2);

   // m = atol(argv[3]);
    m =2;
    /// Array for keeping the query data
    Q = (double *)malloc(sizeof(double)*m);
    if( Q == NULL )
        error(1);

    /// Read the query data from input file and calculate its statistic such as mean, std
    while(fscanf(qp,"%lf",&d) != EOF && i < m)
    {
        ex += d;
        ex2 += d*d;
        Q[i] = d;
        i++;
    }
    mean = ex/m;
    cout<<"Mean "<<mean;
    std = ex2/m;

    std = sqrt(std-mean*mean);
    cout<<"Std "<<std;
    fclose(qp);

    /// Do z_normalixation on query data
    for( i = 0 ; i < m ; i++ )
         {
        Q[i] = (Q[i] - mean)/std;
        cout<<"Q[i] "<<Q[i]<<" ";
    }
    cout<<endl;

    /// Sort the query data
    order = (int *)malloc(sizeof(int)*m);
    if( order == NULL )
        error(1);
    Index *Q_tmp = (Index *)malloc(sizeof(Index)*m);
    if( Q_tmp == NULL )
        error(1);
    for( i = 0 ; i < m ; i++ )
    {
        Q_tmp[i].value = Q[i];
        Q_tmp[i].index = i;
    }
    cout<<endl;
    qsort(Q_tmp, m, sizeof(Index),comp);
    for( i=0; i<m; i++)
    {   Q[i] = Q_tmp[i].value;
        cout<<"After Sorting Q[i] "<<Q[i]<<" ";
        order[i] = Q_tmp[i].index;
        cout<<"order i "<<order[i]<<" ";
    }
    cout<<endl; //blank lines
    free(Q_tmp);



    /// Array for keeping the current data; Twice the size for removing modulo (circulation) in distance calculation
    T = (double *)malloc(sizeof(double)*2*m);
    if( T == NULL )
        error(1);

    double dist = 0;
    i = 0;
    j = 0;
    ex = ex2 = 0;

    /// Read data file, one value at a time
    while(fscanf(fp,"%lf",&d) != EOF )
    {
        ex += d;
        ex2 += d*d;
        T[i%m] = d;
        T[(i%m)+m] = d;

        /// If there is enough data in T, the ED distance can be calculated
        if( i >= m-1 )
        {
            /// the current starting location of T
            j = (i+1)%m;

            /// Z_norm(T[i]) will be calculated on the fly
            mean = ex/m;
            std = ex2/m;
            std = sqrt(std-mean*mean);

            /// Calculate ED distance
            dist = distance(Q,T,j,m,mean,std,order,bsf);
            if( dist < bsf )
            {
                bsf = dist;
                loc = i-m+1;
            }
            ex -= T[j];
            ex2 -= T[j]*T[j];
        }
        i++;
    }
    fclose(fp);
    t2 = clock();

    cout << "Location : " << loc << endl;
    cout << "Distance : " << sqrt(bsf) << endl;
    cout << "Data Scanned : " << i << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
*/
    return 0;
}

int main_bac()
{

//    double d;
//    double t1, t2;  //time variables
/*
    clock_t start = clock();
    calculate();
    cout << "Trillion " << (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;
*/

/*
    GroupOperation *groupObj=new GroupOperation(1.5);     //ST passed as input

    groupObj->readFile("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/ItalyPower.txt", 67,24);
    //groupObj->readFile("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/ShapesAllS.txt", 600,513);
 //   groupObj->readFile("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/Synthetic2.txt", 400,100);
    //groupObj->readFile("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/FaceAll.txt", 560,132);
  //  groupObj->printTS();

    //groupObj->readFile("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/ItalyPower.txt", 67,24);
    clock_t start = clock();
    groupObj->groupOp();
    cout << "Offline Construction " << (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;


    //gettimeofday(&tim, NULL);
   // t1=tim.tv_sec+(tim.tv_usec/1000000.0);

 //   generateData("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/Synthetic.txt");

*/

    //onlineOperation *onlineObj=new onlineOperation(3,67,24);// ST, number of time series in original file, length of time series
 //  onlineOperation *onlineObj=new onlineOperation(1.5,560,132);// ST, number of time series in original file, length of time series
    //onlineOperation *onlineObj=new onlineOperation(1,67,24);// ST, number of time series in original file, length of time series

    //onlineOperation *onlineObj=new onlineOperation(3,3,24);// ST, number of time series in original file, length of time series
    //onlineObj->readTime();
    //read groups genereted in offline step
    //onlineObj->readGroups();
    //onlineObj->readQueryFile(3,"C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/query.txt");



    //onlineObj->readQueryFile(4,"C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/query.txt",true);
    //read time series data from file

    //onlineObj->readTimeSeries("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/ItalyPower.txt",2);
    //onlineObj->readTimeSeries("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/FaceAll.txt",true);
    //onlineObj->readTimeSeries("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/ItalyPower.txt",1);


    //clock_t start1 = clock();
    //onlineObj->kSimilar(true,1,false,0,2);
    //cout << "ONEX" << (clock()-start1)/double(CLOCKS_PER_SEC) << "s]" << endl;
/*
   cout<<"Spring Calculation"<<endl;
    onlineObj->springCalculation(1,1);
    cout << "Spring " << (clock()-start1)/double(CLOCKS_PER_SEC) << "s]" << endl;
*/


 //   onlineObj->dominantOutlier(2,20);
 //   cout << "ONEX "<< (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;
 //   onlineObj->springOutlier(0,90);
 //   cout << "Spring " << (clock()-start)/double(CLOCKS_PER_SEC) << "s]" << endl;

/*
    vector<double> temp1;
    vector<double> temp2;
    temp1.push_back(4);
    temp1.push_back(5);
    temp2.push_back(4);
    temp2.push_back(5);
    temp2.push_back(2);
    cout<<onlineObj->simpleDTW(temp1,temp2);
    //Close the file


    if( fp == NULL )
    {
        cout<<"File not found";
        exit(2);
    }
    int i=0;
    int j=0;
    int temp=fscanf(fp,"%lf",&d);
    while(temp != EOF)
    {
        cout<<"temp "<<temp<<" ";
        if(temp == EOL)
        {
            i++;
            j=0;
        }

        //cout<<d<<" ";
        timeSeries[i][j]=d;
    //    cout<<timeSeries[i][j]<<" ";
        j++;
        temp=fscanf(fp,"%lf",&d);
    }
    cout<<endl;
    cout<<"i "<<i<<" j "<<j;

    double *Q;             // query array
    double *T;             // array of current data
    int *order;            // ordering of query by |z(q_i)|
    double bsf;            // best-so-far
    int m;                 // length of query
    long long loc = 0;     // answer: location of the best-so-far match

    double d;
    long long i , j ;
    double ex , ex2 , mean, std;

    double t1,t2;

    t1 = clock();

    bsf = INF;
    i = 0;
    j = 0;
    ex = ex2 = 0;

  //  if (argc<=3)      error(4);

    //fp = fopen(argv[1],"r");
    //input file
    fp = fopen("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/input.txt","r");
    if( fp == NULL )
    {
        cout<<"File not found";
        exit(2);
    }

    //query file
    //qp = fopen(argv[2],"r");
    qp = fopen("C:/Qt/Tools/QtCreator/bin/build-Onex-Desktop_Qt_5_1_1_MinGW_32bit-Debug/debug/temp.txt","r");
    if( qp == NULL )
        exit(2);

   // m = atol(argv[3]);
    m =2;
    /// Array for keeping the query data
    Q = (double *)malloc(sizeof(double)*m);
    if( Q == NULL )
        error(1);

    /// Read the query data from input file and calculate its statistic such as mean, std
    while(fscanf(qp,"%lf",&d) != EOF && i < m)
    {
        ex += d;
        ex2 += d*d;
        Q[i] = d;
        i++;
    }
    mean = ex/m;
    cout<<"Mean "<<mean;
    std = ex2/m;

    std = sqrt(std-mean*mean);
    cout<<"Std "<<std;
    fclose(qp);

    /// Do z_normalixation on query data
    for( i = 0 ; i < m ; i++ )
         {
        Q[i] = (Q[i] - mean)/std;
        cout<<"Q[i] "<<Q[i]<<" ";
    }
    cout<<endl;

    /// Sort the query data
    order = (int *)malloc(sizeof(int)*m);
    if( order == NULL )
        error(1);
    Index *Q_tmp = (Index *)malloc(sizeof(Index)*m);
    if( Q_tmp == NULL )
        error(1);
    for( i = 0 ; i < m ; i++ )
    {
        Q_tmp[i].value = Q[i];
        Q_tmp[i].index = i;
    }
    cout<<endl;
    qsort(Q_tmp, m, sizeof(Index),comp);
    for( i=0; i<m; i++)
    {   Q[i] = Q_tmp[i].value;
        cout<<"After Sorting Q[i] "<<Q[i]<<" ";
        order[i] = Q_tmp[i].index;
        cout<<"order i "<<order[i]<<" ";
    }
    cout<<endl; //blank lines
    free(Q_tmp);



    /// Array for keeping the current data; Twice the size for removing modulo (circulation) in distance calculation
    T = (double *)malloc(sizeof(double)*2*m);
    if( T == NULL )
        error(1);

    double dist = 0;
    i = 0;
    j = 0;
    ex = ex2 = 0;

    /// Read data file, one value at a time
    while(fscanf(fp,"%lf",&d) != EOF )
    {
        ex += d;
        ex2 += d*d;
        T[i%m] = d;
        T[(i%m)+m] = d;

        /// If there is enough data in T, the ED distance can be calculated
        if( i >= m-1 )
        {
            /// the current starting location of T
            j = (i+1)%m;

            /// Z_norm(T[i]) will be calculated on the fly
            mean = ex/m;
            std = ex2/m;
            std = sqrt(std-mean*mean);

            /// Calculate ED distance
            dist = distance(Q,T,j,m,mean,std,order,bsf);
            if( dist < bsf )
            {
                bsf = dist;
                loc = i-m+1;
            }
            ex -= T[j];
            ex2 -= T[j]*T[j];
        }
        i++;
    }
    fclose(fp);
    t2 = clock();

    cout << "Location : " << loc << endl;
    cout << "Distance : " << sqrt(bsf) << endl;
    cout << "Data Scanned : " << i << endl;
    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;
*/
   return 0;
}
