#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <random>
//#include <boost/program_options.hpp>
#include <utility>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "math.h"

//namespace po = boost::program_options;

using namespace boost::numeric::ublas;
using namespace std;

matrix<double> Covm(const matrix<double>, const double, const double);
std::vector<int> SINRscheduling(const std::vector<int>, const matrix<double>);
const int numMachines = 50;

int main( int argc, char *argv[])
{
    const double ax[15] = {500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000};
    const int corrIndex = atoi(argv[1])-1;
    const int scheAlgorithm = atoi(argv[2]);
    const int numOverSlot = atoi(argv[3]);
    const double corrLevel = 1/(ax[corrIndex]);
    cout << "Correlation level = " << corrLevel << endl;
    cout << "Scheduling Alogrithm Index = " << scheAlgorithm << endl;
    cout << "Number of Overhearing Time Slots = " << numOverSlot << endl;
    //-----------------------------------------------------//
    //load the target topology, stored in matrix myTopology//
    //              topology format: x y C z               //
    //-----------------------------------------------------//
    matrix<double> myTopology(numMachines, 4);
    fstream topologyFile;
    topologyFile.open("topology_x.txt");
    string strLine;
    int countToken = 0;
    while (getline(topologyFile, strLine)) {
        stringstream ss(strLine);
        double x;
        ss >> x;
        myTopology(countToken,0) = x;
        countToken++;
    }
    topologyFile.close();
    topologyFile.open("topology_y.txt");
    countToken = 0;
    while (getline(topologyFile, strLine)) {
        stringstream ss(strLine);
        double y;
        ss >> y;
        myTopology(countToken,1) = y;
        countToken++;
    }
    topologyFile.close();
    topologyFile.open("topology_C.txt");
    countToken = 0;
    while (getline(topologyFile, strLine)) {
        stringstream ss(strLine);
        double C;
        ss >> C;
        myTopology(countToken,2) = C;
        countToken++;
    }
    topologyFile.close();
    topologyFile.open("topology_z.txt");
    countToken = 0;
    while (getline(topologyFile, strLine)) {
        stringstream ss(strLine);
        double z;
        ss >> z;
        myTopology(countToken,3) = z;
        countToken++;
    }
    topologyFile.close();
    //cout << myTopology << endl;
    
    //----------------------------------------------------------------------
    //Create covariance matrix and initialize settings
    matrix<double> covMatrix(numMachines, numMachines);
    covMatrix = Covm(myTopology, 1, corrLevel);
    //cout << covMatrix << endl;
    
    const double totalEntropy[15] = {402.621123, 387.475685, 372.827010, 358.572206, 344.720873, 331.319605, 306.013552, 282.762930, 261.459375, 241.906822, 223.895835, 207.232974, 191.749306, 177.300603, 163.764679};
    const double fidelityRatio = 0.7*totalEntropy[corrIndex];
    matrix<double> residualEnergy(numMachines,1);
    for (int i=0; i!=numMachines; i++) residualEnergy(i,0) = 1;
    int countLifetime = 0;
    
    //---------------------------------------------------------//
    //------------Start the corss entropy algorithm------------//
    //---------------------------------------------------------//
    std::vector<int> vecMachineSelection;
    const int numIter = 50;
    const int baseSize = 8*numIter;
    const double alpha = 0.7;
    std::vector<double> vecCrossEntropyPb;
    std::bernoulli_distribution distribution(0.5);
    
    //for test
    cout << "Machine Selection = ";
    for (int i=0; i!=10; i++){
        vecMachineSelection.push_back(i);
        cout << vecMachineSelection[i] << " ";
    }
    cout << endl;
    //end for test
    int numSelected = vecMachineSelection.size();
    std::vector<int> vecMachineSchedule;
    if (scheAlgorithm == 1) {
        vecMachineSchedule = SINRscheduling(vecMachineSelection, myTopology);
        // cout for test
        cout << "Number of selected machines = " << numSelected << endl;
        cout << "Machine Schedule = " ;
        for (int i=0; i!=numSelected; i++) cout << vecMachineSchedule[i] << " ";
        cout << endl;
        //end cout
    }
    
    return 0;
}

matrix<double> Covm(const matrix<double> m_myTopology, const double m_sigma, const double m_corrLevel)
{
    matrix<double> m_covMatrix(numMachines,numMachines);
    for (int i=0; i!=numMachines; i++) {
        for (int j=0; j!=numMachines; j++) {
            m_covMatrix(i,j) = pow(m_sigma,2)*exp(-m_corrLevel*(pow((m_myTopology(i,0)-m_myTopology(j,0)),2)+pow((m_myTopology(i,1)-m_myTopology(j,1)),2)));
        }
    }
    return m_covMatrix;
}

std::vector<int> SINRscheduling(const std::vector<int> m_vecMachineSelection, const matrix<double> m_myTopology)
{
    int m_numSelected = m_vecMachineSelection.size();
    std::vector<int> m_vecMachineSehedule;
    std::vector< pair<double,int> > m_vecChannelAndMachine;
    for (int i=0; i!=m_numSelected; i++) {
        int m_number = m_vecMachineSelection[i];
        double m_channel = m_myTopology(i,2);
        m_vecChannelAndMachine.push_back(make_pair(m_channel,m_number));
    }
    std::sort(m_vecChannelAndMachine.begin(), m_vecChannelAndMachine.end());
    for (int j=0; j!=m_numSelected; j++) {
        int m_scheduledMachine = m_vecChannelAndMachine[j].second;
        m_vecMachineSehedule.push_back(m_scheduledMachine);
    }
    return m_vecMachineSehedule;
}
