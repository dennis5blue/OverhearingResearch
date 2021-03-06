#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
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
std::vector<int> CrossEntropy(const int, const int, const double, const double, const matrix<double>, const std::vector<double>);
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
    
    //---------------------------------------------------------------------------------//
    //-----Create Covariance matrix, Residual energy vector and set Fidelity ratio-----//
    //---------------------------------------------------------------------------------//
    matrix<double> covMatrix(numMachines, numMachines);
    covMatrix = Covm(myTopology, 1, corrLevel);
    //cout << covMatrix << endl;
    const double totalEntropy[15] = {402.621123, 387.475685, 372.827010, 358.572206, 344.720873, 331.319605, 306.013552, 282.762930, 261.459375, 241.906822, 223.895835, 207.232974, 191.749306, 177.300603, 163.764679};
    const double fidelityRatio = 0.7*totalEntropy[corrIndex];
    std::vector<double> vecResidualEnergy;
    for (int i=0; i!=numMachines; i++) vecResidualEnergy.push_back(1);
    int countLifetime = 0;
    double gatheredEntropy = fidelityRatio;
    
    //---------------------------------------------------------//
    //------------Start the corss entropy algorithm------------//
    //---------------------------------------------------------//
    do {
      std::vector<int> vecMachineSelection;
      const int numIter = 50;
      const int baseSize = 8*numIter;
      const double alpha = 0.7;
      const double beta = 0;
      vecMachineSelection = CrossEntropy(scheAlgorithm, baseSize, alpha, beta, myTopology, vecResidualEnergy);
    } while (gatheredEntropy > fidelityRatio);
    
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

std::vector<int> SINRscheduling(const std::vector<int> m_vecMachineSelection, const std::vector<double> m_ChannelCapacity)
{
    int m_numSelected = m_vecMachineSelection.size();
    std::vector<int> m_vecMachineSchedule;
    std::vector< pair<double,int> > m_vecChannelAndMachine;
    for (int i=0; i!=m_numSelected; i++) {
        int m_number = m_vecMachineSelection[i];
        double m_channel = m_ChannelCapacity[i];
        m_vecChannelAndMachine.push_back(make_pair(m_channel,m_number));
    }
    std::sort(m_vecChannelAndMachine.begin(), m_vecChannelAndMachine.end());
    for (int j=0; j!=m_numSelected; j++) {
        int m_scheduledMachine = m_vecChannelAndMachine[j].second;
        m_vecMachineSchedule.push_back(m_scheduledMachine);
    }
    return m_vecMachineSchedule;
}

std::vector<int> CrossEntropy(const int m_scheAlgorithm, const int m_baseSize, const double m_alpha, const double m_beta, const matrix<double> m_myTopology, const std::vector<double> m_vecResidualEnergy)
{
  std::vector<int>                m_vecMachineSelection;
  std::vector<double>             m_vecCrossEntropyProb;
  std::vector< pair<double,int> > m_vecEnergyAndMachine;
  std::vector<int>                m_vecNoPowerMachine;
  for (int e=0; e!=numMachines; e++){ //record machines which still have energy
    if(m_vecResidualEnergy[e] > 0)
      m_vecEnergyAndMachine.push_back(make_pair(m_vecResidualEnergy[e],e));
    else  m_vecNoPowerMachine.push_back(e);
  }
  double m_maxEnergy = *max_element(m_vecResidualEnergy.begin(), m_vecResidualEnergy.end());
  double m_minEnergy = *min_element(m_vecResidualEnergy.begin(), m_vecResidualEnergy.end());
  
  //Initial probability
  for(int i=0; i!=numMachines; i++){
    if(m_maxEnergy - m_minEnergy > 0){
      double m_deltaTemp = (m_vecResidualEnergy[i]-m_minEnergy)/(m_maxEnergy-m_minEnergy);
      double m_pbTemp = 0.5 + 0.5*(m_deltaTemp-0.5);
      m_vecCrossEntropyProb.push_back(m_pbTemp);
    } else
      m_vecCrossEntropyProb.push_back(0.5);
  }
  
  //--------------------------------------------------//
  //-------for loop for cross entropy algorithm-------//
  //--------------------------------------------------//
  for(int m_crossEntropyCounter=0; m_crossEntropyCounter!=50; m_crossEntropyCounter++){
    
    //Calculate machine selection by bernoulli trail
    std::vector< pair< double, std::vector<bool> > > m_vecTimeAndSelection;
    std::vector<bool> m_vecBernoulliTemp;
    double m_randomSeed = 0;
    bool m_selectSeed = 0;
    for(int j=0; j!=m_baseSize; j++){
      for(int k=0; k!=numMachines; k++){
        m_randomSeed = (double)rand()/(double)RAND_MAX;
        if(m_vecCrossEntropyProb[k] - m_randomSeed > 0) m_selectSeed = 1;
        else  m_selectSeed = 0;
        if(j==0)  m_vecBernoulliTemp.push_back(m_selectSeed);
        else  m_vecBernoulliTemp[k] = m_selectSeed;
      }
      for(int c=0; c!=m_vecNoPowerMachine.size(); c++) //prevent form selecting no power machine
        m_vecBernoulliTemp[m_vecNoPowerMachine[c]] = 0;
      if(j == 0){
        for(int s=0; s!=numMachines; s++)
          if(m_vecBernoulliTemp[s] == 1)  m_vecMachineSelection.push_back(s);
      }
      else{
        m_vecMachineSelection.clear();
        for(int s2=0; s2!=numMachines; s2++)
          if(m_vecBernoulliTemp[s2] == 1)  m_vecMachineSelection.push_back(s2);
      }
      //for(int kk=0; kk!=m_vecMachineSelection.size(); kk++)
        //cout << m_vecMachineSelection[kk] << " ";
      
      //Add machines with highest residual energy for fitting the fidelity ratio
      //TO DO
      
      std::vector<double> m_vecChannelCapacity;
      for(int q=0; q!=numMachines; q++)
        m_vecChannelCapacity.push_back(m_myTopology(q,2));
      
      //SINR scheduling and calculate total time resource needed
      std::vector<int> m_vecMachineSchedule;
      int m_numSelected = m_vecMachineSelection.size();
      if(m_scheAlgorithm == 1){
        m_vecMachineSchedule = SINRscheduling(m_vecMachineSelection, m_vecChannelCapacity);
        cout << "Number of selected machines = " << m_numSelected << endl;
        cout << "Machine schedule = ";
        for(int printsche=0; printsche!=m_numSelected; printsche++)
          cout << m_vecMachineSchedule[printsche] << " ";
        cout << endl;
      }//end SINR scheduling
      
      m_vecTimeAndSelection.push_back(make_pair(0,m_vecBernoulliTemp));
    } //end generating 8*N vector for cross entropy to run
      
  } //end cross entropy

  return m_vecMachineSelection;
}
