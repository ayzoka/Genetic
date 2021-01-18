#include "GA.h"

namespace QKD
{
  // t = which round are we computing
  void EvolveQKD::constructEMOps(CS& cs, MessageSettings& mset, algebra::mat& rho00, algebra::mat& rho01, algebra::mat& rho10, algebra::mat& rho11, std::vector <algebra::mat>& output, int t)
  {
    if(output.size() != 4*mset.totalDim(t))
      output.resize(4*mset.totalDim(t));

    std::vector <int> m;
    m.resize(t);
    for(int i=0; i<t; ++i)
      m[i] = 0;

    for(int ka=0; ka<=1; ++ka){
      for(int kb=0; kb<=1; ++kb){
	do{
	  // i=0, j=0, etc;
	  double pA00 = cs.messages[t-1].prA[getMessageProbIndex(mset, t, m, ka, 0)];
	  double pA01 = pA00;
	  double pA10 = cs.messages[t-1].praA[getMessageProbIndex(mset, t, m, ka, 1)];
	  double pA11 = pA10;

	  double pB00 = cs.messages[t-1].prB[getMessageProbIndex(mset, t, m, kb, 0)];
	  double pB10 = pB00;
	  double pB01 = cs.messages[t-1].prB[getMessageProbIndex(mset, t, m, kb, 1)];
	  double pB11 = pB10;

	  if(t == 1){
	    output[getMessageOpIndex(ka, kb, m, mset)] = rho00*(pA00*pB00)
	      + rho01*(pA01*pB01)
	      + rho10*(pA10*pB10)
	      + rho11*(pA11*pB11);
	  }

	}while(getNextMessage(m, mset));
      }
    }

    }

  void EvolveQKD::calculateFitness(quantum::SimHelper& sim, CS& cs, QKD::Protocol& prot, QKD::Attack& atk, MessageSettings& mset)
  {
    // (rho0 = rho00 + rho01; rho1 = rho10+rho11; rhoE = rho0+rho1)
    algebra::mat rho00, rho01, rho10, rho11, rho0, rho1, rhoE;
    std::vector <algebra::mat> rhoEvec;
    std::vector <double> ev0, ev1, evE;
    // compute key rate:

    if(mset.use_messages)
      rhoeEvec.resize(4*mset.totalDim());

    atk.begin();

    double best = 200;
    while(!atk.getNext(Evec[0], Evec[1], Evec[2], Evec[3])){
      prot.constructEOps(sim, cs, Evec, rho00, rho01, rho10, rho11);

      if(mset.use_messages){
	// construct list of ops for E's system:
	constructEMOps(cs, mset, rho00, rho01, rho10, rho11, rhoEVec, mset.numRounds());
      }
      else{
	// key rate when no messages are sent:
	rho00.add(&rho01, &rho0);
	rho11.add(&rho10, &rho1);
	
	rho0.add(&rho1, &rhoE);
	
	rho0.eigenValues(ev0);
	rho1.eigenValues(ev1);
	rhoE.eigenValues(evE);
	
	double SAE = 0, SE = 0;
	for(int i=0; i<4; ++i){
	  SAE = SAE - (ev0[i]*SafeLog(ev0[i])) - (ev1[i]*SafeLog(ev1[i]));
	  SE = SE - (evE[i]*SafeLog(evE[i]));
	}
	
	double key00 = rho00.traceRe();
	double key01 = rho01.traceRe();
	double key10 = rho10.traceRe();
	double key11 = rho11.traceRe();
	
	double rate = SAE - SE - HAB(key00, key01, key10, key11);
	if(rate < best)
	  best = rate;
      }
    }

    cs.fitness = best;
  }

  void EvolveQKD::evolve(QKD::Protocol& prot, QKD::Attack& atk)
  {
    MessageSettings2 mset;
    const int TOUR = 5;
    const int POPULATION_SIZE = 100;
    const int MAX_ITERATIONS = 100;

    Evec.resize(4);
    for(int i=0; i<4; ++i)
      Evec[i].create(4,1);


    quantum::SimHelper sim;

    std::cout << "Setup Prot...";
    prot.setup(sim);
    std::cout << "done.\n";

    std::list <CS> population, tempPopulation;
    // setup init population:
    std::cout << "Setup Init Population";
    for(int i=0; i<POPULATION_SIZE; ++i){
      CS cs;
      prot.setupCS(cs);
      cs.randomize();
      calculateFitness(sim, cs, prot, atk);
      population.push_back(cs);
      std::cout << ".";
    }

    std::cout << "done\n";
    population.sort();

    std::list <CS>::iterator Iter, iterP1, iterP2;

    for(int iter=0; iter<MAX_ITERATIONS; ++iter){

      tempPopulation.clear();

      
      int count = 0;
      for(Iter = population.begin(); Iter != population.end(); ++Iter){
	if(count <= 1)
	  tempPopulation.push_back(*Iter);

	// tournament:
	int indexP1 = POPULATION_SIZE+1;
	for(int k=0; k<TOUR; ++k){
	  int r = rand()%POPULATION_SIZE;
	  if(r < indexP1)
	    indexP1 = r;
	}
	int indexP2 = POPULATION_SIZE+1;
	for(int k=0; k<TOUR; ++k){
	  int r = rand()%POPULATION_SIZE;
	  if(r < indexP2)
	    indexP2 = r;
	}

	iterP1 = population.begin();
	iterP2 = population.begin();
	for(int k=0; k<indexP1; ++k)
	  ++iterP1;
	for(int k=0; k<indexP2; ++k)
	  ++iterP2;

	CS child;
	child = iterP1->crossover(*iterP2);
	if(rand()%100 < 20)
	  child.mutate();

	calculateFitness(sim, child, prot, atk);
	tempPopulation.push_back(child);

	++count;
      }

      population = tempPopulation;
      population.sort();

      prot.printInfo(sim, *population.begin(), atk);
    }
  }
}
