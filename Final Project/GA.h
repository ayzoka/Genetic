#pragma once

#include "CS.h"
#include "lib-quantum/SimHelper.h"
#include "attack.h"
#include "lib-quantum/Protocol.h"
#include <vector>

namespace QKD
{

struct GAOutput
{
  double rate; // rate of best solution
  CS cs;
};
struct GAInput
{
  int mutation_rate;
  int crossover_type; // 1 = blend style, 2 = pt. wise

  int tournament_size;
  int population_size;
  int max_iterations;
};

class EvolveQKD
{
  double randFloat(double small, double large)
  {
    int r = rand()%RAND_MAX;
    double r2 = r / (double)RAND_MAX;

    return (large-small)*r2 + small;
  }

  // computes H(A|B) assuming kij means A = i, B = j
  double HAB(double key00, double key01, double key10, double key11)
  {
    double first = -key00*SafeLog(key00) - key01*SafeLog(key01) - key10*SafeLog(key10) - key11*SafeLog(key11);

    if(key00+key10 > 1.00001)
      std::cout << "HAB ERROR!\n";
    double second = entropy(key00+key10);

    return first-second;
  }
  // computes H(B|A) assuming kij means A = i, B = j
  double HBA(double key00, double key01, double key10, double key11)
  {
    double first = -key00*SafeLog(key00) - key01*SafeLog(key01) - key10*SafeLog(key10) - key11*SafeLog(key11);

    if(key00+key01 > 1.00001)
      std::cout << "HAB ERROR!\n";
    double second = entropy(key00+key01);

    return first-second;
  }


  double SafeLog(double x)
  {
    if(x < 0.0000001)
      return 0;
    return log(x)/log(2.0);
  }
  double entropy(double x)
  {
    if(x > 1.00001){
      std::cout << "Entropy error: x = " << x << "\n";
      return 0;
    }

    return -x*SafeLog(x) - (1-x)*SafeLog(1-x);
  }


  std::vector <algebra::mat> Evec;

 public:


  void constructEMOps(CS& cs, MessageSettings& mset, algebra::mat& rho00, algebra::mat& rho01, algebra::mat& rho10, algebra::mat& rho11, std::vector <algebra::mat>& output, int t)
  {

    std::vector <algebra::mat> outputPrev;

    if(t > 1){
      constructEMOps(cs, mset, rho00, rho01, rho10, rho11, outputPrev, t-1);
    }

    if(output.size() != 4*mset.totalDim(t))
      output.resize(4*mset.totalDim(t));

    std::vector <int> m;
    m.resize(t);
    for(int i=0; i<t; ++i)
      m[i] = 0;

    double tr = 0;
    for(int ka=0; ka<=1; ++ka){
      for(int kb=0; kb<=1; ++kb){
	do{
	  // i=0, j=0, etc;
	  double pA00 = cs.messages[t-1].prA[mset.getMessageProbIndex(t, m, ka, 0)];
	  double pA01 = pA00;
	  double pA10 = cs.messages[t-1].prA[mset.getMessageProbIndex(t, m, ka, 1)];
	  double pA11 = pA10;

	  double pB00 = cs.messages[t-1].prB[mset.getMessageProbIndex(t, m, kb, 0)];
	  double pB10 = pB00;
	  double pB01 = cs.messages[t-1].prB[mset.getMessageProbIndex(t, m, kb, 1)];
	  double pB11 = pB01;
	  if(t == 1){

	    algebra::mat::addMulti(rho00, pA00*pB00, rho01, pA01*pB01,
				   rho10, pA10*pB10, rho11, pA11*pB11,
				   output[getMessageOpIndex(ka, kb, m, t, mset)]);

	    tr += output[getMessageOpIndex(ka, kb, m, t, mset)].traceRe();
	    //output[getMessageOpIndex(ka, kb, m, t, mset)] = rho00*(pA00*pB00);
	      //+ rho01*(pA01*pB01)
	      //+ rho10*(pA10*pB10)
	      //+ rho11*(pA11*pB11);
	  }
	  else{
	    algebra::mat::addMulti(
		outputPrev[getMessageOpIndex(0,0, m, t-1, mset)], pA00*pB00,
		outputPrev[getMessageOpIndex(0,1,m,t-1,mset)], pA01*pB01,
		outputPrev[getMessageOpIndex(1,0,m,t-1,mset)], pA10*pB10,
		outputPrev[getMessageOpIndex(1,1,m,t-1,mset)], pA11*pB11,
	        output[getMessageOpIndex(ka, kb, m, t, mset)]);

	    tr += output[getMessageOpIndex(ka, kb, m, t, mset)].traceRe();
	  }

	}while(mset.getNextMessage(m));

	
      }
    }
    //std::cout << "TR = " << tr << "\n";
    if(fabs(tr-1) > .000001){
      std::cout << "\t\tWARNING: Tr in construct message ops is " << tr << "\n";
    }

    }

  double partialEntropy(algebra::mat& M, std::vector <double>& ev)
  {
    M.eigenValues(ev);

    double sum=0;
    for(int i=0; i<ev.size(); ++i)
      sum = sum - ev[i]*SafeLog(ev[i]);

    return sum;
  }

  void calculateFitness(quantum::SimHelper& sim, CS& cs, QKD::Protocol* prot, QKD::Attack& atk, MessageSettings& mset)
  {
    bool directRecon = true;

    // (rho0 = rho00 + rho01; rho1 = rho10+rho11; rhoE = rho0+rho1)
    algebra::mat rho00, rho01, rho10, rho11, rho0, rho1, rhoE;
    std::vector <algebra::mat> rhoEvec;
    std::vector <double> ev0, ev1, evE;
    // compute key rate:

    //std::cout << "dim=" << 4*mset.totalDim()<<"\n";
    if(mset.use_messages)
      rhoEvec.resize(4*mset.totalDim());

    atk.begin();

    double best = 200;
    while(!atk.getNext(Evec[0], Evec[1], Evec[2], Evec[3])){
      prot->constructEOps(sim, cs, Evec, rho00, rho01, rho10, rho11);

      if(mset.use_messages){
	// construct list of ops for E's system:
	constructEMOps(cs, mset, rho00, rho01, rho10, rho11, rhoEvec, mset.numRounds());

	// now go through each message:
	std::vector <int> m;
	int t = mset.numRounds();
	m.resize(t);
	for(int i=0; i<t; ++i)
	  m[i] = 0;

	double SAE = 0, SE = 0;
	double key00=0, key01=0, key10=0, key11=0;
	do{
	  int indexKey00 = getMessageOpIndex(0, 0, m, t, mset);
	  int indexKey01 = getMessageOpIndex(0, 1, m, t, mset);
	  int indexKey10 = getMessageOpIndex(1, 0, m, t, mset);
	  int indexKey11 = getMessageOpIndex(1, 1, m, t, mset);

	  key00 += rhoEvec[indexKey00].traceRe();
	  key01 += rhoEvec[indexKey01].traceRe();
	  key10 += rhoEvec[indexKey10].traceRe();
	  key11 += rhoEvec[indexKey11].traceRe();

	  // Direct:
	  if(directRecon){
	    rhoEvec[indexKey00].add(&rhoEvec[indexKey01], &rho0);
	      rhoEvec[indexKey10].add(&rhoEvec[indexKey11], &rho1);
	      rho0.add(&rho1, &rhoE);
	  }
	  else{
	    // Reverse:
	    rhoEvec[indexKey00].add(&rhoEvec[indexKey10], &rho0);
	    rhoEvec[indexKey01].add(&rhoEvec[indexKey11], &rho1);
	    rho0.add(&rho1, &rhoE);
	  }

	  SAE += partialEntropy(rho0, ev0); // get e.v. of rho0 and add ent. fct
	  SAE += partialEntropy(rho1, ev1);

	  SE += partialEntropy(rhoE, evE);
	}while(mset.getNextMessage(m));

	if(key00 < 0 || key10 < 0 || key01 < 0 || key11 < 0 || fabs(key00+key01+key10+key11 - 1) > .000001){
	  std::cout << "WARNING: key sum = " << key00+key01+key10+key11 << "\n";
	  std::cout << "key00 = " << key00 << "\n";
	  std::cout << "key01 = " << key01 << "\n";
	  std::cout << "key10 = " << key10 << "\n";
	  std::cout << "key11 = " << key11 << "\n\n";
	  throw std::string("abc");
	}

	double rate = SAE - SE - HAB(key00, key01, key10, key11);
	if(!directRecon)
	  rate = SAE - SE - HBA(key00, key01, key10, key11);

	if(rate < best)
	  best = rate;
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
	/*
	  std::cout << "key00 = " << key00 << "\n";
	  std::cout << "key01 = " << key01 << "\n";
	  std::cout << "key10 = " << key10 << "\n";
	  std::cout << "key11 = " << key11 << "\n\n";
	*/
	if(key00 < 0 || key10 < 0 || key01 < 0 || key11 < 0 || fabs(key00+key01+key10+key11 - 1) > .000001){
	  std::cout << "WARNING: key sum = " << key00+key01+key10+key11 << "\n";
	  std::cout << "key00 = " << key00 << "\n";
	  std::cout << "key01 = " << key01 << "\n";
	  std::cout << "key10 = " << key10 << "\n";
	  std::cout << "key11 = " << key11 << "\n\n";
	  throw std::string("");
	}

	
	double rate = SAE - SE - HAB(key00, key01, key10, key11);
	if(rate < best)
	  best = rate;
      }
    }

    cs.fitness = best;
  }

  GAOutput evolve(QKD::Protocol* prot, QKD::Attack& atk, MessageSettings& mset, GAInput& input)
  {
    const int TOUR = input.tournament_size;
    const int POPULATION_SIZE = input.population_size;
    const int MAX_ITERATIONS = input.max_iterations;

    Evec.resize(4);
    for(int i=0; i<4; ++i)
      Evec[i].create(4,1);


    quantum::SimHelper sim;

    std::cout << "Setup Prot...";
    prot->setup(sim);
    std::cout << "done.\n";

    std::list <CS> population, tempPopulation;
    // setup init population:
    std::cout << "Setup Init Population";
    for(int i=0; i<POPULATION_SIZE; ++i){
      CS cs;
      prot->setupCS(cs);

      cs.randomize(mset);
      cs.bb84();
      calculateFitness(sim, cs, prot, atk, mset);
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
	if(count <= 1){
	  calculateFitness(sim, *Iter, prot, atk, mset); // part of test
	  tempPopulation.push_back(*Iter);
	}

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
	if(input.crossover_type == 1)
	  child = iterP1->crossover(*iterP2, mset);
	else
	  child = iterP1->crossover_ptwise(*iterP2, mset);

	if(rand()%100 < input.mutation_rate)
	  child.mutate(mset);
	//else
	//  child.mutateMessageOnly(mset);

	child.bb84();

	calculateFitness(sim, child, prot, atk, mset);
	tempPopulation.push_back(child);

	++count;
      }

      population = tempPopulation;
      population.sort();

      // now take best and run HC on messages:
      for(int i=0; i<50; ++i){break;
	std::cout << i << "/" << 5 << "\n";
	CS child = *population.begin();
	child.mutateMessageOnly(mset);
	child.bb84();
	calculateFitness(sim, child, prot, atk, mset);
	population.push_back(child);
      }

      prot->printInfo(sim, *population.begin(), atk);
    }

    GAOutput output;
    output.rate = population.begin()->fitness;
    output.cs = *population.begin();
    return output;
  }



 private:

  int getMessageOpIndex(int ka, int kb, std::vector <int>& m, int t, MessageSettings& mset)
  {
    int index = 0;
    int mult = 1;
    for(int r=0;r<t; ++r){
      index = index + m[r]*mult;
      mult *= mset.dim(r);
    }

    index = index + ka*mult;
    mult *= 2;
    index = index + kb*mult;
    return index;
  }


};
}


