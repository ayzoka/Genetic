#pragma once

namespace QKD
{
  class Domain
  {
  private:
    double minVal;
    double maxVal;

  public:
    Domain()
    {
      minVal = 0;
      maxVal = 1;
    }

    Domain(double min, double max)
    {
      minVal = min;
      maxVal = max;
    }

    double fix(double x)
    {
      if (x < minVal)
        return minVal;
      if (x > maxVal)
        return maxVal;
      else
        return x;
    }

    double getMax() { return maxVal; }
    double getMin() { return minVal; }
  };

  // an actual strategy for message passing
  struct MessageStrategy
  {
    std::vector <double> prA, prB;

    double randFloat(double small, double large)
    {
      int r = rand()%RAND_MAX;
      double r2 = r / (double)RAND_MAX;
      
      return (large-small)*r2 + small;
    }

    // doesn't normalize!
    void crossover(MessageStrategy& other, MessageStrategy& child)
    {
      child = *this;
      double p = randFloat(0,1);

      for(int i=0; i<prA.size(); ++i)
	child.prA[i] = p*prA[i] + (1-p)*other.prA[i];

      p = randFloat(0,1);
      for(int i=0; i<prB.size(); ++i)
	child.prB[i] = p*prB[i] + (1-p)*other.prB[i];
    }

    void crossover_ptwise(MessageStrategy& other, MessageStrategy& child)
    {
      child = *this;
      double p = randFloat(0,1);

      if(prA.size() > 0){
	int p1 = rand()%prA.size();
	for(int i=0; i<prA.size(); ++i){
	  if(i < p1)
	    child.prA[i] = prA[i];
	  else
	    child.prA[i] = other.prA[i];
	}
      }

      if(prB.size() > 0){
	int p2 = rand()%prB.size();
	for(int i=0; i<prB.size(); ++i){
	  if(i < p2)
	    child.prB[i] = prB[i];
	  else
	    child.prB[i] = other.prB[i];
	}
      }
    }

    void mutate()
    {
      for(int i=0; i<prA.size(); ++i){
	if(rand()%100 < 50) // 50
	  prA[i] += randFloat(-.75, .75); // .75
      }

      for(int i=0; i<prB.size(); ++i){
	if(rand()%100 < 50)
	  prB[i] += randFloat(-.75, .75);
      }

    }
  };
  
  struct MessageSettings
  {
    MessageSettings() : use_messages(false) {};
    int numRounds() {return M.size();}
    int totalDim(int t)
    {
      int prod=1;
      for(int i=0; i<t; ++i)
	prod *= M[i];
      return prod;
    }
    int totalDim()
    {
      return totalDim(numRounds());
    }
    int dim(int t){return M[t];}

    int getMessageProbIndex(int t, std::vector<int>& m, int key, int i_or_j)
    {
      int index = 0;
      int mult = 1;
      for(int r=0;r<t; ++r){
	index = index + m[r]*mult;
	mult *= dim(r);
      }
      
      index = index + i_or_j*mult;
      mult *= 2;
      index = index + key*mult;
      return index;
    }
    
    bool getNextMessage(std::vector <int> & m)
    {
      for(int r=0; r<m.size(); ++r){
	m[r] ++;
	if(m[r] >= M[r]){
	  m[r] = 0;
	}
	else
	  return true;
      }

      return false;
    }
    bool getNextMessage(std::vector <int> & m, int t)
    {
      for(int r=0; r<t; ++r){
	m[r] ++;
	if(m[r] >= M[r]){
	  m[r] = 0;
	}
	else
	  return true;
      }

      return false;
    }


    
    std::vector <int> M;
    bool use_messages;
  };

  
  struct CS
  {
    void printMessageStrategy(std::ostream& f, MessageSettings& mset)
    {
      // assumes message dim = 0 or 1
      if(messages.size() == 0)
	return;
      std::vector <int> m;
      m.push_back(0);
      int p00 = mset.getMessageProbIndex(0, m, 0, 0);
      int p01 = mset.getMessageProbIndex(0, m, 0, 1);
      int p10 = mset.getMessageProbIndex(0, m, 1, 0);
      int p11 = mset.getMessageProbIndex(0, m, 1, 1);
      f << "p_0^0 = " << messages[0].prA[p00] << "\t";
      f << "p_1^0 = " << messages[0].prA[p10] << "\t";
      f << "p_0^1 = " << messages[0].prA[p01] << "\t";
      f << "p_1^1 = " << messages[0].prA[p11] << "\t";
    }

    double randFloat(double small, double large)
    {
      int r = rand()%RAND_MAX;
      double r2 = r / (double)RAND_MAX;
      
      return (large-small)*r2 + small;
    }
    
    std::vector <double> variables;
    std::vector <Domain> domain;

    std::vector <MessageStrategy> messages;

    void bb84() // only works with OneWaySimple for now!
    {
      return;
      variables[0] = .5;
      return;
      variables[1] = 1;
      variables[2] = 0;
      variables[3] = 0;
      variables[4] = 1;
    }
    
    double fitness;
    
    // ensures that messages[t] is a valid prob. distribution
    void normalizeMessages(MessageSettings& mset)
    {
      bool testMode = false;
      for(int t=0; t<messages.size(); ++t){

	double N = 0;

	////if(t == 0){
	  
	  std::vector <int> m;
	  for(int i=0; i<=t; ++i)
	    m.push_back(0);
	  do{ // for each previous message possibility
	    double Ni0=0, Ni1=0;
	    m[t] = 0;
	    for(; m[t] < mset.M[t]; ++m[t]){
	      //prA:
	      int index = mset.getMessageProbIndex(t+1, m, 0, 0);
	      if(testMode)	    
		messages[t].prA[index] = 1; // test
	      if(messages[t].prA[index] < 0)
		messages[t].prA[index] = 0.01;
	      Ni0 += messages[t].prA[index];
	      
	      index = mset.getMessageProbIndex(t+1, m, 1, 0);
	      if(messages[t].prA[index] <= 0)
		messages[t].prA[index] = 0.01;
	      if(testMode)
		messages[t].prA[index] = 0; // test
	      Ni0 += messages[t].prA[index];
	      
	      index = mset.getMessageProbIndex(t+1, m, 0, 1);
	      if(messages[t].prA[index] < 0)
		messages[t].prA[index] = 0.01;
	      if(testMode)
		messages[t].prA[index] = 0; // test
	      Ni1 += messages[t].prA[index];
	      
	      index = mset.getMessageProbIndex(t+1, m, 1, 1);
	      if(testMode)
		messages[t].prA[index] = 1; // test
	      if(messages[t].prA[index] <= 0)
		messages[t].prA[index] = 0.01;
	      Ni1 += messages[t].prA[index];
	      
	      
	      //prB:
	      for(int j=0; j<=1; ++j){
		int index0 = mset.getMessageProbIndex(t+1, m, 0, j);
		int index1 = mset.getMessageProbIndex(t+1, m, 1, j);
		if(messages[t].prB[index0] < 0)
		  messages[t].prB[index0] = 0;
		if(messages[t].prB[index0] > 1)
		  messages[t].prB[index0] = 1;
		
		// test (4 lines):
		////if(testMode){
		  if(j==0)	   
		    messages[t].prB[index0] = 1;
		  else
		    messages[t].prB[index0] = 0;
		  ////}
		
		messages[t].prB[index1] = 1-messages[t].prB[index0];
	      }
	    }

	    // now normalize pA:
	    m[t] = 0;

	    for(; m[t] < mset.M[t]; ++m[t]){
	      int index = mset.getMessageProbIndex(t+1, m, 0, 0);
	      messages[t].prA[index] /= Ni0;
	      index = mset.getMessageProbIndex(t+1, m, 1, 0);
	      messages[t].prA[index] /= Ni0;
	      index = mset.getMessageProbIndex(t+1, m, 0, 1);
	      messages[t].prA[index] /= Ni1;
	      index = mset.getMessageProbIndex(t+1, m, 1, 1);
	      messages[t].prA[index] /= Ni1;
	    }
	  }while(mset.getNextMessage(m,t)); // (really prev. transcript)
	  /*}
	else{
	  std::vector <int> m;
	  m.resize(t+1);
	  for(int i=0; i<t+1; ++i)
	    m[i] = 0;
	  
	  do{
	    // pb:
	    for(int j=0; j<=1; ++j){
	      int index0 = mset.getMessageProbIndex(t+1, m, 0, j);
	      int index1 = mset.getMessageProbIndex(t+1, m, 1, j);
	      if(messages[t].prB[index0] < 0)
		messages[t].prB[index0] = 0;
	      if(messages[t].prB[index0] > 1)
		messages[t].prB[index0] = 1;
	      
	      messages[t].prB[index1] = 1-messages[t].prB[index0];
	    }
	  }while(mset.getNextMessage(m));

	  // todo: pa
	  }*/
      }
    }

    void randomize(MessageSettings& mset)
    {
      messages.resize(mset.M.size());
      double prodM = 1;
      for(int t=0; t<messages.size(); ++t){
	prodM *= mset.M[t];
	messages[t].prA.resize(4*prodM);
	messages[t].prB.resize(4*prodM);
      }

      for(int t=0; t<messages.size(); ++t){
	for(int i=0; i<messages[t].prA.size(); ++i)
	  messages[t].prA[i] = randFloat(0, 1);

	for(int i=0; i<messages[t].prB.size(); ++i)
	  messages[t].prB[i] = randFloat(0, 1);
      }

      normalizeMessages(mset);

      for(int i=0; i<variables.size(); ++i){
	variables[i] = randFloat(domain[i].getMin(), domain[i].getMax());
      }
    }
    
    CS crossover(CS& other, MessageSettings& mset)
    {
      CS child = *this;
      
      double p = randFloat(0,1);
      
      // child = p*this + (1-p)* other:
      
      for(int i=0; i<variables.size(); ++i){
	child.variables[i] = variables[i]*p + other.variables[i]*(1.0-p);
	child.variables[i] = child.domain[i].fix(child.variables[i]);
      }

      if(messages.size() > 0 && mset.use_messages){
	child.messages = messages;
	for(int t=0; t<messages.size(); ++t)
	  messages[t].crossover(other.messages[t], child.messages[t]);
	
	child.normalizeMessages(mset);
      }

      return child;
    }

    CS crossover_ptwise(CS& other, MessageSettings& mset)
    {
      CS child = *this;
      
      int pt = rand()%variables.size();
      
      // child = p*this + (1-p)* other:
      
      for(int i=0; i<variables.size(); ++i){
	if(i < pt)
	  child.variables[i] = variables[i];
	else
	  child.variables[i] = other.variables[i];

	child.variables[i] = child.domain[i].fix(child.variables[i]);
      }

      if(messages.size() > 0 && mset.use_messages){
	child.messages = messages;
	for(int t=0; t<messages.size(); ++t)
	  messages[t].crossover_ptwise(other.messages[t], child.messages[t]);
	
	child.normalizeMessages(mset);
      }

      return child;
    }

    
    void mutateMessageOnly(MessageSettings& mset)
    {
      if(mset.use_messages){
	for(int t=0; t<messages.size(); ++t)
	  messages[t].mutate();
	
	normalizeMessages(mset);
      }
    }

    void mutate(MessageSettings& mset)
    {
      for(int i=0; i<variables.size(); ++i){
	if(randFloat(0,1) < .4){ // .4
	  variables[i] += randFloat(-1.25, 1.25);
	  
	  variables[i] = domain[i].fix(variables[i]);
	}
      }

      if(mset.use_messages){
	for(int t=0; t<messages.size(); ++t)
	  messages[t].mutate();
	
	normalizeMessages(mset);
      }
    }
    
    bool operator< (const CS& other) const
    {
      if(fitness > other.fitness)
	return true;
      return false;
    }
  };
}
