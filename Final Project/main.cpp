#include "lib-quantum/SimHelper.h"

#include "lib-quantum/OneWaySimple.h"
#include "lib-quantum/OneWayGeneral.h"
#include "GA.h"

// mode = 0 means sym channel; else random (asym)
void testSym(int mode, double Q, int message_dim, QKD::Protocol* prot, const char* filename, int num_trials, int num_gen, int x_type, int mut_rate)
{
  srand(time(NULL));

  QKD::GAInput input;
  input.mutation_rate = mut_rate;
  input.crossover_type = x_type;
  input.tournament_size = 5; // 5
  input.population_size = 100;
  input.max_iterations = num_gen;

  QKD::MessageSettings mset;

  mset.use_messages = false;
  mset.M.push_back(1);

  if(message_dim > 0){
    mset.use_messages = true;
    mset.M[0] = message_dim;
  }

  QKD::Stats stats;

  
  // int t = time(NULL);
  // t = 1483114561; // asym1
  // t = 1483375172;
  // t = 1483545017;
  // t = 1483546159; // asym3
  // std::cout << "seed = " << t << "\n";
  // srand(t);
  
  // Q = .3;
  // stats.generateRandom(Q,Q);


  if(mode == 0)
    stats.symmetric(Q);
  else{
    srand(Q);
    Q = .25;  // changing this for different channels...
    stats.generateRandom(Q,Q);
  }

  srand(time(NULL));

  std::cout << "\n";
  std::ofstream f;
  f.open(filename, std::ios_base::app);
  stats.printChannel(f);
  stats.printChannel(std::cout);
  f << "\n\n";

  f << "rate\t";
  prot->printCSDataTitle(f);
  f << "\n";
  f.close();

  std::cout << "BB84 Key Rate = " << stats.BB84Rate() << "\n\n";

  for(int i=0; i<num_trials; ++i){
    std::cout << "TRIAL " << i << "/" << num_trials << "\n";
    QKD::Attack atk(stats);

    QKD::EvolveQKD qkd;
    
    try{
      QKD::GAOutput output = qkd.evolve(prot, atk, mset, input);
      f.open(filename, std::ios_base::app);
      f << output.rate << "\t";

      prot->printCSData(output.cs, f);

      output.cs.printMessageStrategy(f, mset);

      f << "\n";
      f.close();
    }catch(std::string& e){
      std::cout << "ERROR: " << e << "\n";
    }
  }
}

// void printUsage()
// {
//   std::cout << "Usage: ./QKDOpt [mode] [noise/seed] [prot-type] [message-dim] [num-trials] [num-gen] [x-type] [mutation-rate] [file-name]\n";

//   std::cout << "\t[mode]: 'sym' or 'asym' for symmetric or random channel\n";
//   std::cout << "\t[noise/seed]: the amount of symmetric noise in the channel\n";
//   std::cout << "\t\tor the seed to use for a random channel.\n";
//   std::cout << "\t[prot-type]: 1 = OneWaySimple; 2 = OneWayGeneral\n";
//   std::cout << "\t[message-dim]: how many messages allowed? 0 = no preproc.\n";
//   std::cout << "\t[num-trials]: number of independent trials to run.\n";
//   std::cout << "\t[num-gen]: number of generations per trial\n";
//   std::cout << "\t[x-type]: 1=blend, 2=one-point\n";
//   std::cout << "\t[mutation-rate]: probability of mutation (sug. 75).\n\n";
//   std::cout << "\t[file-name]: file to save data (WILL BE OVERWRITTEN!)\n";
// }

int main(int argc, char** argv)
{
  srand(time(NULL));

  if(argc != 10){
    // printUsage();
    return 0;
  }

  const char* mode = argv[1];
  double Q = atof(argv[2]);
  int prot_type = atoi(argv[3]);
  int message_dim = atoi(argv[4]);
  int num_trials = atoi(argv[5]);
  int num_gen = atoi(argv[6]);
  int x_type = atoi(argv[7]);
  double mut_rate = atof(argv[8]);
  const char* filename = argv[9];

  int M = 0; // 0 for sym; 1 for random
  if(mode[0] == 's')
    M = 0;
  else if(mode[0] == 'a')
    M = 1;
  else
    std::cout << "Warning: Mode not recognized, using sym\n";

  std::ofstream f;
  f.open(filename);
  for(int i=0; i<argc; ++i)
    f << argv[i] << " ";
  f.close();

  QKD::Protocol* prot;
  QKD::OneWaySimple ows;
  QKD::OneWayGeneral owg;

  if(prot_type == 1)
    prot = &ows;
  else if(prot_type == 2)
    prot = &owg;
  else{
    std::cout << "ILLEGAL PROTOCOL TYPE (should be 1 or 2).\n";
    return 0;
  }
  testSym(M, Q, message_dim, prot, filename, num_trials, num_gen, x_type, mut_rate);

  return 0;
}
