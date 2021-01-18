/*
OneWayGeneral.h
Written by Walter O. Krawec, Michael Nelson, and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#pragma once

#include "Protocol.h"
#include <sstream>

namespace QKD
{
  class OneWayGeneral : public Protocol
  {
  public:
    enum
    {
      PR0 = 0,
      ALPHA0 = 1,
      THETA0 = 2,
      ALPHA1 = 3,
      THETA1 = 4,

      BETA0 = 5,
      PHI0 = 6,
      BETA1 = 7,
      PHI1 = 8
    };

    virtual ~OneWayGeneral()
    {
    }

    virtual void printCSDataTitle(std::ostream& f)
    {
      // p0  a0  t0  a1  t1  |<0|1>|^2
      f << "p0\ta0\tt0\ta1\tt1\tb0\tphi0\tb1\tphi1\n";
    }

    virtual void printCSData(CS& cs, std::ostream& f)
    {
      f << cs.variables[PR0] << "\t";
      f << cs.variables[ALPHA0] << "\t";
      f << cs.variables[THETA0] << "\t";
      f << cs.variables[ALPHA1] << "\t";
      f << cs.variables[THETA1] << "\t";

      f << cs.variables[BETA0] << "\t";
      f << cs.variables[PHI0] << "\t";
      f << cs.variables[BETA1] << "\t";
      f << cs.variables[PHI1] << "\t";
    }



    virtual void setup(quantum::SimHelper& sim)
    {
      sim.setup(1); // two vec-caches: one for AE the other for E
      int EDim[] = {4};

      sim.setVecCacheDim(0, EDim, 1);

      std::string scalarsUsed = "p0,p1";

      std::stringstream rho00,rho01, rho10, rho11;

      // these represent "g" states...
      rho00 << "(p0)*|0><0|";
      rho01 << "(p0)*|1><1|";
      rho10 << "(p1)*|2><2|";
      rho11 << "(p1)*|3><3|";

      sim.addDensityOp(rho00.str(), scalarsUsed);
      sim.addDensityOp(rho01.str(), scalarsUsed);
      sim.addDensityOp(rho10.str(), scalarsUsed);
      sim.addDensityOp(rho11.str(), scalarsUsed);


      // setup scalar cache:
      std::vector <quantum::Complex> sc;
      sc.resize(2);

      sim.setScalarCache(sc);

      scalars.resize(2);

      Fvec.resize(4);
      Gvec.resize(4);

      for(int i=0; i<4; ++i){
	Fvec[i].create(4,1);
	Gvec[i].create(4,1);
      }
    }

    void constructFvec(CS& cs, std::vector <algebra::mat>& Evec)
    {
      quantum::Complex alpha0 = cs.variables[ALPHA0];
      quantum::Complex alpha0hat = quantum::Complex(sqrt(1-cs.variables[ALPHA0]*cs.variables[ALPHA0])*cos(cs.variables[THETA0]), sqrt(1-cs.variables[ALPHA0]*cs.variables[ALPHA0])*sin(cs.variables[THETA0]));

      quantum::Complex alpha1 = cs.variables[ALPHA1];
      quantum::Complex alpha1hat = quantum::Complex(sqrt(1-cs.variables[ALPHA1]*cs.variables[ALPHA1])*cos(cs.variables[THETA1]), sqrt(1-cs.variables[ALPHA1]*cs.variables[ALPHA1])*sin(cs.variables[THETA1]));

      Fvec[0].zero();
      Fvec[0].add(alpha0, Evec[0]);
      Fvec[0].add(alpha0hat, Evec[2]);

      Fvec[1].zero();
      Fvec[1].add(alpha0, Evec[1]);
      Fvec[1].add(alpha0hat, Evec[3]);

      Fvec[2].zero();
      Fvec[2].add(alpha1, Evec[0]);
      Fvec[2].add(alpha1hat, Evec[2]);

      Fvec[3].zero();
      Fvec[3].add(alpha1, Evec[1]);
      Fvec[3].add(alpha1hat, Evec[3]);
    }

    // assumes Fvec already constructed...
    void constructGvec(CS& cs)
    {
      // phi_ij = <phi_i | j>:
      quantum::Complex phi00 = cs.variables[BETA0];
      quantum::Complex phi01 = quantum::Complex(sqrt(1-cs.variables[BETA0]*cs.variables[BETA0])*cos(-cs.variables[PHI0]), sqrt(1-cs.variables[BETA0]*cs.variables[BETA0])*sin(-cs.variables[PHI0]));

      quantum::Complex phi10 = cs.variables[BETA1];
      quantum::Complex phi11 = quantum::Complex(sqrt(1-cs.variables[BETA1]*cs.variables[BETA1])*cos(-cs.variables[PHI1]), sqrt(1-cs.variables[BETA1]*cs.variables[BETA1])*sin(-cs.variables[PHI1]));

      Gvec[0].zero();
      Gvec[0].add(phi00, Fvec[0]);
      Gvec[0].add(phi01, Fvec[1]);

      Gvec[1].zero();
      Gvec[1].add(phi10, Fvec[0]);
      Gvec[1].add(phi11, Fvec[1]);

      Gvec[2].zero();
      Gvec[2].add(phi00, Fvec[2]);
      Gvec[2].add(phi01, Fvec[3]);

      Gvec[3].zero();
      Gvec[3].add(phi10, Fvec[2]);
      Gvec[3].add(phi11, Fvec[3]);
    }

    virtual void setupCS(CS& cs)
    {
      // variables are p0, a0, t0, a1, t1, b0, phi0, b1, phi1:
      cs.variables.resize(9);
      cs.domain.resize(9);

      cs.domain[0] = Domain(0,1);
      cs.domain[1] = Domain(0,1);
      cs.domain[2] = Domain(0,2*3.1415);
      cs.domain[3] = Domain(0,1);
      cs.domain[4] = Domain(0,2*3.1415);

      cs.domain[5] = Domain(0,1);
      cs.domain[6] = Domain(0,2*3.1415);
      cs.domain[7] = Domain(0,1);
      cs.domain[8] = Domain(0,2*3.1415);

    }

    virtual void printInfo(quantum::SimHelper& sim, CS& cs, Attack& atk)
    {
      double kr = cs.fitness;
      std::cout << "key-rate = " << kr << "\n";
      //std::cout << "\tp0 = " << cs.variables[0] << "\n";
    }

    virtual void constructEOps(quantum::SimHelper& sim, CS& cs,
			       std::vector <algebra::mat>& Evec, 
			       algebra::mat& rho00, algebra::mat& rho01,
			       algebra::mat& rho10, algebra::mat& rho11)
    {
      scalars[0] = cs.variables[PR0];
      scalars[1] = 1-cs.variables[PR0];

      sim.setScalarCache(scalars);

      constructFvec(cs, Evec);
      constructGvec(cs);

      // set sim vec cache:
      for(int i=0; i<Gvec.size(); ++i){
	//Evec[i].print(std::cout);
	sim.setVecIn(0, 0, i, Gvec[i]);
      }

      sim.buildOp(0,0, rho00);
      sim.buildOp(1,0, rho01);
      sim.buildOp(2,0, rho10);
      sim.buildOp(3,0, rho11);

      // normalize:
      double N = rho00.traceRe() + rho01.traceRe() + rho10.traceRe() + rho11.traceRe();

      rho00.multiplyScalar(1.0/N);
      rho01.multiplyScalar(1.0/N);
      rho10.multiplyScalar(1.0/N);
      rho11.multiplyScalar(1.0/N);
    }

  private:
    std::vector <quantum::Complex> scalars;
    std::vector <algebra::mat> Fvec, Gvec;
  };
}
