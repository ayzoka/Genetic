/*
OneWaySimple.h
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
  class OneWaySimple : public Protocol
  {
  public:
    enum
    {
      PR0 = 0,
      ALPHA0 = 1,
      THETA0 = 2,
      ALPHA1 = 3,
      THETA1 = 4
    };

    virtual ~OneWaySimple()
    {
    }

    virtual void printCSDataTitle(std::ostream& f)
    {
      // p0  a0  t0  a1  t1  |<0|1>|^2
      f << "p0\ta0\tt0\ta1\tt1\t|<0|1>|^2";
    }

    virtual void printCSData(CS& cs, std::ostream& f)
    {
      f << cs.variables[PR0] << "\t";
      f << cs.variables[ALPHA0] << "\t";
      f << cs.variables[THETA0] << "\t";
      f << cs.variables[ALPHA1] << "\t";
      f << cs.variables[THETA1] << "\t";

      quantum::Complex c = cs.variables[ALPHA0]*cs.variables[ALPHA1];
      quantum::Complex alpha0hat = quantum::Complex(sqrt(1-cs.variables[ALPHA0]*cs.variables[ALPHA0])*cos(-cs.variables[THETA0]), sqrt(1-cs.variables[ALPHA0]*cs.variables[ALPHA0])*sin(-cs.variables[THETA0]));

      quantum::Complex alpha1hat = quantum::Complex(sqrt(1-cs.variables[ALPHA1]*cs.variables[ALPHA1])*cos(cs.variables[THETA1]), sqrt(1-cs.variables[ALPHA1]*cs.variables[ALPHA1])*sin(cs.variables[THETA1]));

      c = c + alpha0hat*alpha1hat;

      f << c.getReal()*c.getReal() + c.getImag()*c.getImag();
    }

    virtual void setup(quantum::SimHelper& sim)
    {
      sim.setup(1); // two vec-caches: one for AE the other for E
      int EDim[] = {4};

      sim.setVecCacheDim(0, EDim, 1);

      std::string scalarsUsed = "p0,p1, a0,a0~,a0hat,a0hat~, a1,a1~,a1hat,a1hat~";

      std::stringstream rho00,rho01, rho10, rho11;

      rho00 << "(p0*a0*a0~)*|0><0| + (p0*a0hat*a0hat~)*|2><2| + "
	    << "(p0*a0*a0hat~)*|0><2| + (p0*a0~*a0hat)*|2><0|";

      rho01 << "(p0*a0*a0~)*|1><1| + (p0*a0hat*a0hat~)*|3><3| + "
	    << "(p0*a0*a0hat~)*|1><3| + (p0*a0~*a0hat)*|3><1|";

      rho10 << "(p1*a1*a1~)*|0><0| + (p1*a1hat*a1hat~)*|2><2| + "
	    << "(p1*a1*a1hat~)*|0><2| + (p1*a1~*a1hat)*|2><0|";

      rho11 << "(p1*a1*a1~)*|1><1| + (p1*a1hat*a1hat~)*|3><3| + "
	    << "(p1*a1*a1hat~)*|1><3| + (p1*a1~*a1hat)*|3><1|";

      sim.addDensityOp(rho00.str(), scalarsUsed);
      sim.addDensityOp(rho01.str(), scalarsUsed);
      sim.addDensityOp(rho10.str(), scalarsUsed);
      sim.addDensityOp(rho11.str(), scalarsUsed);


      // setup scalar cache:
      std::vector <quantum::Complex> sc;
      sc.resize(10);

      sim.setScalarCache(sc);

      scalars.resize(10);
    }

    virtual void setupCS(CS& cs)
    {
      // variables are p0, a0, t0, a1, t1:
      cs.variables.resize(5);
      cs.domain.resize(5);

      cs.domain[0] = Domain(0,1);
      cs.domain[1] = Domain(0,1);
      cs.domain[2] = Domain(0,2*3.1415);
      cs.domain[3] = Domain(0,1);
      cs.domain[4] = Domain(0,2*3.1415);
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
      quantum::Complex alpha0 = cs.variables[ALPHA0];
      quantum::Complex alpha0hat = quantum::Complex(sqrt(1-cs.variables[ALPHA0]*cs.variables[ALPHA0])*cos(cs.variables[THETA0]), sqrt(1-cs.variables[ALPHA0]*cs.variables[ALPHA0])*sin(cs.variables[THETA0]));

      quantum::Complex alpha1 = cs.variables[ALPHA1];
      quantum::Complex alpha1hat = quantum::Complex(sqrt(1-cs.variables[ALPHA1]*cs.variables[ALPHA1])*cos(cs.variables[THETA1]), sqrt(1-cs.variables[ALPHA1]*cs.variables[ALPHA1])*sin(cs.variables[THETA1]));

      scalars[0] = cs.variables[PR0];
      scalars[1] = 1-cs.variables[PR0];
      scalars[2] = alpha0;
      scalars[3] = alpha0.conj();
      scalars[4] = alpha0hat;
      scalars[5] = alpha0hat.conj();
      scalars[6] = alpha1;
      scalars[7] = alpha1.conj();
      scalars[8] = alpha1hat;
      scalars[9] = alpha1hat.conj();

      sim.setScalarCache(scalars);

      // set sim vec cache:
      for(int i=0; i<Evec.size(); ++i){
	sim.setVecIn(0, 0, i, Evec[i]);
      }

      sim.buildOp(0,0, rho00);
      sim.buildOp(1,0, rho01);
      sim.buildOp(2,0, rho10);
      sim.buildOp(3,0, rho11);
    }

  private:
    std::vector <quantum::Complex> scalars;
  };
}
