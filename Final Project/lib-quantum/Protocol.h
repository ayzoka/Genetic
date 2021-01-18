/*
Protocol.h
Written by Walter O. Krawec, Michael Nelson, and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#pragma once
#include <fstream>

#include "../CS.h"
#include "../attack.h"
#include "../lib-quantum/SimHelper.h"

namespace QKD
{
class Protocol
{
 public:
  virtual ~Protocol()
  {
  }

  //virtual double keyRate(quantum::SimHelper& sim, CS& cs, Attack& atk)=0;

  virtual void constructEOps(quantum::SimHelper& sim, CS& cs,
			     std::vector <algebra::mat>& Evec, 
			     algebra::mat& rho00, algebra::mat& rho01,
			     algebra::mat& rho10, algebra::mat& rho11)=0;

  virtual void setup(quantum::SimHelper& sim)=0;

  virtual void setupCS(CS& cs)=0;

  virtual void printInfo(quantum::SimHelper& sim, CS& cs, Attack& atk)=0;

  virtual void printCSDataTitle(std::ostream& f)=0;
  virtual void printCSData(CS& cs, std::ostream& f)=0;

};
}
