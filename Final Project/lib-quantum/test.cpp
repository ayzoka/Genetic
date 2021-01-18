/*
test.cpp
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#include "lib-quantum/SimHelper.h"

int main(int argc, char** argv)
{
  quantum::SimHelper sim;
  sim.setup(1);
  int vecDim[] = {4};
  sim.setVecCacheDim(0, vecDim, 1);

  std::string eq = "(oh)*|0><0| + (oh)*|2><2|";
  std::string sc = "oh,oq,tq";

  sim.addDensityOp(eq, sc);


  std::vector <quantum::Complex> scalars;
  scalars.push_back( quantum::Complex(.5,0) );
  scalars.push_back( quantum::Complex(.25,0) );
  scalars.push_back( quantum::Complex(.75,0) );


  sim.setScalarCache(scalars);


  algebra::mat V;
  V.create(2,1);
  V.ident();
  V.print(std::cout);

  sim.setVecIn(0, 0, 0, V);

  //V.zero();
  //V(1,0) = 1;
  sim.setVecIn(0, 0, 2, V);


  std::cout << "COMPUTE:\n";
  std::cout << "E = " << sim.entropy(0,0) << "\n";
  std::cout << "T = " << sim.trace(0,0) << "\n";

  sim.setScalar(sc, "oh", quantum::Complex(.25,0));
  std::cout << "T = " << sim.trace(0,0) << "\n";
  return 0;
}
