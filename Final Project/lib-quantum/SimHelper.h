/*
SimHelper.h
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#pragma once

#include "quantum.hpp"
#include <string>

namespace quantum
{
class SimHelper
{
 public:
  void setup(int numVecCache)
  {
    vecCache.resize(numVecCache);
    text.addIgnore(' ');
    text.addDelimiter(',');
  }
  void setVecCacheDim(int vecCacheIndex, int* dim, int dim_len)
  {
    vecCache[vecCacheIndex].resize(dim_len);
    for(int i=0; i<dim_len; ++i)
      vecCache[vecCacheIndex][i].resize(dim[i]);
  }

  void addDensityOp(const std::string& equation, std::string& scalarList)
  {
    sim.addOp(equation, scalarList);
  }

  void setScalar(std::string& scalarList, const char* scalarToSet, Complex c)
  {
    text.parse(scalarList);
    for(int i=0; i<scalarList.size(); ++i){
      if(text[i].compare(scalarToSet) == 0){
	scalarCache[i] = c;
	return;
      }
    }
    std::cout << "SimHelper: WARNING: Scalar not found.\n";
  }

  void setScalarCache(std::vector <Complex>& c)
  {
    if(scalarCache.size() != c.size())
      scalarCache.resize(c.size());
    for(int i=0; i<scalarCache.size(); ++i)
      scalarCache[i] = c[i];
  }

  void setVecIn(int vecCacheIndex, int subspaceIndex, int vecIndex, algebra::mat& vec)
  {
    vecCache[vecCacheIndex][subspaceIndex][vecIndex] = vec;
  }

  double entropy(int DOIndex, int vecCacheIndex)
  {
    return sim.entropy(DOIndex, scalarCache, vecCache[vecCacheIndex]);
  }
  double trace(int DOIndex, int vecCacheIndex)
  {
    return sim.trace(DOIndex, scalarCache, vecCache[vecCacheIndex]);
  }
  void buildOp(int DOIndex, int vecCacheIndex, algebra::mat& output)
  {
    sim.buildOp(DOIndex, scalarCache, vecCache[vecCacheIndex], output);
  }


 private:
  wok::Text text;
  System sim;
  // vecCache[i][j][k] is a vector for subspace "j", vector "k"
  std::vector < std::vector < std::vector <algebra::mat> > > vecCache;
  std::vector < Complex > scalarCache;

};
}
