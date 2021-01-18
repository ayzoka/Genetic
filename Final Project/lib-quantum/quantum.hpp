/*
quantum.hpp
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef quantum_hpp
#define quantum_hpp

#include "Complex.hpp"
#include "Algebra.h"
#include "text.h"

#include <iostream>
#include <string>
#include <vector>
#include <list>

namespace quantum
{
    class Ketbra
    {
    public:
        Ketbra();
        // Default Constructor
        
        Ketbra(std::vector <int> &ketIn, std::vector <int> &braIn);
        // Constructor Without Scalar Parameter
        
        Ketbra(std::list <int> &scalarIn, std::vector <int> &ketIn, std::vector <int> &braIn);
        // Constructor With Scalar Parameter
        
        bool buildOp(std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn, algebra::mat &output);
        
        void print();
        // Debug Print
        
        void print(std::vector<Complex> &scalarIn);
        // Print
        
    private:
        std::vector <int> ket;
        std::vector <int> bra;
        std::list <int> scalar;
        Complex constScalar;
    };
    
    class DensityOperator
    {
    public:
        DensityOperator();
        // Default Constructor
        
        void addKb(Ketbra &KB);
        // Add a new Ketbra
        
        bool buildOp(std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn, algebra::mat &output);
        
        void textParser(const std::string& densityOpIn, std::string variablesIn);
        
        void print();
        // Debug Print
        
        void print(std::vector<Complex> &scalarIn);
        //Print
        
    private:
        
        std::list <Ketbra> op;
    };
    
    class System
    {
    public:
        System();
        
        System(std::vector<DensityOperator> vecIn);
        
        bool addOp(const std::string equation, std::string vars);
        
        bool buildOp(int index, std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn, algebra::mat &output);
        
        double entropy(int index, std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn);
        
        double trace(int index, std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn);
        
      private:
        std::vector<DensityOperator> vec;
    };
    
}

#endif /* quantum_hpp */
