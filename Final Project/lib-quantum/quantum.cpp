/*
quantum.cpp
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#include "quantum.hpp"

// BEGIN QUANTUM NAMESPACE----------------------------------------------------------------------
namespace quantum
{
    // BEGIN KETBRA CLASS ----------------------------------------------------------------------
    Ketbra::Ketbra()
    {
        this->ket = std::vector<int>();
        this->bra = std::vector<int>();
        this->scalar = std::list<int>();
    }

    Ketbra::Ketbra(std::vector <int> &ketIn, std::vector <int> &braIn)
    {
        if(ketIn.size() == braIn.size())
        {
            this->ket = ketIn;
            this->bra = braIn;
            this->scalar = std::list<int>();
        }
        else
        {
            std::cout << "Error: Both Ket and bra must have the same number of indices\n";
        }
    }

    Ketbra::Ketbra(std::list <int> &scalarIn, std::vector <int> &ketIn, std::vector <int> &braIn)
    {
        if(ketIn.size() == braIn.size())
        {
            this->ket = ketIn;
            this->bra = braIn;
            this->scalar = scalarIn;
        }
        else
        {
            std::cout << "Error: Both Ket and bra must have the same number of indices";
        }
    }
    
    bool Ketbra::buildOp(std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn, algebra::mat &output)
    {
        algebra::mat temp;
        
        Complex Prod = Complex(1.0, 0.0);
        
        for (std::list<int>::iterator it=this->scalar.begin(); it != this->scalar.end(); it++)
        {
            Prod = Prod.multiply(scalarIn[*it]);
        }
        
        // TO DO; Prox *= constScalar;

        algebra::mat vecK, vecB;
        vecK.create(1);
        vecB.create(1);
        vecK.ident();
        vecB.ident();
        
        for(int i=0; i < ket.size(); ++i)
        {
            vecK.tensor(&vecIn[i][ket[i]], &temp);
            vecK = temp;
            vecB.tensor(&vecIn[i][bra[i]], &temp);
            vecB = temp;
        }

        vecB.transpose(temp);
        vecK.multiply(temp, output);
        
#ifdef USE_LINUX
	ComplexType ProdT = Prod.getReal() + I*Prod.getImaginary();
	output.multiplyScalar(ProdT);
#else
        output.multiplyScalar(Prod);
#endif   
        return true;
    }
    
    void Ketbra::print()
    {
        std::cout << "[";
        for (std::list<int>::iterator it=this->scalar.begin(); it != this->scalar.end(); it++)
        {
            std::cout << *it;
            if(it != --this->scalar.end())
            {
                std::cout << ",";
            }
        }
        std::cout << "]";
        
        std::cout << "|";
        for (int j=0; j < this->ket.size(); ++j)
        {
            std::cout << this->ket[j];
            if(j != this->bra.size()-1)
            {
                std::cout << ",";
            }
        }
        std::cout << ">";
        
        std::cout << "<";
        for (int k=0; k < this->bra.size(); ++k)
        {
            std::cout << this->bra[k];
            if(k != this->bra.size()-1)
            {
                std::cout << ",";
            }
        }
        std::cout << "|";
    }
    
    void Ketbra::print(std::vector<Complex> &scalarIn)
    {
        Complex C = Complex(1.0, 0.0);
        for (std::list<int>::iterator it=this->scalar.begin(); it != this->scalar.end(); it++)
        {
            C = C.multiply(scalarIn[*it]);
        }
        std::cout << "[";
        C.Complex::print();
        std::cout << "]";
        
        std::cout << "|";
        for (int j=0; j < this->ket.size(); ++j)
        {
            std::cout << this->ket[j];
            if(j != this->bra.size()-1)
            {
                std::cout << ",";
            }
        }
        std::cout << ">";
        
        std::cout << "<";
        for (int k=0; k < this->bra.size(); ++k)
        {
            std::cout << this->bra[k];
            if(k != this->bra.size()-1)
            {
                std::cout << ",";
            }
        }
        std::cout << "|";
    }
    
    // END KETBRA CLASS --------------------------------------------------------------------------
    
    // BEGIN DENSITY OPERATOR CLASS --------------------------------------------------------------
    
    DensityOperator::DensityOperator()
    {
        this->op = std::list<Ketbra>();
    }
    
    void DensityOperator::addKb(Ketbra &KB)
    {
        this->op.push_back(KB);
    }
    
    bool DensityOperator::buildOp(std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn, algebra::mat &output)
    {
        algebra::mat temp, temp2;
        algebra::mat Sum;
        
        int size=1;
        
        for(int i=0; i < vecIn.size(); ++i)
        {
            size = size * vecIn[i][0].getSize();
        }
        
	//std::cout << "DEBUG: size = " << size << "\n";
        Sum.create(size);
        Sum.zero();
        

        for (std::list<Ketbra>::iterator it=this->op.begin(); it != this->op.end(); it++)
        {
            it->buildOp(scalarIn, vecIn, temp);
            
            Sum.add(&temp);
        }
        output = Sum;
        return true;
    }
    
    void DensityOperator::textParser(const std::string& densityOpIn, std::string variablesIn)
    {
        Ketbra::Ketbra Kb2;
        
        std::string input = densityOpIn;
        std::string vars = variablesIn;
        
        wok::Text variables;
        
        variables.addDelimiter(',');
        variables.addIgnore(' ');
        variables.parse(vars);
        
        wok::Text text1, text2;
        
        text1.addDelimiter('*');
        text1.addDelimiter('+');
        
        text1.addIgnore(' ');
        text1.addIgnore('\t');
        
        text1.addCollection('|');
        text1.addCollection(wok::Text::Pair('(', ')'));
        text2.addDelimiter(',');
        text2.addDelimiter('*');
        text2.addDelimiter('>');
        text2.addIgnore(' ');
        text2.addIgnore('\t');
        text2.addIgnore('<');
        
        text1.parse(input);
        
        if(text1.size()%2 != 0)
        {
            std::cout << "Equation error 1\n";
            std::cout << "Text size is: " << text1.size() << std::endl << std::endl;
        }
        else
        {
            for(int i=0; i<text1.size(); i+=2)
            {
                text2.parse(text1[i]);
   
                std::list<int> scal2 = std::list<int>();

                for(int j=0; j<text2.size(); ++j)
                {
                    for(int i=0; i < variables.size(); ++i)
                    {
                        if(text2[j] == variables[i])
                        {
                            scal2.push_back(i);
                        }
                    }
                }
            
                text2.parse(text1[i+1]);
            
                if(text2.size()%2 != 0)
                {
                    std::cout << "Equation error 2, '" << text1[i+1] << "'\n";
                }
                else
                {
                    std::vector<int> ket2 = std::vector<int>();
                    for(int j=0; j<text2.size()/2; ++j)
                    {
                        std::string temp = text2[j];
#ifdef USE_LINUX
			ket2.push_back(atoi(temp.c_str()));
#else
                        ket2.push_back(std::stoi(temp));
#endif
                    }
                    
                    std::vector<int> bra2 = std::vector<int>();
                    for(int j=text2.size()/2; j<text2.size(); ++j)
                    {
                        std::string temp = text2[j];
#ifdef USE_LINUX
			bra2.push_back(atoi(temp.c_str()));
#else
                        bra2.push_back(std::stoi(temp));
#endif
                    }
                    
                    Kb2 = Ketbra(scal2, ket2, bra2);

                    op.push_back(Kb2);
                }
            }
        }
    }
    
    void DensityOperator::print()
    {
        for (std::list<Ketbra>::iterator it=this->op.begin(); it != this->op.end(); it++)
        {
            it->Ketbra::print();
            if(it != --this->op.end())
            {
                std::cout << "+";
            }
        }
    }
    
    void DensityOperator::print(std::vector<Complex> &scalarIn)
    {
        for (std::list<Ketbra>::iterator it=this->op.begin(); it != this->op.end(); it++)
        {
            it->Ketbra::print(scalarIn);
            if(it != --this->op.end())
            {
                std::cout << "+";
            }
        }
    }
    // END DENSITY OPERATOR CLASS ----------------------------------------------------------------
    
    // BEGIN SYSTEM CLASS ------------------------------------------------------------------------
    
    System::System()
    {
        this->vec = std::vector<DensityOperator>();
    }
    
    System::System(std::vector<DensityOperator> vecIn)
    {
        this->vec = vecIn;
    }
    
    bool System::addOp(const std::string equation, std::string vars)
    {
        DensityOperator Temp;
        Temp.textParser(equation, vars);
        this->vec.push_back(Temp);
        return true;
    }
    
    bool System::buildOp(int index, std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn, algebra::mat &output)
    {
        return this->vec[index].buildOp(scalarIn, vecIn, output);
    }
    
    double System::entropy(int index, std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn)
    {

        algebra::mat temp;
	try{
        buildOp(index, scalarIn, vecIn, temp);
	//temp.print(std::cout);
        return temp.entropy();
	}
	catch(std::string& e){
	  std::cout << "FATAL ERROR: " << e << "\n";
	  return 0;
	}
    }
    
    double System::trace(int index, std::vector<Complex> &scalarIn, std::vector<std::vector<algebra::mat> > &vecIn)
    {
        algebra::mat temp;
        buildOp(index, scalarIn, vecIn, temp);
#ifdef USE_LINUX
	return creal(temp.trace());
#else
        return temp.trace().r;
#endif
    }
    
    // END SYSTEM CLASS --------------------------------------------------------------------------
}
// END QUANTUM NAMESPACE--------------------------------------------------------------------------
