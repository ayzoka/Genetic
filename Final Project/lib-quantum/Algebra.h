/*
Algebra.h (Mac Version)
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _ALGEBRA_H_
#define _ALGEBRA_H_

#define USE_LAPACK
#define USE_MAC

// remove this later to reactivate complex mats
// but only works in Linux, so keep for testing

//#define USE_REAL

// #include <Accelerate/Accelerate.h> - To use lapack
// #include "Complex.h" - If Our Complex CLass is needed

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <complex.h>

#include "Complex.hpp"

#ifdef USE_LINUX
#define Complex double _Complex
#endif

#ifdef USE_MAC
#include <Accelerate/Accelerate.h>
#endif

#define ComplexType __CLPK_doublecomplex

namespace algebra
{
    // this complex helper business only works on Linux
    // if on Mac, keep USE_REAL defined
#ifdef USE_LINUX
    class ComplexHelper
    {
    public:
        static Complex one() {Complex out; out = 1.0 + 0*I; return out; }
        static Complex zero() {Complex out; out = 0 + 0*I; return out;}
        
        static void sq(Complex& v)
        {
            Complex v2 = creal(v) - cimag(v)*I;
            v = mult(v, v2);
        }
        
        static Complex construct(double theta)
        {
            Complex out = cos(theta) + sin(theta)*I;
            return out;
        }
        static Complex power(Complex& base, int e)
        {
            return cpow(base, e);
            Complex a;
            double real = log( sqrt( creal(base)*creal(base)* + cimag(base)*cimag(base)));
            double imag = atan2(cimag(base), creal(base));
            
            real *= e;
            imag *= e;
            
            a = real + imag*I;
            
            Complex result = construct(cimag(a));
            
            Complex r2 = creal(result)*exp(real) + cimag(result)*real*I;
            
            return r2;
        }
        
        static Complex add(Complex& a, Complex& b)
        {
            return a+b;
        }
        
        static Complex subtract(Complex& a, Complex& b)
        {
            return a-b;
        }
        
        static Complex mult(Complex& a, Complex& b)
        {
            return a*b;
        }
        
        static Complex invert(Complex& a)
        {
            Complex out;
            double d = creal(a)*creal(a) + cimag(a)*cimag(a);
            
            out = creal(a) / d - cimag(a)/d *I;
            
            return out;
        }
    };
#endif
    
    // actual matrix class (also for vectors)
    class matrix
    {
    public:
        matrix()
        {
            elements = NULL;
            noCleanup = false;
            r = 0;
            c = 0;
            n = 0;
        }
        
        matrix(const char* _name)
        {
            elements = NULL;
            noCleanup = false;
            r = 0;
            c = 0;
            n = 0;
            varName = std::string(_name);
        }
        
        matrix(const matrix& m)
        {
            elements = NULL;
            r = m.r;
            c = m.c;
            n = r*c;
            
            if(r == 0 || c == 0 || !m.elements)
                return;
#ifdef USE_LINUX
            elements = new double[n];
#elif defined(USE_MAC)
            elements = new __CLPK_doublecomplex[n];
#else
            elements = new Complex[n];
#endif
            noCleanup = false;
            
            for(int i=0; i<n; ++i)
                elements[i] = m.elements[i];
        }
        
        ~matrix() {cleanup();}
        
        // sets this to be "m" correctly (i.e., a copy is made)
        matrix& operator=(const matrix& m)
        {
            if(!elements || r != m.r || c != m.c)
            {
                cleanup();
                r = m.r;
                c = m.c;
                n = r*c;
                
                if(r == 0 || c == 0 || !m.elements)
                {
                    throw std::string("operator= error.\n");
                    return *this;
                }
#ifdef USE_LINUX
                elements = new double[n];
#elif defined(USE_MAC)
                elements = new __CLPK_doublecomplex[n];
#else
                elements = new Complex[n];
#endif
            }
            
            if(!m.elements)
            {
                throw std::string("operator= error 2.\n");
                return *this;
            }
            
            for(int i=0; i<n; ++i)
                elements[i] = m.elements[i];
            
            noCleanup = false;
            
            return *this;
        }
        
        // destroy the matrix (also done automatically)
        void cleanup()
        {
            if(elements && !noCleanup)
            {
                delete[] elements;
            }
            elements = NULL;
            r = 0;
            c = 0;
            n = 0;
        }
        
        // sets this to be a zero matrix
        void zero()
        {
            if(!elements)
                return;
            
            for(int i=0; i<n; ++i)
            {
#ifdef USE_LINUX
                elements[i] = 0;
#elif defined(USE_MAC)
                elements[i].r = 0;
                elements[i].i = 0;
#else
                elements[i] = 0.0+0.0*I;
#endif
            }
        }
        
        // sets this to be the identity matrix
        void ident()
        {
            if(!elements)
                return;
            
            zero();
            int dim = getRows();
            if(getCols() < getRows())
                dim = getCols();
#ifdef USE_LINUX
            double one = 1;
#elif defined(USE_MAC)
            __CLPK_doublecomplex one;
            one.r = 1;
            one.i = 0;
#else
            Complex one = 1.0+0.0*I;
#endif
            for(int i=0; i<dim; ++i)
                this->setDirect(i,i,one);
        }
        
        // compute eigenvalues (needs LAPACK!)
        bool eigenValues(std::vector <double>& values);
        
        // set ouput = sum_i pi*Mi
        static void addMulti(matrix& M1, double p1, matrix& M2, double p2, matrix& M3, double p3, matrix& M4, double p4, matrix& output)
        {
            // todo: check sizes correct
            output.create(M1.getRows(), M1.getCols());
            for(int i=0; i<M1.n; ++i){
                output.elements[i].r = p1*M1.elements[i].r + p2*M2.elements[i].r
                + p3*M3.elements[i].r + p4*M4.elements[i].r;
                output.elements[i].i = p1*M1.elements[i].i + p2*M2.elements[i].i
                + p3*M3.elements[i].i + p4*M4.elements[i].i;
            }
        }
        
        double traceRe()
        {
            __CLPK_doublecomplex tr = trace();
            return tr.r;
        }
        
        // compute the trace of the matrix
#ifdef USE_LINUX
        double trace()
#elif defined(USE_MAC)
        __CLPK_doublecomplex trace()
#else
        Complex trace()
#endif
        {
#ifdef USE_LINUX
            double c = 0;
#elif defined(USE_MAC)
            __CLPK_doublecomplex c;
            c.r = 0;
            c.i = 0;
#else
            Complex c = 0.0 + 0.0*I;
#endif
            for(int i=0; i<getRows(); ++i)
            {
#ifdef USE_MAC
                c.r = c.r + (*this)(i,i).r;
                c.i = c.i + (*this)(i,i).i;
#else
                c = c + (*this)(i,i);
#endif
            }
            return c;
        }
        
        // set this to be all elements = value
#ifdef USE_LINUX
        void clear(double value)
#elif defined(USE_MAC)
        void clear(__CLPK_doublecomplex value)
#else
        void clear(Complex value)
#endif
        {
            if(!elements)
                return;
            for(int i=0; i<n; ++i)
                elements[i] = value;
        }
        
        // create a matrix with _r rows and _c cols
        bool create(int _r, int _c)
        {
            if(_r <= 0 || _c <= 0)
                return false;
            if(elements && r == _r && c == _c)
                return true;
            
            cleanup();
            r = _r;
            c = _c;
            n = r*c;
            
#ifdef USE_LINUX
            elements = new double[r*c];
#elif defined(USE_MAC)
            elements = new __CLPK_doublecomplex[r*c];
#else
            elements = new Complex[r*c];
#endif
            
            noCleanup = false;
            
            return true;
        }
        
        // create a square _r by _r matrix
        bool create(int _r)
        {
            return create(_r, _r);
        }
        
        // get reference to element in row "row" and col "col"
#ifdef USE_LINUX
        double& operator() (int row, int col)
#elif defined(USE_MAC)
        __CLPK_doublecomplex& operator() (int row, int col)
#else
        Complex& operator() (int row, int col)
#endif
        {
            int _r = row;
            int _c = col;
            if(!elements || _r >= r || _c >= c || _r < 0 || _c < 0)
            {
                std::cout << "Error 5 (" << elements << ", " << _r << ", " << r << ", " << _c << ", " << c << ")\n";
                std::cout << "Name = " << varName << "\n";
                throw(std::string("Algebra Error 5.1"));
            }
            return elements[_c*r + _r];
        }
        
        // similar to above, but not a reference
#ifdef USE_LINUX
        double getDirect(int row, int col)
#elif defined(USE_MAC)
        __CLPK_doublecomplex getDirect(int row, int col)
#else
        Complex getDirect(int row, int col)
#endif
        {
            int _r = row;
            int _c = col;
            if(!elements || _r >= r || _c >= c || _r < 0 || _c < 0)
            {
                std::cout << "Error 5 (" << elements << ", " << _r << ", " << r << ", " << _c << ", " << c << ")\n";
                throw(std::string("Algebra Error 5.2"));
            }
            
            return elements[_c*r + _r];
        }
        
        // similar to above, but not a reference
#ifdef USE_LINUX
        void setDirect(int _n, double data)
#elif defined(USE_MAC)
        void setDirect(int _n, __CLPK_doublecomplex data)
#else
        void setDirect(int _n, Complex data)
#endif
        {
            if(!elements || _n >= n)
            {
                throw(std::string("Algebra Error 5.3"));
            }
            
            elements[_n] = data;
        }
        
#ifdef USE_LINUX
        void setDirect(int row, int col, double data)
#elif defined(USE_MAC)
        void setDirect(int row, int col, __CLPK_doublecomplex data)
#else
        void setDirect(int row, int col, Complex data)
#endif
        {
            int _r = row;
            int _c = col;
            if(!elements || _r >= r || _c >= c || _r < 0 || _c < 0)
            {
                std::cout << "Error 5 (" << elements << ", " << _r << ", " << r << ", " << _c << ", " << c << ")\n";
                throw(std::string("Algebra Error 5.4"));
            }
            elements[_c*r + _r] = data;
        }
        
        void setDirect(int row, int col, double re, double im)
        {
            int _r = row;
            int _c = col;
            if(!elements || _r >= r || _c >= c || _r < 0 || _c < 0)
            {
                std::cout << "Error 5 (" << elements << ", " << _r << ", " << r << ", " << _c << ", " << c << ")\n";
                throw(std::string("Algebra Error 5.4"));
            }
            
#ifdef USE_LINUX
            elements[_c*r + _r] = re;
#elif defined(USE_MAC)
            elements[_c*r + _r].r = re;
            elements[_c*r + _r].i = im;
#else
            elements[_c*r + _r] = re + im*I;
#endif
        }
        
#ifdef USE_LINUX
        double get(int _n)
#elif defined(USE_MAC)
        __CLPK_doublecomplex get(int _n)
#else
        Complex get(int _n)
#endif
        {
            if(!elements || _n >= n)
            {
                std::cout << "Error 2 (" << elements << ", " << _n << ", " << n << ")\n";
                throw std::string("error 2.");
            }
            
            return elements[_n];
        }
        
#ifdef USE_LINUX
        double* getElements()
#elif defined(USE_MAC)
        __CLPK_doublecomplex* getElements()
#else
        Complex* getElements()
#endif
        {
            return elements;
        }
        
        // prints to a stream (e.g., std::cout or a file)
        void print(std::ostream& f, bool printAll=false)
        {
            if(!elements)
                return;
            int rStart = 0, cStart=0;
            int rows = alg_min(10, r);
            int cols = alg_min(10, c);
            
            if(printAll)
            {
                rows = r;
                cols = c;
            }
            for(int _r=rStart; _r < rows; ++_r)
            {
                for(int _c=cStart; _c<cols; ++_c)
                {
#ifdef USE_LINUX
                    f << (*this)(_r, _c) << "\t";
#elif defined(USE_MAC)
                    f << (*this)(_r, _c).r << " + " << (*this)(_r, _c).i << "I\t";
                    
#else
                    f << creal((*this)(_r,_c)) << "+" << cimag((*this)(_r,_c)) << "i\t";
#endif
                }
                f << "\n";
            }
            f << "\n";
        }
        
        /** TODO For USE_MAC, cast to quantum::Complex, multiply, then cast back **/
#ifdef USE_MAC
        void multiplyScalar(quantum::Complex &c)
        {
            // may need USE_MAC
            // // like this with complex multiplication elements[i] = elements[i] * c;
            for( int i = 0; i <n; ++i)
            {
                quantum::Complex Temp = quantum::Complex(elements[i].r, elements[i].i);
                Temp = Temp.multiply(c);
                elements[i].r = Temp.getReal();
                elements[i].i = Temp.getImaginary();
            }
        }
        
        void multiplyScalar(double c)
        {
            // may need USE_MAC
            // // like this with complex multiplication elements[i] = elements[i] * c;
            for( int i = 0; i <n; ++i)
            {
                elements[i].r *= c;
                elements[i].i *= c;
            }
        }
        
#endif
        // multiply all elements by given scalar
#ifdef USE_LINUX
        void multiplyScalar(double c)
#elif defined(USE_MAC)
        void multiplyScalar(__CLPK_doublecomplex& c)
#else
        void multiplyScalar(Complex& c)
#endif
        {
            if(!elements)
                return;
            
            quantum::Complex cTemp = quantum::Complex(c.r, c.i);
            for(int i=0; i<n; ++i)
            {
                ////elements[i] = elements[i] * c;
                //elements[i] = elements[i].multiply(c);
                quantum::Complex Temp = quantum::Complex(elements[i].r, elements[i].i);
                Temp = Temp.multiply(cTemp);
                elements[i].r = Temp.getReal();
                elements[i].i = Temp.getImaginary();
            }
            
        }
        
        
        // transpose
        void sqTranspose()
        {
            if(!elements || c != r)
            {
                throw(std::string("Algebra error 1.2"));
                return;
            }
            
#ifdef USE_LINUX
            double temp;
#elif defined(USE_MAC)
            __CLPK_doublecomplex temp;
#else
            Complex temp;
#endif
            for(long row=0; row<r; ++row)
            {
                for(long col=row; col<c; ++col)
                {
                    if(row == col)
                        continue;
                    temp = getDirect(row,col);
                    setDirect(row,col, getDirect(col,row));
                    setDirect(col,row, temp);
                }
            }
        }
        
        /** TODO: for USE_MAC **/
        
        
        // conjugate transpose
        void transpose(matrix& B)
        {
            if(!elements)
                return;
            B.create(c,r);
            for(int _r=0; _r<r; ++_r)
            {
                for(int _c=0; _c<c; ++_c)
                {
#ifdef USE_LINUX
                    double c = B(_r, _c);
                    B(_c, _r) = c;
#elif defined(USE_MAC)
                    // TODO
                    //quantum::Complex c = B(_r, _c);
                    //B(_c, _r) = creal(c).subtract(cimag(c));
                    quantum::Complex C = quantum::Complex((*this)(_r, _c).r, (*this)(_r, _c).i);
                    C = C.conjugate();
                    B(_c, _r).r = C.getReal();
                    B(_c, _r).i = C.getImaginary();
#else
                    Complex c = B(_r, _c);
                    B(_c, _r) = creal(c) - cimag(c)*I;
#endif
                }
            }
        }
        
        
        /** TODO USE_MAC **/
        
        
#ifdef USE_LINUX
        void subtract(double c)
#elif defined(USE_MAC)
        void subtract(__CLPK_doublecomplex& c)
#else
        void subtract(Complex& c)
#endif
        {
            if(!elements)
                return;
            
            quantum::Complex cTemp = quantum::Complex(c.r, c.i);
            for(int i=0; i<n; ++i)
            {
                ////elements[i] = elements[i] - c;
                //elements[i] = elements[i].subtract(c);
                quantum::Complex Temp = quantum::Complex(elements[i].r, elements[i].i);
                Temp = Temp.subtract(cTemp);
                elements[i].r = Temp.getReal();
                elements[i].i = Temp.getImaginary();
            }
            
        }
        
        void subtract(matrix* B, matrix* output)
        {
            if(!B || !output)
                return;
            
            if(!elements || !B->getElements() || getRows() != B->getRows() || getCols() != B->getCols())
            {
                throw std::string("Error - cannot subtract matrix\n");
                return;
            }
            
            output->create(getRows(), getCols());
            
            for(int i=0; i<n; ++i)
            {
                ////output->elements[i] = elements[i] - B->getElements()[i];
                //output->elements[i] = elements[i].subtract(B->getElements()[i];
                quantum::Complex Temp = quantum::Complex(elements[i].r, elements[i].i);
                quantum::Complex BTemp = quantum::Complex(B->getElements()[i].r, B->getElements()[i].i);
                Temp = Temp.subtract(BTemp);
                output->elements[i].r = Temp.getReal();
                output->elements[i].i = Temp.getImaginary();
            }
        }
        
        void add(quantum::Complex& c2, matrix& B)
        {
#ifdef USE_LINUX
            ComplexType c = c2.getReal() + I*c2.getImaginary();
#endif
            matrix* output = this;
            
            if(!elements || !B.getElements() || getRows() != B.getRows() || getCols() != B.getCols())
            {
                throw std::string("Error - cannot add matrix\n");
                return;
            }
            
            //output->create(getRows(), getCols());
            
            for(int i=0; i<n; ++i)
            {
                quantum::Complex z(B.getElements()[i].r, B.getElements()[i].i);;
                z = c2.multiply(z);
                output->elements[i].r += z.getReal();
                output->elements[i].i += z.getImaginary();
            }
        }
        
        void add(matrix* B, matrix* output)
        {
            if(!B || !output)
                return;
            
            if(!elements || !B->getElements() || getRows() != B->getRows() || getCols() != B->getCols())
            {
                std::cout << "Error 5 (" << elements << ", " << getRows() << ", " << B->getRows() << ", " << getCols() << ", " << B->getCols() << ")\n";
                throw std::string("Error - cannot add matrix\n");
                return;
            }
            
            output->create(getRows(), getCols());
            
            for(int i=0; i<n; ++i)
            {
                ////output->elements[i] = elements[i] + B->getElements()[i];
                //output->elments[i] = elements[i].add(B->getElements()[i];
                quantum::Complex Temp = quantum::Complex(elements[i].r, elements[i].i);
                quantum::Complex BTemp = quantum::Complex(B->getElements()[i].r, B->getElements()[i].i);
                Temp = Temp.add(BTemp);
                output->elements[i].r = Temp.getReal();
                output->elements[i].i = Temp.getImaginary();
            }
        }
        
        
        void add(matrix* B)
        {
            matrix* output = this;
            if(!B)
                return;
            
            if(!elements || !B->getElements() || getRows() != B->getRows() || getCols() != B->getCols())
            {
                std::cout << "Error 5 (" << elements << ", " << getRows() << ", " << B->getRows() << ", " << getCols() << ", " << B->getCols() << ")\n";
                throw std::string("Error - cannot add matrix\n");
                return;
            }
            
            
            for(int i=0; i<n; ++i)
            {
                ////output->elements[i] = elements[i] + B->getElements()[i];
                //output->elments[i] = elements[i].add(B->getElements()[i];
                quantum::Complex Temp = quantum::Complex(elements[i].r, elements[i].i);
                quantum::Complex BTemp = quantum::Complex(B->getElements()[i].r, B->getElements()[i].i);
                Temp = Temp.add(BTemp);
                elements[i].r = Temp.getReal();
                elements[i].i = Temp.getImaginary();
            }
        }
        
        
#ifdef USE_LINUX
        void dot(matrix* B, Complex& output);	// returns this <dot> B
#endif
        // set Trans* = 'T' for *^T
        bool multiply(matrix& B, matrix& output, char TRANSB='N');
        
        /** TODO USE_MAC **/
        
        // create this tensor with "B" and put it in output
        void tensor(matrix* B, matrix* output)
        {
            if(!B || !output)
                return;
            
            output->create(getRows()*B->getRows(), getCols()*B->getCols());
            
            for(int r=0; r<output->getRows(); ++r)
            {
                for(int c=0; c<output->getCols(); ++c)
                {
                    int aRow = r / B->getRows();
                    int aCol = c / B->getCols();
                    int bRow = r % B->getRows();
                    int bCol = c % B->getCols();
                    
                    //output->setDirect(r, c, getDirect(aRow, aCol)* B->getDirect(bRow, bCol));
                    quantum::Complex Temp1 = quantum::Complex(getDirect(aRow, aCol).r, getDirect(aRow, aCol).i);
                    quantum::Complex Temp2 = quantum::Complex(B->getDirect(bRow, bCol).r, B->getDirect(bRow, bCol).i);
                    Temp1 = Temp1.multiply(Temp2);
                    __CLPK_doublecomplex res;
                    res.r = Temp1.getReal();
                    res.i = Temp1.getImaginary();
                    output->setDirect(r, c, res);
                }
            }
        }
        
        // Von Nueman Entropy
        double entropy()
        {
            std::vector<double> ev;
            double sum = 0;
            
            eigenValues(ev);
            
            for(int i=0; i < ev.size(); ++i)
            {
                sum = sum - ev[i] * safeLog(ev[i]);
            }
            return sum;
        }
        
        int getRows()
        {
            return r;
        }
        int getCols()
        {
            return c;
        }
        int getSize()
        {
            return n;
        }
        
        // Matrix Additions 8/21**************************************************************************
#ifdef USE_MAC
        __CLPK_doublecomplex dotProduct(matrix& other)
#else
        Complex dotProduct(matrix& other)
#endif
        {
#ifdef USE_MAC
            __CLPK_doublecomplex sum;
            sum.r = 0;
            sum.i = 0;
#else
            Complex sum = 0 + 0*I;
#endif
            
            for(int i=0; i<getRows(); ++i)
            {
#ifdef USE_MAC
                sum.r += (*this)(i,0).r*other(i,0).r + (*this)(i,0).i*other(i,0).i;
                sum.i += (*this)(i,0).i*other(i,0).r - (*this)(i,0).r*other(i,0).i;
#else
                Complex c;
                c = creal(other(i,0)) - cimag(other(i,0))*I;
                sum = sum + (*this)(i,0)*c;
#endif
            }
            
            return sum;
        }
        
        double norm()
        {
#ifdef USE_MAC
            __CLPK_doublecomplex n = dotProduct(*this);
            return sqrt(n.r);
#else
            Complex n = dotProduct(*this);
            return sqrt(creal(n));
#endif
        }
        
        static void orthogonalize(std::vector < matrix >& B, std::vector <matrix>& output)
        {
            if(B.empty()) return;
            
            std::vector < matrix > u;
            u.push_back(B[0]);
            for(int i=1; i<(int)B.size(); ++i)
            {
                u.push_back(B[i]);
                
#ifdef USE_MAC
                __CLPK_doublecomplex top, bottom, p, c;
#else
                Complex top, bottom, p, c;
#endif
                for(int k=0; k<i; ++k)
                {
                    top = B[i].dotProduct(u[k]);
                    p = u[k].dotProduct(u[k]);
                    
                    // now compute bottom = 1/p:
                    // and then p = top*bottom
#ifdef USE_MAC
                    bottom.r = p.r / (p.r*p.r + p.i*p.i);
                    bottom.i = -p.i /(p.r*p.r + p.i*p.i);
                    
                    p.r = top.r*bottom.r - top.i*bottom.i;
                    p.i = top.i*bottom.r + top.r*bottom.i;
#else
                    bottom = (creal(p) - I*cimag(p)) / (creal(p)*creal(p) + cimag(p)*cimag(p));
                    
                    p = top*bottom;
#endif
                    
                    for(int j=0; j<(int)B[i].getRows(); ++j)
                    {
                        // c = p*u[k](j,0)
                        // then u[i](j,0) -= c
#ifdef USE_MAC
                        c.r = p.r*u[k](j,0).r - p.i*u[k](j,0).i;
                        c.i = p.i*u[k](j,0).r + p.r*u[k](j,0).i;
                        
                        u[i](j,0).r = u[i](j,0).r - c.r;
                        u[i](j,0).i = u[i](j,0).i - c.i;
#else
                        c = p*u[k](j,0);
                        u[i](j,0) = u[i](j,0) - c;
#endif
                    }
                }
            }
            
            output.resize(B.size());
            for(int i=0; i<(int)B.size(); ++i){
                double n = u[i].norm();
                if(n > 0.0000001)
                    n = 1.0/n;
                
                output[i] = u[i];
                for(int j=0; j<B[i].getRows(); ++j){
#ifdef USE_MAC
                    output[i](j,0).r *= n;
                    output[i](j,0).i *= n;
#else
                    output[i](j,0) *= n;
#endif
                }
            }
        }
        
        static double safeLog(double x)
        {
            if(x < 0.0000001)
                return 0;
            else
                return log(x)/log(2.0);
        }
        // End Matrix Additions 8/21**********************************************************************
        
    private:
        int alg_min(int a, int b)
        {
            return (a < b) ? a : b;
        }
        int alg_max(int a, int b)
        {
            return (a > b) ? a : b;
        }
        
#ifdef USE_LINUX
        double* elements;
#elif defined(USE_MAC)
        __CLPK_doublecomplex* elements;
#else
        Complex* elements;
#endif
        int r, c, n;
        
        bool noCleanup;
        
        std::string varName;
    };
    
    typedef matrix mat;
};

#endif
