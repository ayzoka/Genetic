/*
Complex.hpp
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef Complex_hpp
#define Complex_hpp

#include <iostream>
#include <cmath>

namespace quantum
{
    class Complex
    {
    public:
        Complex();
        // Default Constructor
        // Sets real and imaginary parts to zero
        
        Complex(double realIn);
        // For Defining A Real Number
        
        Complex(double realIn, double imagIn);
        // For Defining A Complex Number

        Complex add(const Complex &complexIn);
        // Adds paramatized user input to object
        // Returns a complex number
        // Does not change the state of the object
        
        Complex subtract(const Complex &complexIn);
        // Subtracts the paramaterized user input from object
        // Returns a complex number
        // Does not change the state of the object
        
        Complex multiply(const Complex &complexIn);
        // Multiplies object with the paramaterized user input
        // Returns a complex number
        // Does not change the state of the object
        
        Complex divide(const Complex &complexIn);
        // Divides this by the paramaterized user input
        // Returns a complex number
        // Does not change the state of the object

        Complex square();
        // Squares object
        // Returns a complex number
        // Does not change the state of the object

      // return |this|^2
      double sq()
      {
	return real*real + imaginary*imaginary;
      }

      // returns 1/this
      Complex invert()
      {
	double d = getReal()*getReal() + getImaginary()*getImaginary();

	quantum::Complex out(getReal()/d, -getImaginary()/d);
	return out;
      }

      Complex operator*(const Complex& complexIn) const
      {
        Complex C;
        C.real = this->real * complexIn.real - this->imaginary * complexIn.imaginary;
        C.imaginary = this->imaginary * complexIn.real + this->real * complexIn.imaginary;
        return C;
      }
      Complex operator+ (const Complex& other) const
      {
	Complex C = Complex(real + other.real, imaginary + other.imaginary);
	return C;
      }
      Complex operator- (const Complex& other) const
      {
	Complex C = Complex(real - other.real, imaginary - other.imaginary);
	return C;
      }

      Complex conj()
      {
	Complex C = Complex(real, -imaginary);
	return C;
      }

      Complex operator* (double other) const
      {
	Complex C = Complex(real*other, imaginary*other);
	return C;
      }

        
        Complex pow(double exponent);
        // Returns a complex number that equals object^exponent
        // Does not change the state of the object
        
        Complex conjugate();
        // Flips sign of imaginary side
        // Returns a complex number
        // Does not change the state of the object
        
        double getReal();
        
        double getImaginary();
      double getImag(){return getImaginary();};

      void setReal(double r){real = r;}
      void setImag(double im){imaginary = im;}

        void print();
        // Prints this in complex number format i.e. 3+2i
        
    private:
        double real;
        double imaginary;
    };
}

#endif /* Complex_hpp */
