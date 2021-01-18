/*
Complex.cpp
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "Complex.hpp"

namespace quantum
{
    Complex::Complex()
    {
        this->real = 0.0;
        this->imaginary = 0.0;
    }
    
    Complex::Complex(double realIn)
    {
        this->real = realIn;
        this->imaginary = 0.0;
    }

    
    Complex::Complex(double realIn, double imagIn)
    {
        this->real = realIn;
        this->imaginary = imagIn;
    }
    
    Complex Complex::add(const Complex &complexIn)
    {
        Complex C;
        C.real = this->real + complexIn.real;
        C.imaginary = this->imaginary + complexIn.imaginary;
        return C;
    }

    Complex Complex::subtract(const Complex &complexIn)
    {
        Complex C;
        C.real = this->real - complexIn.real;
        C.imaginary = this->imaginary - complexIn.imaginary;
        return C;
    }
    
    Complex Complex::multiply(const Complex &complexIn)
    {
        Complex C;
        C.real = this->real * complexIn.real - this->imaginary * complexIn.imaginary;
        C.imaginary = this->imaginary * complexIn.real + this->real * complexIn.imaginary;
        return C;
    }
    
    Complex Complex::divide(const Complex &complexIn)
    {
        Complex C;
        C.real = (this->real * complexIn.real + this->imaginary * complexIn.imaginary)/(complexIn.real * complexIn.real + complexIn.imaginary * complexIn.imaginary);
        C.imaginary = (this->imaginary * complexIn.real - this->real * complexIn.imaginary)/(complexIn.real * complexIn.real + complexIn.imaginary * complexIn.imaginary);
        return C;
    }
    
    Complex Complex::square()
    {
        Complex C;
        C.real = this->real;
        C.imaginary = this->imaginary;
        return C.multiply(C);
    }
    
    Complex Complex::pow(double exponent)
    {
        Complex C;
        Complex temp;
        if(exponent == 0)
        {
            C.real = 1.0;
            C.imaginary = 0.0;
            return C;
        }
        else if(exponent == 1)
        {
            C.real = this->real;
            C.imaginary = this->imaginary;
            return C;
        }
        else
        {
            /*
             (a+b*i)^n =
             Real = (a^2+b^2)^(n/2)*cos(n*arctan(b/a))
             Imaginary = (a^2+b^2)^(n/2)*sin(n*arctan(b/a))
             */
            
            C.real = (std::pow( (std::pow(this->real,2) + std::pow(this->imaginary,2)),(exponent/2.0)) * cos(exponent * atan(this->imaginary/this->real)) );
            C.imaginary = (std::pow( (std::pow(this->real,2) + std::pow(this->imaginary,2)),(exponent/2.0)) * sin(exponent * atan(this->imaginary/this->real)) );
            return C;
        }
    }
    
    Complex Complex::conjugate()
    {
        Complex C;
        C.real = this->real;
        C.imaginary = -(this->imaginary);
        return C;
    }
    
    double Complex::getReal()
    {
        return this->real;
    }
    
    double Complex::getImaginary()
    {
        return this->imaginary;
    }
    
    void Complex::print()
    {
        if(this->imaginary > 0.0)
        {
            if(this->imaginary != 1.0)
            {
                if(this->real != 0)
                {
                    std::cout << this->real << "+" << this->imaginary << "i";
                }
                else
                {
                    std::cout << this->imaginary << "i";
                }
            }
            else
            {
                if(this->real != 0.0)
                {
                    std::cout << this->real << "+i";
                }
                else
                {
                    std::cout << "i";
                }
            }
        }
        else if(this->imaginary < 0.0)
        {
            if(this->imaginary != -1.0)
            {
                if(this->real != 0.0)
                {
                    std::cout << this->real << this->imaginary << "i";
                }
                else
                {
                    std::cout << this->imaginary << "i";
                }
            }
            else
            {
                std::cout << this->real << "-i";
            }
        }
        else
        {
            std::cout << this->real;
        }
    }
}
