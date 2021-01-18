/*
Algebra.cpp (Mac Version)
Written by Walter O. Krawec and Eric Geiss
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "Algebra.h"

#ifdef USE_LINUX
extern "C"
{
#include <cblas.h>
#include <lapacke.h>
}
#endif


namespace algebra
{
    bool matrix::eigenValues(std::vector <double>& values)
    {

#ifdef USE_LAPACK
        int n, lda, info, lwork, lrwork, liwork;
        n = getRows();
        lda = n;
        int iwkopt;
        int* iwork;
        double rwkopt;
        double* rwork;
#ifdef USE_MAC
        __CLPK_doublecomplex wkopt, *work;
#else
        Complex wkopt, *work;
#endif
        
        double *w = new double[n];
    
        matrix a;
        a.create(n);
        for(int r=0; r<n; ++r)
        {
            for(int c=0; c<n; ++c)
            {
#ifdef USE_MAC
                a(r,c).r = 0;
                a(r,c).i = 0;
#else
                a(r,c) = 0;
#endif
	    
                if(c >= r)
                    a(r,c) = (*this)(r,c);
            }
        }
        //a.print(std::cout);
    
        lwork = -1;
        lrwork = -1;
        liwork = -1;
    
#ifdef USE_MAC
        zheevd_("N", "Upper", &n, a.getElements(), &lda, w, &wkopt, &lwork, &rwkopt, &lrwork, &iwkopt, &liwork, &info);
        
        lwork = (int)wkopt.r;
        
        work = new __CLPK_doublecomplex[lwork];
#else
        LAPACK_zheevd("N", "Upper", &n, a.getElements(), &lda, w, &wkopt, &lwork, &rwkopt, &lrwork, &iwkopt, &liwork, &info);

        lwork = (int)creal(wkopt);
        work = new Complex[lwork];
#endif
        
        
        lrwork = (int)rwkopt;
        rwork = new double[lrwork];
        liwork = iwkopt;
        iwork = new int[liwork];
    
#ifdef USE_MAC
        zheevd_("N", "Upper", &n, a.getElements(), &lda, w, work, &lwork, rwork, &lwork, iwork, &liwork, &info);
#else
        LAPACK_zheevd("N", "Upper", &n, a.getElements(), &lda, w, work, &lwork, rwork, &lwork, iwork, &liwork, &info);
#endif
        if(info == 0)
        {
            values.clear();
            for(int i=0; i<n; ++i)
                values.push_back(w[i]);
        }
    
        delete[] w;
        delete[] work;
        delete[] rwork;
        delete[] iwork;
    
        return (info==0);
#endif
        return false;
    }

    
    // TO DO: change with USE_MAC definition to get this working
    // Test in main method
    
    bool matrix::multiply(matrix& B, matrix& output, char TRANSB)
    {
#ifdef USE_MAC
        if(TRANSB != 'N')
            std::cout << "Warning: matrix transpose not implemented in multiply function; result will be wrong\n";
    
        __CLPK_doublecomplex alpha;
        alpha.r=1.0;
        alpha.i=0.0;
        
        char TRANSA = 'N';
        //char TRANSB = 'N';
        bool printMats = false;
    
        if(!elements || !B.elements)
            return false;
    
        int testThisC = c;
        int testBR = B.getRows();
    
        if(TRANSA == 'T')
            testThisC = r;
        if(TRANSB == 'T')
            testBR = B.getCols();
        if(testThisC != testBR)
        {
            throw(std::string("Multiply sizes wrong"));
            return false;
        }
    
        //std::cout << "Mult stage 0\n";
    
        //std::cout << "Mult stage 1\n";
        //char TRANSA = 'N', TRANSB = 'N';
        __CLPK_doublecomplex beta;
        beta.r = 0.0;
        beta.i = 0.0;
    
    
        if(printMats)
        {
            std::cout << "this:";
            print(std::cout);
            std::cout << "\nB:";
            B.print(std::cout);
        }
    
        int M = r;
        int N = B.c;
        int K = c;
    
        int LDA = r, LDB = c, LDC = r;
    
        if(TRANSA == 'T')
        {
            M = c;
            K = r;
        }
        if(TRANSB == 'T')
        {
            LDB = B.r;
            N = B.r;
        }
    
        if(TRANSA == 'N')
            LDA = M;
        else
            LDA = K;
    
        if(TRANSB == 'N')
            LDB = K;
        else
            LDB = N;
    
        LDC = M;
    
        output.create(M, N);
        if(!output.getElements())
        {
            std::cout << "Mult error!\n";
            throw(1);
        }
    
    
    
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, &alpha, getElements(), LDA, B.getElements(), LDB, &beta, output.getElements(), LDC);

        if(printMats){
            std::cout << "\nC:";
            output.print(std::cout);}

        return true;
#endif
        //std::cout << "Multipy error = " << result << "\n";
        return false;
    }
}
