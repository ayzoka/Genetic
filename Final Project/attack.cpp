#include "attack.h"

namespace QKD
{
  void Basis::orthogonalize()
  {
    if(B.empty()) return;

    std::vector < std::vector < quantum::Complex > > u;
    u.push_back(B[0]);

    for(int i=1; i<B.size(); ++i){
      u.push_back(B[i]);
      quantum::Complex top, bottom, p, c;

      for(int k=0; k<i; ++k){
	top = innerProduct(B[i], u[k]);
	p = innerProduct(u[k], u[k]);

	bottom = p.invert();
	p = top*bottom;
	for(int j=0; j<B[i].size(); ++j){
	  c = p*u[k][j];
	  u[i][j] = u[i][j] - c;
	}
      }
    }

    for(int i=0; i<B.size(); ++i){
      double n = norm(u[i]);
      if(n > 0.00000001)
	n = 1.0/n;
      B[i] = u[i];
      for(int j=0; j<B[i].size(); ++j)
	B[i][j] = B[i][j] * n;
    }
  }

  void Stats::chooseRandomVector(double minScale, double maxScale, std::vector <quantum::Complex>& v, int startIndex, int endIndex)
  {
    double norm = 0;
    for(int i=startIndex; i<endIndex; ++i){
      v[i] = quantum::Complex(randFloat(-1,1), randFloat(-1,1));

      norm += v[i].getReal()*v[i].getReal() + v[i].getImaginary()*v[i].getImaginary();
    }

    if(norm < 0.0000001)
      throw(std::string("Choose Random Vector - an unlikely thing has happened!"));

    double scale = randFloat(minScale, maxScale);
    norm = sqrt(scale/norm);
    for(int i=startIndex; i<endIndex; ++i)
      v[i] = v[i]*norm;
  }


  double Stats::BB84Rate()
  {
    double R01 = p0a - .5;
    double R23 = p1a - .5;
    double R02 = pa0 - .5*(1-p01+p10);
    double R13 = -R02;

    double I01 = p0b - .5;
    double I23 = p1b - .5;
    double I02 = .5*(1-p01+p10) - pb0;
    double I13 = -I02;

    double R03 = 1 - Qx - Qy - .5*(R01+R23+I01+I23);
    double R12 = Qy - Qx + .5*(I01+I23-R01-R23);

    double l01 = .5 + sqrt( (p10-p01)*(p10-p01) + 4*R03*R03 ) / (2.0*(2-p01-p10));
    double l23 = .5 + sqrt( (p01-p10)*(p01-p10) + 4*R12*R12 ) / (2.0*(p01+p10));

    double first = entropy( (1-p01) / (2.0-p01-p10) ) - entropy(l01);
    double second = entropy( p01 / (p01+p10) ) - entropy(l23);

    if(p01+p10 <= .0000001)
      second = 0;
    if(2.0-p01-p10 <= 0.0000001)
      first = 0;

    double SAE = (2-p01-p10)/2.0 * first + (p01+p10)/2.0 * second;

    double key00 = .5*(1-p01);
    double key11 = .5*(1-p10);
    double key01 = .5*p01;
    double key10 = .5*p10;

    double H = HAB(key00,key01,key10,key11);

    return SAE - H;
  }

  void Stats::symmetric(double Q)
  {
    p00 = 1-Q;
    p01 = Q;
    p10 = Q;
    p11 = 1-Q;
    
    Qx = Q;
    Qy = Q;
    
    p0a = .5;
    p1a = .5;
    pa0 = .5;
    
    p0b = .5;
    p1b = .5;
    pb0 = .5;
  }

  void Stats::generateRandom(double Q01, double Q10)
  {
    Basis B;
    B.B.resize(2);

    for(int i=0; i<B.B.size(); ++i)
      B.B[i].resize(8);

    chooseRandomVector(1-Q01, 1, B.B[0], 0, 4);
    chooseRandomVector(0, Q01, B.B[0], 4, 8);
    chooseRandomVector(0, Q10, B.B[1], 0, 4);
    chooseRandomVector(0, Q10, B.B[1], 4, 8);

    // e3 = e0 + z:
    for(int i=4; i<8; ++i)
      B.B[1][i] = B.B[1][i] + B.B[0][i-4];

    B.orthogonalize();

    Basis E;
    E.B.resize(4);
    for(int i=0; i<4; ++i)
      E.B[i].resize(4);

    // pack solution:
    for(int i=0; i<4; ++i){
      E.B[0][i] = B.B[0][i];
      E.B[1][i] = B.B[0][4+i];
      E.B[2][i] = B.B[1][i];
      E.B[3][i] = B.B[1][4+i];
    }
  
    // now set stats:
    p01 = E.innerProduct(E.B[1], E.B[1]).getReal();
    p10 = E.innerProduct(E.B[2], E.B[2]).getReal();

    p00 = 1-p01;
    p11 = 1-p10;

    p0a = .5 + E.innerProduct(E.B[0], E.B[1]).getReal();
    p1a = .5 + E.innerProduct(E.B[2], E.B[3]).getReal();
    pa0 = .5*(1-p01) + .5*p10 + E.innerProduct(E.B[0], E.B[2]).getReal();

    p0b = .5 + E.innerProduct(E.B[0], E.B[1]).getImag();
    p1b = .5 + E.innerProduct(E.B[2], E.B[3]).getImag();
    pb0 = .5*(1-p01) + .5*p10 - E.innerProduct(E.B[0], E.B[2]).getImag();

    double R01 = E.innerProduct(E.B[0], E.B[1]).getReal();
    double R03 = E.innerProduct(E.B[0], E.B[3]).getReal();
    double R12 = E.innerProduct(E.B[1], E.B[2]).getReal();
    double R23 = E.innerProduct(E.B[2], E.B[3]).getReal();

    double I01 = E.innerProduct(E.B[0], E.B[1]).getImag();
    double I23 = E.innerProduct(E.B[2], E.B[3]).getImag();

    std::cout << "E12:\n";
    E.innerProduct(E.B[1], E.B[2]).print();
    std::cout << "E03:\n";
    E.innerProduct(E.B[0], E.B[3]).print();
    std::cout << "\n\n";
    Qx = .5 - .5*(R01+R03+R12+R23);
    Qy = .5 - .5*(I01+R03-R12+I23);

    std::cout << "STATS:\n\n";
    std::cout << "\nE01:\n";
    E.innerProduct(E.B[0], E.B[1]).print();
    std::cout << "\nE23:\n";
    E.innerProduct(E.B[2], E.B[3]).print();
    std::cout << "\nE02:\n";
    E.innerProduct(E.B[0], E.B[2]).print();
    std::cout << "\nE13:\n";
    E.innerProduct(E.B[1], E.B[3]).print();
  }
}
