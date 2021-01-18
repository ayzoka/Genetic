#pragma once
#include <vector>
#include "lib-quantum/Algebra.h"
namespace QKD
{
struct Basis
{
  std::vector < std::vector < quantum::Complex > > B;
  void orthogonalize();

  quantum::Complex innerProduct(std::vector < quantum::Complex >& a, std::vector < quantum::Complex >& b)
  {
    quantum::Complex sum(0,0);

    for(int i=0; i<a.size(); ++i)
      sum = sum + a[i] * b[i].conj();

    return sum;
  }

  double norm(std::vector <quantum::Complex> & a)
  {
    double sum = 0;
    quantum::Complex b;
    for(int i=0; i<a.size(); ++i){
      b = a[i] * a[i].conj();
      sum += b.getReal();
    }

    return sqrt(sum);
  }
};

struct Stats
{
  double p00, p01, p10, p11;
  double Qx, Qy; // Qx = 1-paa; Qy = 1-pbb
  double p0a, p1a, pa0;
  double p0b, p1b, pb0;

  double BB84Rate();
  double SafeLog(double x)
  {
    if(x < 0.0000001)
      return 0;
    return log(x)/log(2.0);
  }
  double entropy(double x)
  {
    if(x > 1.00001){
      std::cout << "Entropy error: x = " << x << "\n";
      return 0;
      //throw std::string("Entropy Error.");
    }

    return -x*SafeLog(x) - (1-x)*SafeLog(1-x);
  }
  double HAB(double key00, double key01, double key10, double key11)
  {
    double first = -key00*SafeLog(key00) - key01*SafeLog(key01) - key10*SafeLog(key10) - key11*SafeLog(key11);

    if(key00+key10 > 1.00001)
      std::cout << "HAB ERROR!\n";
    double second = entropy(key00+key10);

    return first-second;
  }


  void cnot();
  void symmetric(double Q);
  void printChannel(std::ostream& f)
  {
    f << "CHANNEL STATS:\n";
    f << "p00 = " << p00 << "\n";
    f << "p01 = " << p01 << "\n";
    f << "p10 = " << p10 << "\n";
    f << "p11 = " << p11 << "\n\n";

    f << "Qx = " << Qx << "\n";
    f << "Qy = " << Qy << "\n\n";

    f << "p0+ = " << p0a << "\n";
    f << "p1+ = " << p1a << "\n";
    f << "p+0 = " << pa0 << "\n\n";

    f << "p0Y = " << p0b << "\n";
    f << "p1Y = " << p1b << "\n";
    f << "pY0 = " << pb0 << "\n\n";
  }

  double randFloat(double small, double large)
  {
    int r = rand()%RAND_MAX;
    double r2 = r / (double)RAND_MAX;
    
    return (large-small)*r2 + small;
  }

  void generateRandom(double p01, double p10);
  void chooseRandomVector(double minScale, double maxScale, std::vector <quantum::Complex>& v, int startIndex, int endIndex);
};

class Attack
{
  static constexpr double precision = 0.01;

public:
  Attack(Stats& _stats) : stats(_stats)
  {
    setLocalVars();

    v0.create(4,1);
    v1.create(4,1);
    v2.create(4,1);
    v3.create(4,1);

    v0.zero();
    v1.zero();
    v2.zero();
    v3.zero();

#ifdef USE_LINUX
    v0(0,0) = 1;
    v1(1,0) = 1;
    v2(2,0) = 1;
    v3(3,0) = 1;
#else
    v0(0,0).r = 1;
    v1(1,0).r = 1;
    v2(2,0).r = 1;
    v3(3,0).r = 1;
#endif
  }

  // sets E03 and E12 to their start values
  void begin()
  {
    double im03lb = -SafeSqrt(E00*E33 - ReE03*ReE03);
    double im12lb = -SafeSqrt(E11*E22 - ReE12*ReE12);

    E03.setImag(im03lb);
    E12.setImag(im12lb);

    debugStart = true;
  }

  // gets the next possible attack; returns TRUE when this is the last legal attack; FALSE o/w
  bool getNext(algebra::mat & e0, algebra::mat& e1, algebra::mat& e2, algebra::mat& e3)
  {
    bool done = false;
    do{
      //E12 = quantum::Complex(-.0361068, .00601651);
      //E03 = quantum::Complex(.73823, -.230657);

      done = buildAttackVec(e0, e1, e2, e3);

      // add to E12.imag:
      E12.setImag(E12.getImaginary()+precision);
      if(E12.sq() > E11*E22){
	double im12lb = -SafeSqrt(E11*E22 - ReE12*ReE12);
	E12.setImag(im12lb);
	
	E03.setImag(E03.getImaginary()+precision);
	if(E03.sq() > E00*E33){
	  if(debugStart)
	    throw(std::string("NO ATTACKS"));
	  return true; // no more attacks to consider
	}
      }
    }while(!done);

    debugStart= false;
    return false;
  }

    bool debugStart;

  double SafeSqrt(double x)
  {
    if(x < -0.00001){
      std::cout << "\n\nSafeSQRT Error: x = " << x << "\n\n";
      throw(std::string("SafeSQRT Error; X negative."));
    }
    if(x < .00000001)
      return 0;

    return sqrt(x);
  }

  bool buildAttackVec(algebra::mat & e0, algebra::mat& e1, algebra::mat& e2, algebra::mat& e3)
  {
    // e0:
    quantum::Complex a0 = SafeSqrt(stats.p00);

    // e3:
    quantum::Complex a3 = quantum::Complex(E03.getReal()/a0.getReal(), E03.getImaginary() / a0.getReal());


    if(stats.p11 - a3.getReal()*a3.getReal() - a3.getImaginary()*a3.getImaginary() < 0)
      return false;
    quantum::Complex b3 = SafeSqrt(stats.p11 - a3.getReal()*a3.getReal() - a3.getImaginary()*a3.getImaginary());


    // e1:
    quantum::Complex a1 = quantum::Complex(E01.getReal()/a0.getReal(), E01.getImaginary() / a0.getReal());
    quantum::Complex b1 = E13 - a1.conj()*a3;
    b1 = b1 * (1.0/b3.getReal());
    b1 = b1.conj();
    if(b3.getReal() < .0000001){
      b1 = 0;
    }

    if(stats.p01 - a1.sq() - b1.sq() < 0)
      return false;

    quantum::Complex c1 = SafeSqrt(stats.p01 - a1.sq() - b1.sq());


    // e2:
    quantum::Complex a2 = quantum::Complex(E02.getReal()/a0.getReal(), E02.getImaginary() / a0.getReal());
    quantum::Complex b2 = E23 - a2.conj()*a3;
    b2 = b2 * (1.0/b3.getReal());
    b2 = b2.conj();

    if(b3.getReal() < .0000001){
      // check to make sure this is a valid solution:
      quantum::Complex check = a1.conj()*a3;
      if(fabs(check.getReal() - E13.getReal()) > .0000001 ||
	 fabs(check.getImag() - E13.getImag()) > .0000001){
	//std::cout << "INVALID\n";
	return false;
      }
      b2 = 0;
    }

    quantum::Complex c2 = E12 - a1.conj()*a2 - b1.conj()*b2;


    c2 = c2 * (1.0/c1.getReal());

    if(c1.getReal() < .0000001){
      // check:
      quantum::Complex check = a1.conj()*a2 + b1.conj()*b2;
      if(fabs(check.getReal() - E12.getReal()) > .0000001 ||
	 fabs(check.getImag() - E12.getImag()) > .0000001){
	//std::cout << "INVALID12\n";
	return false;
      }
      c2 = 0;
    }

    

    double test = stats.p10 - a2.sq() - b2.sq() - c2.sq();
    if(test < 0)
      return false;
    quantum::Complex d2 = SafeSqrt(stats.p10 - a2.sq() - b2.sq() - c2.sq());


    /*
    // a test:
    std::cout << "\n--test--\n";
    // e0e2:
    quantum::Complex t = a0.conj()*a2;
    std::cout << "E02: ";
    t.print();
    //e1e3:
    t = a1.conj()*a3 + b1.conj()*b3;
    std::cout << "\nE13: ";
    t.print();
    std::cout << "--end test--\n\n";
    */

    e0.zero();
    e1.zero();
    e2.zero();
    e3.zero();

    e0.add(a0, v0);

    e3.add(a3, v0);
    e3.add(b3, v1);

    e1.add(a1, v0);
    e1.add(b1, v1);
    e1.add(c1, v2);

    e2.add(a2, v0);
    e2.add(b2, v1);
    e2.add(c2, v2);
    e2.add(d2, v3);

    return true;
  }

 private:
  void setLocalVars()
  {
    E00 = stats.p00;
    E11 = stats.p01;
    E22 = stats.p10;
    E33 = stats.p11;

    E01 = quantum::Complex(stats.p0a-.5, stats.p0b-.5);
    E23 = quantum::Complex(stats.p1a-.5, stats.p1b-.5);
    E02 = quantum::Complex(stats.pa0-.5*(E00+E22), .5*(E00+E22)-stats.pb0);
    E13 = E02*(-1);

    ReE03 = 1 - stats.Qx - stats.Qy - .5*(E01.getReal() + E01.getImaginary() + E23.getReal() + E23.getImaginary());
    ReE12 = stats.Qy - stats.Qx + .5*(E01.getImaginary() - E01.getReal() + E23.getImaginary() - E23.getReal());

    E03 = quantum::Complex(ReE03, 0);
    E12 = quantum::Complex(ReE12, 0);

    /*std::cout << "FOUND:\n";
    std::cout << "E00:" << E00 << "\n";

    std::cout << "E11:" << E11 << "\n";

    std::cout << "E22:" << E22 << "\n";

    std::cout << "E33:" << E33 << "\n";


    std::cout << "E01:\n";
    E01.print();
    std::cout << "\nE23:\n";
    E23.print();
    std::cout << "\nE02:\n";
    E02.print();
    std::cout << "\nE13:\n";
    E13.print();
    std::cout << "\n\n\n";*/
  }

  Stats stats;
  algebra::mat v0, v1, v2, v3; // standard basis vectors

  double E00, E11, E22, E33;
  quantum::Complex E01, E23, E02, E13;

  double ReE03, ReE12;
  quantum::Complex E03, E12;
};
}
