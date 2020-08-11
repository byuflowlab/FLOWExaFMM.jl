/***
 *  For the complex step derivative method:
 *  f'(x) ~  Im [ f(x+ih) ] / h
 *  Define a double complex class that inherits from the
 *  library complex type and overloads appropriate operators.
 *  Mon Jan  8 22:42:20 PST 2001
 *      CODE DOWNLOADED FROM PROF. MARTIN'S WEBSITE
 *          http://mdolab.engin.umich.edu/content/complex-step-guide-cc
 ***/

#ifndef HDRcomplexify
#define HDRcomplexify

#include <complex>
using namespace std;


#if EXAFMM_SINGLE
typedef float real_t;                         //!< Floating point type is single precision
#else
typedef double real_t;                        //!< Floating point type is double precision
#endif


#ifndef HDRderivify
inline real_t real(const real_t& r) {
  /***
   *  So the real() statement can be used even with
   *  the double version of the code to be complexified.
   *  Most useful inside printf statements.
   ***/
return r;
}

inline real_t imag(const real_t& r) {
  return 0.;
}
#endif // HDRderivify


class cplx : public complex<real_t> {
public:
  cplx() : complex<real_t>() {};
  cplx(const real_t& d) : complex<real_t>(d) {};
  cplx(const real_t& r, const real_t& i) : complex<real_t>(r,i) {};
  cplx(const complex<double>& z) : complex<real_t>(z) {};
  cplx(const complex<float>& z) : complex<real_t>(z) {};
  operator real_t() {return this->real();}
  operator int() {return int(this->real());}
  // relational operators
  // Conversion constructor should be able to take care of the
  // operator== and != calls with double, but MIPS compiler
  // complains of ambiguous inheritance.  This should be more
  // efficient anyway.  (A hint of what comes below.)
  friend inline bool operator==(const cplx&,const cplx&);
  friend inline bool operator==(const cplx&,const real_t&);
  friend inline bool operator==(const real_t&,const cplx&);
  friend inline bool operator!=(const cplx&,const cplx&);
  friend inline bool operator!=(const cplx&,const real_t&);
  friend inline bool operator!=(const real_t&,const cplx&);
  friend inline bool operator>(const cplx&,const cplx&);
  friend inline bool operator>(const cplx&,const real_t&);
  friend inline bool operator>(const real_t&,const cplx&);
  friend inline bool operator<(const cplx&,const cplx&);
  friend inline bool operator<(const cplx&,const real_t&);
  friend inline bool operator<(const real_t&,const cplx&);
  friend inline bool operator>=(const cplx&,const cplx&);
  friend inline bool operator>=(const cplx&,const real_t&);
  friend inline bool operator>=(const real_t&,const cplx&);
  friend inline bool operator<=(const cplx&,const cplx&);
  friend inline bool operator<=(const cplx&,const real_t&);
  friend inline bool operator<=(const real_t&,const cplx&);
  // here's the annoying thing:
  // Every function in class complex<double> that returns a
  // complex<double> causes ambiguities with function overloading
  // resolution because of the mix of types cplx and
  // complex<double> and double and int in math expressions.
  // So, although they are inherited, must redefine them
  // to return type cplx:
  // basic arithmetic
  inline cplx operator+() const;
  inline cplx operator+(const cplx&) const;
  inline cplx operator+(const real_t&) const;
  inline cplx operator+(const int&) const;
  inline friend cplx operator+(const real_t&, const cplx&);
  inline friend cplx operator+(const int&, const cplx&);
  inline cplx operator-() const;
  inline cplx operator-(const cplx&) const;
  inline cplx operator-(const real_t&) const;
  inline cplx operator-(const int&) const;
  inline friend cplx operator-(const real_t&, const cplx&);
  inline friend cplx operator-(const int&, const cplx&);
  inline cplx operator*(const cplx&) const;
  inline cplx operator*(const real_t&) const;
  inline cplx operator*(const int&) const;
  inline friend cplx operator*(const real_t&, const cplx&);
  inline friend cplx operator*(const int&, const cplx&);
  inline cplx operator/(const cplx&) const;
  inline cplx operator/(const real_t&) const;
  inline cplx operator/(const int&) const;
  inline friend cplx operator/(const real_t&, const cplx&);
  inline friend cplx operator/(const int&, const cplx&);
  // from <math.h>
  inline friend cplx sin(const cplx&);
  inline friend cplx sinh(const cplx&);
  inline friend cplx cos(const cplx&);
  inline friend cplx cosh(const cplx&);
  inline friend cplx tan(const cplx&);
  inline friend cplx tanh(const cplx&);
  inline friend cplx log10(const cplx&);
  inline friend cplx log(const cplx&);
  inline friend cplx sqrt(const cplx&);
  inline friend cplx exp(const cplx&);
  inline friend cplx pow(const cplx&, const cplx&);
  inline friend cplx pow(const cplx&, const real_t&);
  inline friend cplx pow(const cplx&, const int&);
  inline friend cplx pow(const real_t&, const cplx&);
  inline friend cplx pow(const int&, const cplx&);
  // complex versions of these are not in standard library
  // or they need to be redefined:
  // (frexp, modf, and fmod have not been dealt with)
  inline friend cplx fabs(const cplx&);
  inline friend cplx asin(const cplx&);
  inline friend cplx acos(const cplx&);
  inline friend cplx atan(const cplx&);
  inline friend cplx atan2(const cplx&, const cplx&);
  inline friend cplx ceil(const cplx&);
  inline friend cplx floor(const cplx&);
  inline friend cplx ldexp(const cplx&, const int&);
};

inline bool operator==(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) == real(rhs);
}

inline bool operator==(const cplx& lhs, const real_t& rhs)
{
  return real(lhs) == rhs;
}

inline bool operator==(const real_t& lhs, const cplx& rhs)
{
  return lhs == real(rhs);
}

inline bool operator!=(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) != real(rhs);
}

inline bool operator!=(const cplx& lhs, const real_t& rhs)
{
  return real(lhs) != rhs;
}

inline bool operator!=(const real_t& lhs, const cplx& rhs)
{
  return lhs != real(rhs);
}

inline bool operator>(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) > real(rhs);
}

inline bool operator>(const cplx& lhs, const real_t& rhs)
{
  return real(lhs) > rhs;
}

inline bool operator>(const real_t& lhs, const cplx& rhs)
{
  return lhs > real(rhs);
}

inline bool operator<(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) < real(rhs);
}

inline bool operator<(const cplx& lhs, const real_t& rhs)
{
  return real(lhs) < rhs;
}

inline bool operator<(const real_t& lhs, const cplx& rhs)
{
  return lhs < real(rhs);
}

inline bool operator>=(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) >= real(rhs);
}

inline bool operator>=(const cplx& lhs, const real_t& rhs)
{
  return real(lhs) >= rhs;
}

inline bool operator>=(const real_t& lhs, const cplx& rhs)
{
  return lhs >= real(rhs);
}

inline bool operator<=(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) <= real(rhs);
}

inline bool operator<=(const cplx& lhs, const real_t& rhs)
{
  return real(lhs) <= rhs;
}

inline bool operator<=(const real_t& lhs, const cplx& rhs)
{
  return lhs <= real(rhs);
}

inline cplx cplx::operator+() const
{
  return +complex<real_t>(*this);
}

inline cplx cplx::operator+(const cplx& z) const
{
  return complex<real_t>(*this)+complex<real_t>(z);
}

inline cplx cplx::operator+(const real_t& r) const
{
  return complex<real_t>(*this)+r;
}

inline cplx cplx::operator+(const int& i) const
{
  return complex<real_t>(*this)+real_t(i);
}

inline cplx operator+(const real_t& r, const cplx& z)
{
  return r+complex<real_t>(z);
}

inline cplx operator+(const int& i, const cplx& z)
{
  return real_t(i)+complex<real_t>(z);
}

inline cplx cplx::operator-() const
{
  return -complex<real_t>(*this);
}

inline cplx cplx::operator-(const cplx& z) const
{
  return complex<real_t>(*this)-complex<real_t>(z);
}

inline cplx cplx::operator-(const real_t& r) const
{
  return complex<real_t>(*this)-r;
}

inline cplx cplx::operator-(const int& i) const
{
  return complex<real_t>(*this)-real_t(i);
}

inline cplx operator-(const real_t& r, const cplx& z)
{
  return r-complex<real_t>(z);
}

inline cplx operator-(const int& i, const cplx& z)
{
  return real_t(i)-complex<real_t>(z);
}

inline cplx cplx::operator*(const cplx& z) const
{
  return complex<real_t>(*this)*complex<real_t>(z);
}

inline cplx cplx::operator*(const real_t& r) const
{
  return complex<real_t>(*this)*r;
}

inline cplx cplx::operator*(const int& i) const
{
  return complex<real_t>(*this)*real_t(i);
}

inline cplx operator*(const real_t& r, const cplx& z)
{
  return r*complex<real_t>(z);
}

inline cplx operator*(const int& i, const cplx& z)
{
  return real_t(i)*complex<real_t>(z);
}

inline cplx cplx::operator/(const cplx& z) const
{
  return complex<real_t>(*this)/complex<real_t>(z);
}

inline cplx cplx::operator/(const real_t& r) const
{
  return complex<real_t>(*this)/r;
}

inline cplx cplx::operator/(const int& i) const
{
  return complex<real_t>(*this)/real_t(i);
}

inline cplx operator/(const real_t& r, const cplx& z)
{
  return r/complex<real_t>(z);
}

inline cplx operator/(const int& i, const cplx& z)
{
  return real_t(i)/complex<real_t>(z);
}

inline cplx sin(const cplx& z)
{
  return sin(complex<real_t>(z));
}

inline cplx sinh(const cplx& z)
{
  return sinh(complex<real_t>(z));
}

inline cplx cos(const cplx& z)
{
  return cos(complex<real_t>(z));
}

inline cplx cosh(const cplx& z)
{
  return cosh(complex<real_t>(z));
}

#ifdef __GNUC__ // bug in gcc ?? get segv w/egcs-2.91.66 and 2.95.2
inline cplx tan(const cplx& z)
{
  return sin(complex<real_t>(z))/cos(complex<real_t>(z));
}

inline cplx tanh(const cplx& z)
{
  return sinh(complex<real_t>(z))/cosh(complex<real_t>(z));
}

inline cplx log10(const cplx& z)
{
  return log(complex<real_t>(z))/log(real_t(10.));
}
#else
inline cplx tan(const cplx& z)
{
  return tan(complex<real_t>(z));
}

inline cplx tanh(const cplx& z)
{
  return tanh(complex<real_t>(z));
}

inline cplx log10(const cplx& z)
{
  return log10(complex<real_t>(z));
}
#endif

inline cplx log(const cplx& z)
{
  return log(complex<real_t>(z));
}

inline cplx sqrt(const cplx& z)
{
  return sqrt(complex<real_t>(z));
}

inline cplx exp(const cplx& z)
{
  return exp(complex<real_t>(z));
}

inline cplx pow(const cplx& a, const cplx& b)
{
  return pow(complex<real_t>(a),complex<real_t>(b));
}

inline cplx pow(const cplx& a, const real_t& b)
{
  return pow(complex<real_t>(a),b);
}

inline cplx pow(const cplx& a, const int& b)
{
  return pow(complex<real_t>(a),real_t(b));
}

inline cplx pow(const real_t& a, const cplx& b)
{
  return pow(a,complex<real_t>(b));
}

inline cplx pow(const int& a, const cplx& b)
{
  return pow(real_t(a),complex<real_t>(b));
}

inline cplx fabs(const cplx& z)
{
  return (real(z)<0.0) ? -z:z;
}

#define surr_TEENY (1.e-24) /* machine zero compared to nominal magnitude of
			       the real part */

inline cplx asin(const cplx& z)
{
  // derivative trouble if imag(z) = +/- 1.0
  return cplx(asin(real(z)),imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

inline cplx acos(const cplx& z)
{
  // derivative trouble if imag(z) = +/- 1.0
  return cplx(acos(real(z)),-imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

#undef surr_TEENY

inline cplx atan(const cplx& z)
{
  return cplx(atan(real(z)),imag(z)/(1.0+real(z)*real(z)));
}

inline cplx atan2(const cplx& z1, const cplx& z2)
{
  return cplx(atan2(real(z1),real(z2)),
	      (real(z2)*imag(z1)-real(z1)*imag(z2))
	      /(real(z1)*real(z1)+real(z2)*real(z2)));
}

inline cplx ceil(const cplx& z)
{
  return cplx(ceil(real(z)),0.);
}

inline cplx floor(const cplx& z)
{
  return cplx(floor(real(z)),0.);
}

inline cplx ldexp(const cplx& z, const int& i)
{
  return cplx(ldexp(real(z),i),ldexp(imag(z),i));
}

#endif
