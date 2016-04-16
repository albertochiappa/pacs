#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <iosfwd>
#include <string>
struct parameters
{
  //! max number of iteration for Gauss-Siedel
  int   itermax;
  //! Tolerance for stopping criterion
  double  toler;
  //! Select the norm used to evaluate the increment 
  //! for the stopping criterion: 1 for L2 norm, 0 
  //! for H1 norm
  int useL2norm;
  //! Bar length
  double L;
  //! First longitudinal dimension
  double a1;
 //! Second longitudinal dimension
  double a2;
  //! Dirichlet condition
  double To;
  //! External temperature 
  double Te;
  //! Conductivity
  double k;
  //! Convection coefficient
  double hc;
  //! Number of elements
  int M;
  //! Name of the output file
  std::string outname;
  //! Where the result are stored
  int where_result;
  //! Constructor takes default values
  parameters():
    itermax(1000000),
    toler(1e-8),
    useL2norm(1),
    L(40.),
    a1(4.),
    a2(50.),
    To(46.),
    Te(20.),
    k(0.164),
    hc(1.e-6*200.),
    M(100),
    outname("result.dat"),
    where_result(2)
  {}
};
//! Prints parameters
std::ostream & operator << (std::ostream &,const parameters &);
#endif
