////// This file has been modified by Alberto Chiappa, matr. 8854468 \\\\\\\

// Here is the solution to challenge 1.3. The tridiagonal matrix representing 
// the finite elements problem has been built directly and stored inside three vectors.
// Thomas algorithm has been implemented to solve the system directly.




#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
    
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  
  //const int&    itermax= param.itermax;   //max number of iteration for Gauss-Siedel
  //const double& toler=param.toler;   // Tolerance for stopping criterion
  //const int& useL2norm=param.useL2norm; // Norm used to compare increment between iterations
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto&    M=param.M; // Number of grid elements
  const auto& outname=param.outname; // Name of the output
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1);
     
  // Thomas algorithm for tridiagonal matrices
  
  // Explicitly creating the matrix: as it is tridiagonal it is
  // possible to store it in 3 vectors, in order to save mamory
  // and to have an efficient code
  
  std::vector<double> a(M,2.+h*h*act),b(M-1,-1.),c(M-1,-1.);
  a[M-1]=1.;
  // LU decomposing the matrix. The result will be stored in 
  // the vectors of the matrix (to save memory)
  for(int m=0; m<M-1; ++m){
  	c[m]=c[m]/a[m];
  	a[m+1]-=b[m]*c[m];
  		}
  // Solving the system
  std::vector <double> y(M);
  y[0]=To-Te;
  for(int m=0;m<M-1;++m){
  	y[m+1]=-c[m]*y[m];
  	}
  theta[0]=To-Te;
  theta[M]=y[M-1]/a[M-1];
  for(int m=M-1;m>0;--m){
  	theta[m]=(y[m-1]-b[m-1]*theta[m+1])/a[m-1];
  		}
  
  	
  /*
  // Gauss-Seidel
  // epsilon=||x^{k+1}-x^{k}||
  // Stopping criteria epsilon<=toler
  
  int iter=0;
  double xnew, epsilon;
     do
       { epsilon=0.;

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   epsilon += (xnew-theta[m])*(xnew-theta[m]);
	   theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1]; 
	 epsilon += (xnew-theta[M])*(xnew-theta[M]);
	 theta[M]=  xnew; 

	 iter=iter+1;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );

    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }
 */
 /*if(useL2norm){
   // Gauss-Seidel
  // epsilon=sqrt(int((x^{k+1}-x^{k})^2) (a first order approximation of L2 norm has been used)
  // Stopping criteria epsilon<=toler
  
  int iter=0;
  double xnew, epsilon;
     do
       { epsilon=0.;

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   
	   //internal nodes have weight h
	   epsilon += h*(xnew-theta[m])*(xnew-theta[m]);
	   theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1]; 
	 
	 //last node has weight h/2
	 epsilon += h/2*(xnew-theta[M])*(xnew-theta[M]);
	 theta[M]=  xnew; 

	 iter=iter+1;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );

    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }
 } else {
 
 // Gauss-Seidel
  // epsilon=sqrt(int((x^{k+1}-x^{k})^2+int((x'^{k+1}-x'^{k})^2) (a first order approximation of H1 norm has been used)
  // Stopping criteria epsilon<=toler
  
  int iter=0;
  double epsilon;
  auto theta_prov=theta;
     do
       { epsilon=0.;

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   theta_prov[m]  = (theta_prov[m-1]+theta_prov[m+1])/(2.+h*h*act);
	   
	   //internal nodes have weight h
	   epsilon += h*(theta_prov[m]-theta[m])*(theta_prov[m]-theta[m]);
         }

	 //Last row
	 theta_prov[M] = theta[M-1]; 
	//last node has weight h/2
	 epsilon += h/2*(theta_prov[M]-theta[M])*(theta_prov[M]-theta[M]);
	  
	 
	 // computation of gradient
	 std::vector<double> der(M);
	 for(int i=0; i<M; ++i){
	 	der[i] = ((theta_prov[i+1]-theta_prov[i])-(theta[i+1]-theta[i]))/h;
	 	}
	 // using rectangle quadrature formula
	          
	 // all medium points have weight h
	 for(int m=0;m<M;++m){
	 	epsilon+=h*der[m]*der[m];
	 	}
	 theta=theta_prov;
	 iter=iter+1;
       }while((sqrt(epsilon) > toler) && (iter < itermax) );

    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }
 }*/

 // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 

     Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);

     cout<<"Result file: "<<outname<<endl;
     //const char* outname="result.dat";
     ofstream f(outname);
     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         //f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<endl;
         f<<m*h*L<<"\t"<<theta[m]+Te<<"\t"<<thetaa[m]<<endl;
	 // An example of use of tie and tuples!
         
	 std::tie(coor[m],sol[m],exact[m])=
	   //std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);
	   std::make_tuple(m*h*L,theta[m]+Te,thetaa[m]);
       }
     // Using temporary files (another nice use of tie)
     gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<std::endl;
     f.close();
     return status;
}
