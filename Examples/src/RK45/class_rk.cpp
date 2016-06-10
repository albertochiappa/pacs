#include <cmath>
#include <algorithm>
#include "rk45.hpp"
template<typename rkType>

class generic_rk: public rkType{ // vedere se funziona private
	
	public:
			std::vector<std::pair<double,double>> solve(
				std::function<double (double const &, double const &)> const & dy,
				double const & t0,
				double const & T,
				double const & y0,
				double const & h_initial, 
				double const & h_max, 
				double const & final_error,
				int & status,
				std::size_t const & maxSteps){
				
					status=0;
					const std::size_t maxReduction=maxSteps;
					// parameters for decreasing/increasing time step
					double const c1=1.0;
					// I need to have a sufficient decrease of the local error
					// to allow time step coarsening
					double const c2=1./64.;

					double length=T-t0;
					//! Make sure that h allows to reach T
					std::size_t initialNSteps=std::max(static_cast<size_t>(1),static_cast<size_t>(length/h_initial));
					double h=length/initialNSteps;
					// To avoid underflow we need in any case to limit the time step to a positive number
					// Here I allow h to become 128 time smaller than that giving the maximal number of steps
					double h_min = length/(128*maxSteps);
					// SOme counters
					std::size_t stepsCounter(0);
					// Initial data
					double time(t0);
					double y(y0);
					double errorPerTimeStep=final_error/initialNSteps;
					if (initialNSteps>=maxSteps) throw std::runtime_error("RK45: initial time step h too small!");
					std::vector<std::pair<double,double>> solution;
					solution.emplace_back(std::make_pair(t0,y0));
					double localError;
					double newy;
					while (time<T && stepsCounter <maxSteps)
				 {
				//Do a step
				//adjust h if needed for the last step
				if (time + h > T) h = T-time;
				newy = rk45_step(dy,y,time,h,localError);
				while (h> h_min && localError > c1*errorPerTimeStep)
					{
						// half time step
						h /=2;
						errorPerTimeStep /=2;
						newy = rk45_step(dy,y,time,h,localError);
					}
				if (localError>errorPerTimeStep)status=1;
				//! advance
				y = newy;
				time +=h;
				++stepsCounter;
				solution.emplace_back(std::make_pair(time,y));
				//! check if we reached end
				if(localError<c2*errorPerTimeStep && h<h_max)
					{
						// Double step
						h *=2;
						errorPerTimeStep *=2;
					}
				 }
					//handle exceptions
					if(stepsCounter>=maxSteps && time < T)
				 {
				status=2;
				throw std::runtime_error("RK45: Max number of time steps exceeded");
				 }
					return solution;
				}
}
