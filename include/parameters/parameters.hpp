#ifndef NAST_PARAMETERS_SIMULATION_PARAMETERS_HPP_
#define NAST_PARAMETERS_SIMULATION_PARAMETERS_HPP_

#include "util/defines.hpp"

namespace nast { namespace parameters {

enum solver_type {
	Jacobi,
	SOR
};

struct parameters
{
	parameters(Real reynolds_ = 1000, Real initial_dt_ = 0.002, Real t_end_ = 1, std::size_t max_timesteps_ = 100, Real tau_ = 0.5,
		solver_type solver_ = solver_type::SOR, std::size_t max_solver_iterations_ = 100, Real eps_ = 1e-3, Real alpha_ = 0.9, Real omega_ = 1.7, bool verbose_ = false) 
	: reynolds(reynolds_), initial_dt(initial_dt_), t_end(t_end_), max_timesteps(max_timesteps_), tau(tau_),
		solver(solver_), max_solver_iterations(max_solver_iterations_), eps(eps_), alpha(alpha_), omega(omega_), verbose(verbose_) {}
	
	Real reynolds;
	
	Real initial_dt;
	Real t_end;
	std::size_t max_timesteps;
	
	// for stability
	Real tau;
	
	solver_type solver;
	
	std::size_t max_solver_iterations;
	Real eps;
	
	// for derivatives
	Real alpha;
	
	// relaxation for sor
	Real omega;
	
	bool verbose;
	
	friend std::ostream& operator<<(std::ostream& os, parameters const& params)
    {
		std::string solver_string = (params.solver == solver_type::Jacobi) ? "Jacobi" : "SOR";
		
        os  << "Parameters:"
			<< "\n\t Reynolds " << params.reynolds 
			<< "\n Time integration"
			<< "\n\t initial_dt " << params.initial_dt
			<< "\n\t t_end " << params.t_end
			<< "\n\t max_timesteps " << params.max_timesteps
			<< "\n\t safety factor " << params.tau
			<< "\n Solver"
			<< "\n\t type " << solver_string
			<< "\n\t max_solver_iterations " << params.max_solver_iterations 
			<< "\n\t tolerance " << params.eps
			<< "\n\t relaxation factor (SOR) " << params.omega
			<< "\n\t smoothing factor (stencils) " << params.alpha
			<< "\n Misc"
			<< "\n\t verbose " << params.verbose
			<< "\n";
			
        return os;
    }

};

} //namespace parameters
} //namespace nast

#endif
