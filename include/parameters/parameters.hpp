#ifndef NAST_PARAMETERS_SIMULATION_PARAMETERS_HPP_
#define NAST_PARAMETERS_SIMULATION_PARAMETERS_HPP_

#include "util/defines.hpp"

#include <iostream>
#include <sstream>

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
		solver(solver_), max_solver_iterations(max_solver_iterations_), eps(eps_), alpha(alpha_), omega(omega_), verbose(verbose_)
	{}
	
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
		
	friend std::ostream& operator<<(std::ostream& os, parameters const& params);
	
	std::string to_string() const
	{
		std::ostringstream stream;
		
		stream << *this;
		
		return stream.str();
	}
};

std::ostream& operator<<(std::ostream& os, parameters const& params)
{
	std::string solver_string = (params.solver == solver_type::Jacobi) ? "Jacobi" : "SOR";
	
	os  << "Parameters:"
		<< "\n\tReynolds " << params.reynolds 
		<< "\nTime integration"
		<< "\n\tinitial_dt " << params.initial_dt
		<< "\n\tt_end " << params.t_end
		<< "\n\tmax_timesteps " << params.max_timesteps
		<< "\n\tsafety factor " << params.tau
		<< "\nSolver"
		<< "\n\ttype " << solver_string
		<< "\n\tmax_solver_iterations " << params.max_solver_iterations 
		<< "\n\ttolerance " << params.eps
		<< "\n\trelaxation factor (SOR) " << params.omega
		<< "\n\tsmoothing factor (stencils) " << params.alpha
		<< "\nMisc"
		<< "\n\tverbose " << params.verbose
		<< "\n";
		
	return os;
}

} //namespace parameters
} //namespace nast

#endif
