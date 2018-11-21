/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//------------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include "typedef.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include "compute.hpp"
#include "iterator.hpp"
#include "solver.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>


/** abstract base class for an iterative solver
*/

	/// Constructor of the abstract Solver class
Solver::Solver(const Geometry *geom){
	_geom = geom;
}

	/// Destructor of the Solver Class
Solver::~Solver(){}

	/// Returns the residual at [it] for the pressure-Poisson equation
real_t Solver::localRes(const Iterator &it, const Grid *grid, const Grid *rhs) const{
	  return grid->dxx(it) + grid->dyy(it) - rhs->Cell(it);
}

	/// Constructs an actual SOR solver
SOR::SOR(const Geometry *geom, const real_t &omega): Solver(geom) {
	_omega = omega;
	if(_omega == 0) {
	// optimized omega after lecture (middle over h_x and h_y)
	_omega = 2 / (1 + sin(M_PI*((_geom->Mesh()[0] + _geom->Mesh()[1]) /2)));
	}
}

  /// Destructor
SOR::~SOR(){

}

  /// Returns the total residual and executes a solver cycle
  // @param grid current pressure values
  // @param rhs right hand side
real_t SOR::Cycle(Grid *grid, const Grid *rhs) const{
	//pre compute const. factor
	real_t _const_cor = (_geom->Mesh()[0] * _geom->Mesh()[0] * _geom->Mesh()[1] * _geom->Mesh()[1] /
			(2.0 * (_geom->Mesh()[0] * _geom->Mesh()[0] + _geom->Mesh()[1] * _geom->Mesh()[1])));
	real_t _global_normed_res = 0.0;
	real_t _l_res = 0.0;
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()) {
		_l_res = localRes(intIterator, grid, rhs);
		// norm residuum
		_global_normed_res += _l_res * _l_res;
		grid->Cell(intIterator) += _omega * _const_cor * _l_res;
	}
	return _global_normed_res;
}

	//Constructs an RedOrBlackSOR-Solver
RedOrBlackSOR::RedOrBlackSOR(const Geometry *geom, const real_t &omega) : SOR(geom, omega){

}

	//Destructor
RedOrBlackSOR::~RedOrBlackSOR(){

}

real_t RedOrBlackSOR::RedCycle(Grid *grid, const Grid *rhs) const{
	//pre compute const. factor
	real_t _const_cor = (_geom->Mesh()[0] * _geom->Mesh()[0] * _geom->Mesh()[1] * _geom->Mesh()[1] /
			(2.0 * (_geom->Mesh()[0] * _geom->Mesh()[0] + _geom->Mesh()[1] * _geom->Mesh()[1])));
	real_t _global_normed_res = 0.0;
	real_t _l_res = 0.0;
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.DoubleNext()) {
		_l_res = localRes(intIterator, grid, rhs);
		// norm residuum
		_global_normed_res += _l_res * _l_res;
		grid->Cell(intIterator) += _omega * _const_cor * _l_res;
	}
	return _global_normed_res;
}

real_t RedOrBlackSOR::BlackCycle(Grid *grid, const Grid *rhs) const{
	//pre compute const. factor
	real_t _const_cor = (_geom->Mesh()[0] * _geom->Mesh()[0] * _geom->Mesh()[1] * _geom->Mesh()[1] /
			(2.0 * (_geom->Mesh()[0] * _geom->Mesh()[0] + _geom->Mesh()[1] * _geom->Mesh()[1])));
	real_t _global_normed_res = 0.0;
	real_t _l_res = 0.0;
	InteriorIterator intIterator(_geom);
	intIterator.First();
	for(intIterator.Next(); intIterator.Valid(); intIterator.DoubleNext()) {
		_l_res = localRes(intIterator, grid, rhs);
		// norm residuum
		_global_normed_res += _l_res * _l_res;
		grid->Cell(intIterator) += _omega * _const_cor * _l_res;
	}
	return _global_normed_res;
}
