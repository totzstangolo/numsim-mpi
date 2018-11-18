
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
#include <ctime>
using namespace std;

Compute::Compute(const Geometry *geom, const Parameter *param,
				 const Communicator *comm){
	_geom = geom;
	_param = param;

	// Initializing time, dtlimit and epslimit_
	_t = 0.0;
	multi_real_t h = _geom->Mesh();
#ifdef USE_OPT_DT
	_dtlimit = h[0] * h[0] * h[1] * h[1] * _param->Re() /
			( 2 * ((h[0] * h[0]) + (h[1] * h[1])));
#else
	_dtlimit = _param->Dt();
#endif  // USE_OPT_DT
	_epslimit = _param->Eps() * _param->Eps() * _geom->Size()[0] * _geom->Size()[1];

	// creating grids with offset
	multi_real_t compute_offset_x;
	compute_offset_x[0] = 0.0;
	compute_offset_x[1] = -0.5 * h[1];
	_u = new Grid(_geom, compute_offset_x);
	_F = new Grid(_geom, compute_offset_x);
	_u->Initialize(0);
	_geom->Update_U(_u);

	multi_real_t compute_offset_y;
	compute_offset_y[0] = -0.5 * h[0];
	compute_offset_y[1] = 0.0;
	_v = new Grid(_geom, compute_offset_y);
	_G = new Grid(_geom, compute_offset_y);
	_v->Initialize(0);
	_geom->Update_V(_v);

	multi_real_t compute_offset_p;
	compute_offset_p[0] = -0.5 * h[0];
	compute_offset_p[1] = -0.5 * h[1];
	_p = new Grid(_geom, compute_offset_p);
	_rhs = new Grid(_geom, compute_offset_p);
	_tmp = new Grid(_geom, compute_offset_p);
	_p->Initialize(0);

	//create solver (used script omega, not param omega)
	_solver = new SOR(_geom, _param->Omega());


}

Compute::~Compute(){
	delete _u;
	delete _F;
	delete _v;
	delete _G;
	delete _p;
	delete _rhs;
	delete _tmp;
	delete _solver;
}

void Compute::TimeStep(bool printInfo){


	// Compute like in script page 23
	//compute dt
#ifdef USE_OPT_DT
	multi_real_t max;
	max[0] = _u->Max();
	max[1] = _v->Max();
	real_t dt;
	real_t limit = _dtlimit;
	if(max[0] > 0 && max[1] > 0){
	   if(_geom->Mesh()[0]/max[0] < _geom->Mesh()[0]/max[0]){
	 	limit = _geom->Mesh()[0]/max[0];
	   }else{
		limit = _geom->Mesh()[1]/max[1];
	   }
	}
	if(_dtlimit > limit){
		dt = limit * _param->Tau();
	}else{
		dt = _dtlimit * _param->Tau();
	}
#else
	real_t dt = _dtlimit * _param->Tau();
#endif  // USE_OPT_DT

	if(printInfo) {
		printf("dt: %f \n", dt);
	}


	// compute FG and update bound.
	MomentumEqu(dt);
	_geom->Update_U(_F);
	_geom->Update_V(_G);

	// compute rhs and update bound.
	RHS(dt);
	_geom->Update_P(_rhs);

	// Solver, relative eps had bad performance
	index_t index = 0;
	real_t res = 1;
	while(index < _param->IterMax() && res > _epslimit) {
		index++;
		res = _solver->Cycle(_p, _rhs);
		_geom->Update_P(_p);
	}
	if(printInfo) {
		printf("iterations: %d / %d \n", index, _param->IterMax());
	printf("pmax: %f\n", _p->Max());
	}

	//compute uv and update bound.
	NewVelocities(dt);
	_geom->Update_U(_u);
	_geom->Update_V(_v);

	_t += dt;
	if(printInfo)
		printf("time: %f \n", _t);

}

  /// Returns the simulated time in total
const real_t& Compute::GetTime() const{
    return _t;
}

  /// Returns the pointer to U
const Grid* Compute::GetU() const{
    return _u;
}
  /// Returns the pointer to V
const Grid* Compute::GetV() const{
    return _v;
}
  /// Returns the pointer to P
const Grid* Compute::GetP() const{
    return _p;
}
  /// Returns the pointer to RHS
const Grid* Compute::GetRHS() const{
    return _rhs;
}

/// Computes and returns the absolute velocity
const Grid* Compute::GetVelocity(){
	multi_real_t mid;
	Iterator iterator(_geom);
	for (iterator.First(); iterator.Valid(); iterator.Next()){
		mid[0] = _u->Cell(iterator) + _u->Cell(iterator.Left());
		mid[1] = _v->Cell(iterator) + _v->Cell(iterator.Down());
		_tmp->Cell(iterator) = sqrt((mid[0] * mid[0]) + (mid[1] * mid[1]));
	}
	return _tmp;
}

/// Computes and returns the vorticity
const Grid* Compute::GetVorticity(){
	// TODO real implementation needed
    Grid *grid = new Grid(_geom);
	return grid;
}

/// Computes and returns the stream line values
const Grid* Compute::GetStream(){
	// TODO real implementation needed
    Grid *grid = new Grid(_geom);
    return grid;
}

void Compute::NewVelocities(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_u->Cell(intIterator) = _F->Cell(intIterator) - (dt * _p->dx_r(intIterator));
		_v->Cell(intIterator) = _G->Cell(intIterator) - (dt * _p->dy_r(intIterator));
	}
}

/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_F->Cell(intIterator) = _u->Cell(intIterator) +
			(
				 (_u->dxx(intIterator) + _u->dyy(intIterator))/ _param->Re()
				- _u->DC_udu_x(intIterator, _param->Alpha())
				- _u->DC_vdu_y(intIterator, _param->Alpha(), _v
			)) * dt;
		_G->Cell(intIterator) = _v->Cell(intIterator) +
			(
			 (_v->dyy(intIterator) + _v->dxx(intIterator))/ _param->Re()
				- _v->DC_vdv_y(intIterator, _param->Alpha())
				- _v->DC_udv_x(intIterator, _param->Alpha(), _u)
				) * dt;
	}

}
/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_rhs->Cell(intIterator) = (_F->dx_l(intIterator) + _G->dy_l(intIterator)) / dt;
	}
}
