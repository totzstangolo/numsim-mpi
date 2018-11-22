#include "grid.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include <iostream>
#include <cmath>
using namespace std;
/// Constructs a grid based on a geometry
Grid::Grid(const Geometry *geom){
	_geom = geom;
	_data = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_offset[0] = 0;
	_offset[1] = 0;
}

/// Constructs a grid based on a geometry with an offset
// @param geom   Geometry information
// @param offset distance of staggered grid point to cell's anchor point
//               (anchor point = lower left corner)
Grid::Grid(const Geometry *geom, const multi_real_t &offset): Grid(geom) {
	_geom = geom;
	_data = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_offset = offset;
}

/// Deletes the grid
Grid::~Grid(){
	if(_data!=nullptr){
		delete[] _data;
	}
}

///     Initializes the grid with a value
void Grid::Initialize(const real_t &value){
	index_t gridsize = (_geom->Size()[0]+2)*(_geom->Size()[1]+2);
	for(index_t i = 0; i < gridsize; i++){
		_data[i] = value;
	}
}

/// Write access to the grid cell at position [it]
real_t &Grid::Cell(const Iterator &it){
	return _data[it];
}
/// Read access to the grid cell at position [it]
const real_t &Grid::Cell(const Iterator &it) const{
	return _data[it];
}

/// Interpolate the value at a arbitrary position
// pos a value in the domain e.g. between 0 and 1
//Offset: p = {0.0,0.0}, u = {0.0,0.5}, v = {0.5,0.0}
real_t Grid::Interpolate (const multi_real_t& pos) const {
	index_t idx[DIM][2];
	real_t alpha[DIM];
	multi_index_t _size;
	_size[0] = 1;
	_size[1] = _geom->Size()[0];
	for (uint32_t i=0;i<DIM;++i) {
		if (pos[i] > _offset[i]) {
			alpha[i]  = (pos[i] - _offset[i])/_geom->Mesh()[i];
			idx[i][0] = floor(alpha[i]);
			alpha[i] -= idx[i][0];
		}
		else {
			alpha[i] = 0;
			idx[i][0] = 0;
		}
		idx[i][1] = idx[i][0] + 1;
		if (idx[i][1] >= _geom->Size()[i])
			idx[i][1] = idx[i][0];
    }
	real_t res = 0.0;
    real_t weight;
    index_t node;
    for (uint32_t i=0;i<(1<<DIM);++i) {
		weight = 1.0;
		node = 0;
		for (uint32_t d=0;d<DIM;++d) {
			if (i&(1<<d)) {
				weight *= alpha[d];
				node += idx[d][1]*_size[d];
			} else {
				weight *= (1.0 - alpha[d]);
				node += idx[d][0]*_size[d];
			}
		}
		res += _data[node]*weight;
    }
    return res;
}


/// Computes the left-sided difference quotient in x-dim at [it]
real_t Grid::dx_l(const Iterator &it) const{
	return (Cell(it)-Cell(it.Left()))/_geom->Mesh()[0];
}
/// Computes the right-sided difference quotient in x-dim at [it]
real_t Grid::dx_r(const Iterator &it) const{
	return (Cell(it.Right())-Cell(it))/_geom->Mesh()[0];
}
/// Computes the left-sided difference quotient in y-dim at [it]
real_t Grid::dy_l(const Iterator &it) const{
	return (Cell(it)-Cell(it.Down()))/_geom->Mesh()[1];
}
/// Computes the right-sided difference quotient in x-dim at [it]
real_t Grid::dy_r(const Iterator &it) const{
	return (Cell(it.Top())-Cell(it))/_geom->Mesh()[1];
}
/// Computes the central difference quotient of 2nd order in x-dim at [it]
real_t Grid::dxx(const Iterator &it) const{
	return (Cell(it.Right())-2*Cell(it)+Cell(it.Left()))/(_geom->Mesh()[0]*_geom->Mesh()[0]);
}
/// Computes the central difference quotient of 2nd order in y-dim at [it]
real_t Grid::dyy(const Iterator &it) const{
	return (Cell(it.Top())-2*Cell(it)+Cell(it.Down()))/(_geom->Mesh()[1]*_geom->Mesh()[1]);
}

/// Computes u*du/dx with the donor cell method
real_t Grid::DC_udu_x(const Iterator &it, const real_t &alpha) const{
	real_t cell_it = Cell(it);
	real_t term_1 = cell_it + Cell(it.Right());
	real_t term_2 = Cell(it.Left()) + cell_it;
	real_t term_3 = cell_it - Cell(it.Right());
	real_t term_4 = Cell(it.Left()) - cell_it;
	return (0.25*(((term_1*term_1)-(term_2*term_2))
		+ alpha*(abs(term_1)*term_3-abs(term_2)*term_4))/_geom->Mesh()[0]);
}
/// Computes v*du/dy with the donor cell method
real_t Grid::DC_vdu_y(const Iterator &it, const real_t &alpha, const Grid *v) const{
	real_t cell_it_u = Cell(it);
	real_t cell_it_v = v->Cell(it);
	real_t term_1 = cell_it_v + v->Cell(it.Right());
	real_t term_2 = cell_it_u + Cell(it.Top());
	real_t term_3 = v->Cell(it.Down()) + v->Cell(it.Right().Down());
	real_t term_4 = Cell(it.Down()) + cell_it_u;
	real_t term_5 = cell_it_u - Cell(it.Top());
	real_t term_6 = Cell(it.Down()) - cell_it_u;
	return (0.25*(((term_1*term_2)-(term_3*term_4))
		+ alpha*(abs(term_1)*term_5-abs(term_3)*term_6))/_geom->Mesh()[1]);
}
/// Computes u*dv/dx with the donor cell method
real_t Grid::DC_udv_x(const Iterator &it, const real_t &alpha, const Grid *u) const{
	real_t cell_it_v = Cell(it);
	real_t cell_it_u = u->Cell(it);
	real_t term_1 = cell_it_u + u->Cell(it.Top());
	real_t term_2 = cell_it_v + Cell(it.Right());
	real_t term_3 = u->Cell(it.Left()) + u->Cell(it.Left().Top());
	real_t term_4 = Cell(it.Left()) + cell_it_v;
	real_t term_5 = cell_it_v - Cell(it.Right());
	real_t term_6 = Cell(it.Left()) - cell_it_v;
	return (0.25*(((term_1*term_2)-(term_3*term_4))
		+ alpha*(abs(term_1)*term_5-abs(term_3)*term_6))/_geom->Mesh()[0]);
}
/// Computes v*dv/dy with the donor cell method
real_t Grid::DC_vdv_y(const Iterator &it, const real_t &alpha) const{
	real_t cell_it = Cell(it);
	real_t term_1 = cell_it + Cell(it.Top());
	real_t term_2 = Cell(it.Down()) + cell_it;
	real_t term_3 = cell_it - Cell(it.Top());
	real_t term_4 = Cell(it.Down()) - cell_it;
	return (0.25*(((term_1*term_1)-(term_2*term_2))
		+ alpha*(abs(term_1)*term_3-abs(term_2)*term_4))/_geom->Mesh()[1]);
}

/// Returns the maximal value of the grid
real_t Grid::Max() const{
	InteriorIterator *iter = new InteriorIterator(_geom);
	real_t max = Cell(*iter);
	for(iter->First();iter->Valid();iter->Next()){
		if (max < Cell(*iter)) max = Cell(*iter);
	}
	return max;
}
/// Returns the minimal value of the grid
real_t Grid::Min() const{
	InteriorIterator *iter = new InteriorIterator(_geom);
	real_t min = Cell(*iter);
	for(iter->First();iter->Valid();iter->Next()){
		if (min > Cell(*iter)) min = Cell(*iter);
	}
	return min;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const{
	InteriorIterator *iter = new InteriorIterator(_geom);
	real_t max = Cell(*iter);
	for(iter->First();iter->Valid();iter->Next()){
		if (max < abs(Cell(*iter))) max = abs(Cell(*iter));
	}
	return max;
}

/// Returns a pointer to the raw data
real_t *Grid::Data(){
	return _data;
}

const multi_real_t &Grid::getOffset() const{
	return _offset;
}

/// Return a pointer to the Geometry
const Geometry *Grid::getGeometry() const{
	return _geom;
}
