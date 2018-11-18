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
real_t Grid::Interpolate(const multi_real_t &pos) const{
	if (pos[0]>_geom->Length()[0] || pos[1]>_geom->Length()[1]){
		cout << "Position request out of bounds. Exiting." << endl;
		exit(0);
	}
	multi_real_t length = _geom->Length();
	multi_index_t size = _geom->Size();
	multi_real_t h = _geom->Mesh();

	index_t one = (pos[0]/h[0]-_offset[0])+1;
	index_t two = (pos[1]/h[1]-_offset[1])+1;

	index_t iterValue = one+(size[0]+2)*two;

	Iterator iter = Iterator(_geom,iterValue);

	multi_real_t dist;

	// (possibly negative) x- and y-distance from sampling point
	// to interpolation point pos
	dist[0] = pos[0]-(iter.Right().Pos()[0]-1)*h[0]-_offset[0]*h[0];
	dist[1] = pos[1]-(iter.Top().Pos()[1]-1)*h[1]-_offset[1]*h[1];

	multi_real_t interpx;
	real_t interpy;
	interpx[0] = (Cell(iter.Right())-Cell(iter))/h[0]*dist[0]+Cell(iter.Right());
	interpx[1] = (Cell(iter.Top().Right())-Cell(iter.Top()))/h[0]*dist[0]+Cell(iter.Top().Right());
	interpy = (interpx[1]-interpx[0])/h[1]*dist[1]+interpx[1];
	return interpy;
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
	Iterator *iter = new Iterator(_geom);
	real_t max = Cell(*iter);
	for(iter->First();iter->Valid();iter->Next()){
		if (max < Cell(*iter)) max = Cell(*iter);
	}
	return max;
}
/// Returns the minimal value of the grid
real_t Grid::Min() const{
	Iterator *iter = new Iterator(_geom);
	real_t min = Cell(*iter);
	for(iter->First();iter->Valid();iter->Next()){
		if (min > Cell(*iter)) min = Cell(*iter);
	}
	return min;
}

/// Returns the absolute maximal value
real_t Grid::AbsMax() const{
	Iterator *iter = new Iterator(_geom);
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
