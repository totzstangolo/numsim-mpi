#include "iterator.hpp"
#include "geometry.hpp"
//------------------------------------------------------------------------------
Iterator::Iterator (const Geometry* geom) : _geom(geom) {
	_valid = true;
	_value = 0;
}
//------------------------------------------------------------------------------
Iterator::Iterator (const Geometry* geom, const index_t& value) : _geom(geom) {
	_value = value;
	_valid = true;
	if (_value >= _geom->Size()[0]*_geom->Size()[1])
		_valid = false;
}
//------------------------------------------------------------------------------
const index_t& Iterator::Value () const {
	return _value;
}
//------------------------------------------------------------------------------
Iterator::operator const index_t&() const {
	return _value;
}
//------------------------------------------------------------------------------
multi_index_t Iterator::Pos () const {
	multi_index_t pos;
	pos[0] = _value%_geom->Size()[0];
	pos[1] = _value/_geom->Size()[0];
	return pos;
}
//------------------------------------------------------------------------------
void Iterator::First () {
	_value = 0;
	_valid = true;
}
//------------------------------------------------------------------------------
void Iterator::Next () {
	if (!_valid) return;
	++_value;
	if (_value >= _geom->Size()[0]*_geom->Size()[1])
		_valid = false;
}
//------------------------------------------------------------------------------
bool Iterator::Valid () const {
	return _valid;
}
//------------------------------------------------------------------------------
Iterator Iterator::Left () const {
    if (_value%_geom->Size()[0] == 0) return *this;
    return Iterator(_geom,_value - 1);
}
//------------------------------------------------------------------------------
Iterator Iterator::Right () const {
	if (_value%_geom->Size()[0] == (_geom->Size()[0] - 1)) return *this;
    return Iterator(_geom,_value + 1);
}
//------------------------------------------------------------------------------
Iterator Iterator::Top () const {
	if (_value/_geom->Size()[0] == (_geom->Size()[1] - 1)) return *this;
    return Iterator(_geom,_value + _geom->Size()[0]);
}
//------------------------------------------------------------------------------
Iterator Iterator::Down () const {
	if (_value/_geom->Size()[0] == 0) return *this;
    return Iterator(_geom,_value - _geom->Size()[0]);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
InteriorIterator::InteriorIterator (const Geometry* geom) : Iterator(geom) {
}
//------------------------------------------------------------------------------
void InteriorIterator::First () {
	_value = _geom->Size()[0] + 1;
	_valid = true;
}
//------------------------------------------------------------------------------
void InteriorIterator::Next  () {
    if (!_valid) return;
    ++_value;
    if (_value%_geom->Size()[0] >= (_geom->Size()[1] - 1)) _value += 2;
    if (_value >= _geom->Size()[0]*(_geom->Size()[1]-1)) _valid = false;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
BoundaryIterator::BoundaryIterator (const Geometry* geom) : Iterator(geom) {
	_boundary = 0;
}
//------------------------------------------------------------------------------
void BoundaryIterator::SetBoundary (const index_t& boundary) {
	_boundary = boundary%4;
	_valid = false;
}
//------------------------------------------------------------------------------
void BoundaryIterator::First () {
	_valid = true;
	switch (_boundary) {
	case 0:
		_value = 0;
		break;
	case 1:
		_value = 0;
		break;
	case 2:
		_value = _geom->Size()[0]*(_geom->Size()[1] - 1);
		break;
	case 3:
		_value = _geom->Size()[0] - 1;
		break;
	default:
		_boundary = 0;
		_value = 0;
		break;
	};
}
//------------------------------------------------------------------------------
void BoundaryIterator::Next () {
	if (!_valid) return;
	switch (_boundary) {
	case 0:
		++_value;
		if (_value >= _geom->Size()[0]) _valid = false;
		break;
	case 1:
		_value += _geom->Size()[0];
		if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false;
		break;
	case 2:
		++_value;
		if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false;
		break;
	case 3:
		_value += _geom->Size()[0];
		if (_value >= _geom->Size()[0]*_geom->Size()[1]) _valid = false;
		break;
	default:
		++_value;
		if (_value >= _geom->Size()[0]) _valid = false;
		break;
	};
}

void InteriorIterator::DoubleNext() {
	multi_index_t pos = this->Pos();
	if( pos[0] >= (_geom->Size()[0] - 2)){
			_value+= 3;
	}else {
	_value+= 2;
	}
	_valid = (_value < (((_geom->Size()[1]) - 1)*(_geom->Size()[0])));
}


void Iterator::DoubleNext(){
	_value+= 2;
	_valid = ((_geom->Size()[0])*(_geom->Size()[1]) > _value); //Check if the value is even valid
}
//------------------------------------------------------------------------------
