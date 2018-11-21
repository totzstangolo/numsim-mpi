#include "iterator.hpp"
#include "geometry.hpp"


// The Iterator Class
Iterator::Iterator(const Geometry *geom){
	_geom = geom;
	_value = 0;
	_valid = true; //no check needed, because the iterator has at least 4 fields
}

Iterator::Iterator(const Geometry *geom, const index_t &value): Iterator(geom){
	_value = value;
	_valid = ((_geom->Size()[0]+2)*(_geom->Size()[1]+2) > value); //Check if the value is even valid
}

///     Returns the current position value
const index_t &Iterator::Value() const{
	return _value;
}
/// Cast operator to convert Iterators to integers
Iterator::operator const index_t &() const{
	return _value;
}
/// Returns the position coordinates
multi_index_t Iterator::Pos() const{
	return { _value % (_geom->Size()[0]+2), _value / (_geom->Size()[0]+2) };
}

/// Sets the iterator to the first element
void Iterator::First(){
	_value = 0;
	_valid = true;
}
/// Goes to the next element of the iterator, disables it if position is end
void Iterator::Next(){
	_value++;
	_valid = ((_geom->Size()[0]+2)*(_geom->Size()[1]+2) > _value); //Check if the value is even valid
}
/// Goes to the second next element of the iterator, disables it if position is end
void Iterator::DoubleNext(){
	_value+= 2;
	_valid = ((_geom->Size()[0]+2)*(_geom->Size()[1]+2) > _value); //Check if the value is even valid
}

/// Checks if the iterator still has a valid value
bool Iterator::Valid() const{
    return _valid;
}

/// Returns an Iterator that is located left from this one.
// if we are at the left boundary, the cell sees itself
Iterator Iterator::Left() const{
	if(_value % (_geom->Size()[0]+2) == 0){
		return *this;
	} else {
		return Iterator(_geom, _value-1);
	}
}

/// Returns an Iterator that is located right from this one
// If we are at the right boundary, the cell sees itself
Iterator Iterator::Right() const{
	if(_value % (_geom->Size()[0]+2) == (_geom->Size()[0]+2) - 1){
		return *this;
	} else {
		return Iterator(_geom,_value+1);
	}
}

/// Returns an Iterator that is located above this one
// If we are at the upper domain boundary, the cell sees itself
Iterator Iterator::Top() const{
	if(_value / (_geom->Size()[0]+2) == (_geom->Size()[1]+2) - 1){
		return *this;
	} else {
		return Iterator(_geom,_value+(_geom->Size()[0]+2));
	}
}

/// Returns an Iterator that is located below this one
// If we are at the lower domain boundary, the cell sees itself
Iterator Iterator::Down() const{
	if(_value / (_geom->Size()[0]+2) == 0){
		return *this;
	} else {
		return Iterator(_geom,_value-(_geom->Size()[0]+2));
	}
}

InteriorIterator::InteriorIterator(const Geometry *geom): Iterator(geom){
	_value = 1 + (_geom->Size()[0]+2);
	_valid = (_value < (((_geom->Size()[1]+2) - 1)*(_geom->Size()[0]+2)));
}

void InteriorIterator::First() {
	// Calling the First method in the constructor might not work properly
	// during construction because it is a part of a virtual
	_value = 1 + (_geom->Size()[0]+2);
	_valid = (_value < (((_geom->Size()[1]+2) - 1)*(_geom->Size()[0]+2)));
}
/// Goes to the next element of the iterator, disables it if position is end
void InteriorIterator::Next() {
	multi_index_t pos = this->Pos();
	if( pos[0] == ((_geom->Size()[0]+2) - 2)){ //are we at the border?
		_value = _value + 3;
	} else { //not at the border
		_value++;
	}
	_valid = (_value < (((_geom->Size()[1]+2) - 1)*(_geom->Size()[0]+2)));
}

void InteriorIterator::DoubleNext() {
	multi_index_t pos = this->Pos();
	if( pos[0] >= (_geom->Size()[0] - 1)){
			_value+= 4;
	}else {
	_value+= 2;
	}
	_valid = (_value < (((_geom->Size()[1]+2) - 1)*(_geom->Size()[0]+2)));
}

BoundaryIterator::BoundaryIterator(const Geometry *geom) : Iterator(geom) {
	_boundary = 0;
}

/// Sets the boundary to iterate
// This also sets the _value to a valid starting point
// Boundary values are as follows:
// bottom: 0 (default)
// left: 1
// top: 2
// right: 3
void BoundaryIterator::SetBoundary(const index_t &boundary){
	_boundary = boundary; //TODO Exception werfen wenn boundary > 3?
	this->First();
}

/// Sets the iterator to the first element
// The first element is always the bottom/left most one
void BoundaryIterator::First(){
	switch(_boundary) {
		case 0:
			_value = 0;
			break;
		case 1:
			_value = 0;
			break;
		case 2:
			_value = ((_geom->Size()[1]+2) - 1)*(_geom->Size()[0]+2);
			break;
		case 3:
			_value = (_geom->Size()[0]+2) - 1;
			break;
	}
	_valid = (_value < (_geom->Size()[0]+2)*(_geom->Size()[1]+2));
}

/// Goes to the next element of the iterator, disables it if position is end
void BoundaryIterator::Next(){
	if (_boundary % 2 == 0){
		_value ++;
	} else {
		_value = _value + (_geom->Size()[0]+2);
	}
	if (_boundary == 0){
		_valid = (_value < (_geom->Size()[0]+2));
	} else {
		_valid = (_value < (_geom->Size()[0]+2)*(_geom->Size()[1]+2));
	}
}
