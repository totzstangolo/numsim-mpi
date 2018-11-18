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

#include "geometry.hpp"
#include "iterator.hpp"
#include "grid.hpp"
#include "communicator.hpp"
#include <stdio.h>
#include <string.h>
//------------------------------------------------------------------------------
Geometry::Geometry() : _comm(NULL) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;

  _bsize = _size;
  _blength = _length;
}
//------------------------------------------------------------------------------
Geometry::Geometry(const Communicator *comm) : _comm(comm) {
  _length[0] = 1.0;
  _length[1] = 1.0;
  _size[0] = 128;
  _size[1] = 128;
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];
  _pressure = 0.0;
  _velocity[0] = 1.0;
  _velocity[1] = 0.0;

  _bsize[0] = _size[0] / _comm->ThreadDim()[0] + 2;
  _bsize[1] = _size[1] / _comm->ThreadDim()[1] + 2;
  if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
    _bsize[0] += _size[0] % _comm->ThreadDim()[0];
  if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
    _bsize[1] += _size[1] % _comm->ThreadDim()[0];

  _blength[0] = _h[0] * (_bsize[0] - 2);
  _blength[1] = _h[1] * (_bsize[1] - 2);

  // create boundary halo
  _size[0] += 2;
  _size[1] += 2;
}
//------------------------------------------------------------------------------
void Geometry::Load(const char *file) {
  FILE *handle = fopen(file, "r");
  double inval[2];
  char name[20];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s =", name))
      continue;
    if (strcmp(name, "size") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _size[0] = inval[0];
        _size[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "length") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _length[0] = inval[0];
        _length[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "velocity") == 0) {
      if (fscanf(handle, " %lf %lf\n", &inval[0], &inval[1])) {
        _velocity[0] = inval[0];
        _velocity[1] = inval[1];
      }
      continue;
    }
    if (strcmp(name, "pressure") == 0) {
      if (fscanf(handle, " %lf\n", &inval[0]))
        _pressure = inval[0];
      continue;
    }
  }
  fclose(handle);
  _h[0] = _length[0] / _size[0];
  _h[1] = _length[1] / _size[1];

  _blength = _length;

  if (_comm) {
    _bsize[0] = _size[0] / _comm->ThreadDim()[0] + 2;
    _bsize[1] = _size[1] / _comm->ThreadDim()[1] + 2;
    if (_comm->ThreadIdx()[0] == _comm->ThreadDim()[0] - 1)
      _bsize[0] += _size[0] % _comm->ThreadDim()[0];
    if (_comm->ThreadIdx()[1] == _comm->ThreadDim()[1] - 1)
      _bsize[1] += _size[1] % _comm->ThreadDim()[0];

    _blength[0] = _h[0] * (_bsize[0] - 2);
    _blength[1] = _h[1] * (_bsize[1] - 2);
  }

  _size[0] += 2;
  _size[1] += 2;

  if (!_comm)
    _bsize = _size;
}
//------------------------------------------------------------------------------
const multi_index_t &Geometry::Size() const { return _bsize; }
//------------------------------------------------------------------------------
const multi_index_t &Geometry::TotalSize() const { return _size; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Length() const { return _blength; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::TotalLength() const { return _length; }
//------------------------------------------------------------------------------
const multi_real_t &Geometry::Mesh() const { return _h; }
//------------------------------------------------------------------------------
void Geometry::Update_U(Grid *u) const {
  // to be implemented
}
//------------------------------------------------------------------------------
void Geometry::Update_V(Grid *v) const {
  // to be implemented
}
//------------------------------------------------------------------------------
void Geometry::Update_P(Grid *p) const {
  // to be implemented
}
//------------------------------------------------------------------------------
