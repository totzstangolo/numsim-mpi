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

#include "communicator.hpp"
#include "vtk.hpp"
#include <cstring>
#include <cstdio>
//------------------------------------------------------------------------------
uint32_t VTK::_cnt = 0;
//------------------------------------------------------------------------------
VTK::VTK(const multi_real_t &h, const multi_real_t &fieldlength,
      const multi_real_t &globalFieldLength, const int &rank, const int &size,
         const multi_index_t &fieldDims)
    : _h(h), _fieldDims(fieldDims), _rank(rank), _size(size) {
  _fieldsize[0] = fieldlength[0]/h[0] + 2;
  _fieldsize[1] = fieldlength[1]/h[1] + 2;
#if DIM == 3
  _fieldsize[2] = fieldlength[2]/h[2] + 2;
#endif
  _globalFieldSize[0] = globalFieldLength[0]/h[0] + 2;
  _globalFieldSize[1] = globalFieldLength[1]/h[1] + 2;
#if DIM == 3
  _globalFieldSize[2] = globalFieldLength[2]/h[2] + 2;
#endif
  _offset = multi_real_t(0.0);
  _handle = NULL;
  _masterHandle = NULL;
}
//------------------------------------------------------------------------------
VTK::VTK(const multi_real_t &h, const multi_real_t &fieldlength,
      const multi_real_t &globalFieldLength, const multi_real_t &offset,
         const int &rank, const int &size, const multi_index_t &fieldDims)
    : _h(h), _fieldDims(fieldDims), _offset(offset), _rank(rank), _size(size) {
  _fieldsize[0] = fieldlength[0]/h[0] + 2;
  _fieldsize[1] = fieldlength[1]/h[1] + 2;
#if DIM == 3
  _fieldsize[2] = fieldlength[2]/h[2] + 2;
#endif
  _globalFieldSize[0] = globalFieldLength[0]/h[0] + 2;
  _globalFieldSize[1] = globalFieldLength[1]/h[1] + 2;
#if DIM == 3
  _globalFieldSize[2] = globalFieldLength[2]/h[2] + 2;
#endif
  _handle = NULL;
  _masterHandle = NULL;
}

//------------------------------------------------------------------------------

void VTK::Init(const char *path) {
  // Init may only be done once
  if (_handle)
    return;

  // create filename string and build file handles from this.
  int flength = strlen(path) + 20;
  char *filename;
  filename = new char[flength];
  if (strlen(path)) {
    sprintf(filename, "%s_%05i_%04i.vtr", path, _cnt, _rank);
  } else {
    sprintf(filename, "%s_%05i_%04i.vtr", "field", _cnt, _rank);
  }
  _handle = fopen(filename, "w");

  if (_rank == 0) {
    if (strlen(path)) {
      sprintf(filename, "%s_%05i.pvtr", path, _cnt);
    } else {
      sprintf(filename, "%s_%05i.pvtr", "field", _cnt);
    }
    _masterHandle = fopen(filename, "w");
  }

  delete[] filename;

  // write header information in slave file
  fprintf(_handle, "<?xml version=\"1.0\"?>\n");
  fprintf(_handle, "<VTKFile type=\"RectilinearGrid\">\n");

  // calculate shifting to get piece extents right
  multi_index_t shift;
  shift[0] = _rank % _fieldDims[0];
  shift[1] = _rank / _fieldDims[0];
  fprintf(_handle, "<RectilinearGrid WholeExtent=\"%i %i %i %i %i %i\""
                   " GhostLevel=\"0\">\n",
          shift[0] * (_fieldsize[0] - 2), (shift[0] + 1) * (_fieldsize[0] - 2),
          shift[1] * (_fieldsize[1] - 2), (shift[1] + 1) * (_fieldsize[1] - 2),
          (DIM == 3 ? _fieldsize[2] - 1 : 0),
          (DIM == 3 ? _fieldsize[2] - 1 : 0));
  fprintf(_handle, "<Piece Extent=\"%i %i %i %i %i %i\">\n",
          shift[0] * (_fieldsize[0] - 2), (shift[0] + 1) * (_fieldsize[0] - 2),
          shift[1] * (_fieldsize[1] - 2), (shift[1] + 1) * (_fieldsize[1] - 2),
          (DIM == 3 ? _fieldsize[2] - 1 : 0),
          (DIM == 3 ? _fieldsize[2] - 1 : 0));

  // print coordinates
  fprintf(_handle, "<Coordinates>\n");
  fprintf(_handle, "<DataArray type=\"Float64\" format=\"ascii\">\n");
  for (uint32_t x = 0; x < _fieldsize[0] - 1; ++x) {
    fprintf(_handle, "%le ", (double)x * _h[0] + _offset[0]);
  }
  fprintf(_handle, "</DataArray>\n");

  fprintf(_handle, "<DataArray type=\"Float64\" format=\"ascii\">\n");
  for (uint32_t y = 0; y < _fieldsize[1] - 1; ++y) {
    fprintf(_handle, "%le ", (double)y * _h[1] + _offset[1]);
  }
  fprintf(_handle, "</DataArray>\n");

  if (DIM == 3) {
    fprintf(_handle, "<DataArray type=\"Float64\" format=\"ascii\">\n");
    for (uint32_t z = 0; z < (DIM == 3 ? _fieldsize[2] : 1); ++z) {
      fprintf(_handle, "%le ", (DIM == 3 ? (double)z * _h[2] + _offset[2] : 0));
    }
    fprintf(_handle, "</DataArray>\n");
  } else {
    fprintf(_handle, "<DataArray type=\"Float64\" format=\"ascii\">\n");
    fprintf(_handle, "0 0 </DataArray>\n");
  }
  fprintf(_handle, "</Coordinates>\n");

  fprintf(_handle, "<CellData>\n");

  // write header information in master file
  if (_rank == 0) {
    fprintf(_masterHandle, "<?xml version=\"1.0\"?>\n");
    fprintf(_masterHandle, "<VTKFile type=\"PRectilinearGrid\">\n");

    // define whole domain extent
    fprintf(
        _masterHandle,
        "<PRectilinearGrid WholeExtent=\"0 %i 0 %i 0 %i\" GhostLevel=\"0\">\n",
        _globalFieldSize[0] - 2, _globalFieldSize[1] - 2,
        (DIM == 3 ? _globalFieldSize[2] - 2 : 0));

    // announce arrays in 3 dimensions
    fprintf(_masterHandle, "<PCoordinates>\n");
    for (uint32_t i = 0; i < 3; ++i) {
      fprintf(_masterHandle, "<PDataArray type=\"Float64\"/>\n");
    }
    fprintf(_masterHandle, "</PCoordinates>\n");

    // announce pieces and their respective extents
    int procCtr = 0;
    for (uint32_t y = 0; y < _fieldDims[1]; ++y) {
      for (uint32_t x = 0; x < _fieldDims[0]; ++x) {
        // create filename string and build file handle from this.
        filename = new char[flength];
        sprintf(filename, "%s_%05i_%04i.vtr", path, _cnt, procCtr);
        char *fName = &filename[4];
        fprintf(_masterHandle,
                "<Piece Extent=\"%i %i %i %i %i %i\" Source=\"%s\"/>\n",
                x * (_fieldsize[0] - 2), (x + 1) * (_fieldsize[0] - 2),
                y * (_fieldsize[1] - 2), (y + 1) * (_fieldsize[1] - 2), 0, 0,
                fName);
        ++procCtr;
        delete[] filename;
      }
    }

    // begin announcing payload
    fprintf(_masterHandle, "<PCellData>\n");
  }
}

//------------------------------------------------------------------------------

void VTK::Finish() {
  if (_rank == 0) {
    if (!_handle || !_masterHandle) {
      return;
    }
  } else {
    if (!_handle) {
      return;
    }
  }

  // close open xml tags
  fprintf(_handle, "</PointData>\n");
  fprintf(_handle, "</Piece>\n");
  fprintf(_handle, "</RectilinearGrid>\n");
  fprintf(_handle, "</VTKFile>\n");

  fclose(_handle);
  _handle = NULL;

  if (_rank == 0) {
    fprintf(_masterHandle, "</PPointData>\n");
    fprintf(_masterHandle, "</PRectilinearGrid>\n");
    fprintf(_masterHandle, "</VTKFile>\n");

    fclose(_masterHandle);
    _masterHandle = NULL;
  }

  _cnt++;
}

//------------------------------------------------------------------------------

void VTK::SwitchToPointData() {
  if (_rank == 0) {
    if (!_handle || !_masterHandle) {
      return;
    }
  } else {
    if (!_handle) {
      return;
    }
  }

  fprintf(_handle, "</CellData>\n");
  fprintf(_handle, "<PointData>\n");

  if (_rank == 0) {
    fprintf(_masterHandle, "</PCellData>\n");
    fprintf(_masterHandle, "<PPointData>\n");
  }
}

//------------------------------------------------------------------------------

void VTK::AddCellScalar(const char *title, const Grid *grid) {
  if (_rank == 0) {
    if (!_handle || !_masterHandle) {
      return;
    }
  } else {
    if (!_handle) {
      return;
    }
  }

  // slave file
  // header
  fprintf(_handle,
          "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\">\n", title);

  multi_real_t pos;
  for (uint32_t y = 1; y < (_fieldsize[1] - 1); ++y) {
    pos[1] = (double)y * _h[1];
    for (uint32_t x = 1; x < (_fieldsize[0] - 1); ++x) {
      pos[0] = (double)x * _h[0];
      fprintf(_handle, "%le ", grid->Interpolate(pos));
    }
    fprintf(_handle, "\n");
  }

  fprintf(_handle, "</DataArray>\n");

  // print info about data-segment into parallel master file
  if (_rank == 0) {
    fprintf(_masterHandle,
            "<PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\"/>\n",
            title);
  }
}

//------------------------------------------------------------------------------

void VTK::AddCellField(const char *title, const Grid *v1, const Grid *v2) {
  if (_rank == 0) {
    if (!_handle || !_masterHandle) {
      return;
    }
  } else {
    if (!_handle) {
      return;
    }
  }

  // slave file
  // header
  fprintf(_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n",
          title);

  multi_real_t pos;
  for (uint32_t y = 1; y < (_fieldsize[1] - 1); ++y) {
    pos[1] = (double)y * _h[1];
    for (uint32_t x = 1; x < (_fieldsize[0] - 1); ++x) {
      pos[0] = (double)x * _h[0];
      fprintf(_handle, "%le %le 0 ", v1->Interpolate(pos), v2->Interpolate(pos));
    }
    fprintf(_handle, "\n");
  }

  fprintf(_handle, "</DataArray>\n");

  // print info about data-segment into parallel master file
  if (_rank == 0) {
    fprintf(_masterHandle,
            "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" "
            "format=\"ascii\"/>\n",
            title);
  }
}

//------------------------------------------------------------------------------

void VTK::AddPointScalar(const char *title, const Grid *grid) {
  if (_rank == 0) {
    if (!_handle || !_masterHandle) {
      return;
    }
  } else {
    if (!_handle) {
      return;
    }
  }

  // slave file
  // header
  fprintf(_handle,
          "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\">\n", title);

  // grab payload and print values to file
  multi_real_t pos;
#if DIM == 3
  pos[2] = 0;
#endif // DIM
  for (uint32_t y = 0; y < (_fieldsize[1] - 1); ++y) {
    pos[1] = (double)y * _h[1];
    for (uint32_t x = 0; x < (_fieldsize[0] - 1); ++x) {
      pos[0] = (double)x * _h[0];
      fprintf(_handle, "%le ", (double)grid->Interpolate(pos));
    }
#if DIM == 2
    fprintf(_handle, "\n");
#endif // DIM
  }

  fprintf(_handle, "</DataArray>\n");

  // print info about data-segment into parallel master file
  if (_rank == 0) {
    fprintf(_masterHandle,
            "<PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\"/>\n",
            title);
  }
}

//------------------------------------------------------------------------------

void VTK::AddPointField(const char *title, const Grid *v1, const Grid *v2) {
  if (_rank == 0) {
    if (!_handle || !_masterHandle) {
      return;
    }
  } else {
    if (!_handle) {
      return;
    }
  }

  // slave file
  // header
  fprintf(_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n",
          title);

  // grab payload and print values to file
  multi_real_t pos;
#if DIM == 3
  pos[2] = 0;
#endif // DIM
  for (uint32_t y = 0; y < (_fieldsize[1] - 1); ++y) {
    pos[1] = (double)y * _h[1];
    for (uint32_t x = 0; x < (_fieldsize[0] - 1); ++x) {
      pos[0] = (double)x * _h[0];
      fprintf(_handle, "%le %le 0\n", (double)v1->Interpolate(pos),
              (double)v2->Interpolate(pos));
    }
  }

  // close xml tag
  fprintf(_handle, "</DataArray>\n");

  // master file
  // print info about DataField into masterfile
  if (_rank == 0) {
    fprintf(_masterHandle,
            "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" "
            "format=\"ascii\"/>\n",
            title);
  }
}

//------------------------------------------------------------------------------

void VTK::AddPointField(const char *title, const Grid *v1, const Grid *v2,
                        const Grid *v3) {
  if (!_handle)
    return;

  fprintf(_handle, "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\" "
                   "NumberOfComponents=\"3\">\n",
          title);

  multi_real_t pos;
#if DIM == 3
  pos[2] = 0;
#endif // DIM
  for (uint32_t y = 0; y < _fieldsize[1]; ++y) {
    pos[1] = (double)y * _h[1] + _offset[1];
    for (uint32_t x = 0; x < _fieldsize[0]; ++x) {
      pos[0] = (double)x * _h[0] + _offset[0];
      fprintf(_handle, "%le %le %le\n", (double)v1->Interpolate(pos),
              (double)v2->Interpolate(pos), (double)v3->Interpolate(pos));
    }
  }

  fprintf(_handle, "</DataArray>\n");
}

//------------------------------------------------------------------------------

void VTK::AddRank() {
  if (_rank == 0) {
    if (!_handle || !_masterHandle) {
      return;
    }
  } else {
    if (!_handle) {
      return;
    }
  }

  // slave file
  fprintf(_handle,
          "<DataArray Name=\"MPI-Rank\" type=\"Int32\" format=\"ascii\">\n");

  for (uint32_t y = 0; y < _fieldsize[1]; ++y) {
    for (uint32_t x = 0; x < _fieldsize[0]; ++x) {
      fprintf(_handle, "%i ", _rank);
    }
    fprintf(_handle, "\n");
  }
  fprintf(_handle, "</DataArray>\n");

  // master file
  if (_rank == 0) {
    fprintf(
        _masterHandle,
        "<PDataArray type=\"Int32\" Name=\"MPI-Rank\" format=\"ascii\"/>\n");
  }
}

//------------------------------------------------------------------------------
