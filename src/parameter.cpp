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

#include "parameter.hpp"
#include <cstdio>
#include <cstring>
#include <cstdlib>
//------------------------------------------------------------------------------
Parameter::Parameter() {
  _re = 1000;
  _omega = 1.7;
  _alpha = 0.9;
  _dt = 0.0;
  _tend = 50;
  _itermax = 100;
  _eps = 0.001;
  _tau = 0.5;
}
//------------------------------------------------------------------------------
void Parameter::Load(const char *file) {
  FILE *handle = fopen(file, "r");
  double inval;
  char name[20];
  while (!feof(handle)) {
    if (!fscanf(handle, "%s = %lf\n", name, &inval))
      continue;
    if (strcmp(name, "re") == 0)
      _re = inval;
    else if (strcmp(name, "omg") == 0)
      _omega = inval;
    else if (strcmp(name, "alpha") == 0)
      _alpha = inval;
    else if (strcmp(name, "dt") == 0)
      _dt = inval;
    else if (strcmp(name, "tend") == 0)
      _tend = inval;
    else if (strcmp(name, "iter") == 0)
      _itermax = inval;
    else if (strcmp(name, "eps") == 0)
      _eps = inval;
    else if (strcmp(name, "tau") == 0)
      _tau = inval;
    else
      printf("Unknown parameter %s\n", name);
  }
  fclose(handle);
}
//------------------------------------------------------------------------------
const real_t &Parameter::Re() const { return _re; }
//------------------------------------------------------------------------------
const real_t &Parameter::Omega() const { return _omega; }
//------------------------------------------------------------------------------
const real_t &Parameter::Alpha() const { return _alpha; }
//------------------------------------------------------------------------------
const real_t &Parameter::Dt() const { return _dt; }
//------------------------------------------------------------------------------
const real_t &Parameter::Tend() const { return _tend; }
//------------------------------------------------------------------------------
const index_t &Parameter::IterMax() const { return _itermax; }
//------------------------------------------------------------------------------
const real_t &Parameter::Eps() const { return _eps; }
//------------------------------------------------------------------------------
const real_t &Parameter::Tau() const { return _tau; }
//------------------------------------------------------------------------------
