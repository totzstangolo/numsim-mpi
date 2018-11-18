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

#include "visu.hpp"
#include <cmath>
#include <limits>
//------------------------------------------------------------------------------
uint8_t HueR (real_t value, real_t min, real_t max) {
  if (value < min) value = min;
  if (value > max) value = max;
  double dval = max - min;
  
  if (value < (min + 0.25 * dval))
    return 0;
  if (value < (min + 0.5  * dval))
    return 0;
  if (value < (min + 0.75 * dval))
    return  255*(4 * (value - min - 0.5 * dval) / dval);
  return 255;
}
//------------------------------------------------------------------------------
uint8_t HueG (real_t value, real_t min, real_t max) {
  if (value < min) value = min;
  if (value > max) value = max;
  double dval = max - min;
  
  if (value < (min + 0.25 * dval))
    return 255*(4 * (value - min) / dval);
  if (value < (min + 0.5  * dval))
    return 255;
  if (value < (min + 0.75 * dval))
    return  255;
  return 255*(1 + 4 * (min + 0.75 * dval - value) / dval);
}
//------------------------------------------------------------------------------
uint8_t HueB (real_t value, real_t min, real_t max) {
  if (value < min) value = min;
  if (value > max) value = max;
  double dval = max - min;
  
  if (value < (min + 0.25 * dval))
    return 255;
  if (value < (min + 0.5  * dval))
    return 255*(1 + 4 * (min + 0.25 * dval - value) / dval);;
  if (value < (min + 0.75 * dval))
    return  0;
  return 0;
}
//------------------------------------------------------------------------------
void setpixelhue(SDL_Surface *screen, int x, int y, real_t value, real_t min,
                 real_t max) {
  uint32_t *pixmem32;
  uint32_t colour;

  colour = SDL_MapRGB(screen->format, HueR(value, min, max),
                      HueG(value, min, max), HueB(value, min, max));

  pixmem32 = (uint32_t *)screen->pixels + y * screen->w + x;
  *pixmem32 = colour;
}
//------------------------------------------------------------------------------
void setpixelrgb(SDL_Surface *screen, int x, int y, uint8_t r, uint8_t g,
                 uint8_t b) {
  uint32_t *pixmem32;
  uint32_t colour;

  colour = SDL_MapRGB(screen->format, r, g, b);

  pixmem32 = (uint32_t *)screen->pixels + y * screen->w + x;
  *pixmem32 = colour;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
uint32_t Renderer::_count = 0;
//------------------------------------------------------------------------------
Renderer::Renderer(const multi_real_t &length, const multi_real_t &h)
    : _length(length), _h(h) {
  if (_count == 0)
    SDL_Init(SDL_INIT_VIDEO);
  _count++;
  _state = 0;
  for (uint32_t i = 0; i < DIM; ++i)
    _orig[i] = 0;
  _x = 0;
  _y = 1;
  _click_x = 0;
  _click_y = 1;
  _grid = true;
  _min = std::numeric_limits<real_t>::max();
  _max = std::numeric_limits<real_t>::min();
}
//------------------------------------------------------------------------------
Renderer::~Renderer() {
  if (_window)
    SDL_DestroyWindow(_window);
  if (_count == 1)
    SDL_Quit();
  _count--;
}
//------------------------------------------------------------------------------
void Renderer::Init(const index_t &width, const index_t &height,
                    const int &idx) {
  // SDL_WM_SetCaption("Grid Renderer","Grid Renderer");
  _width = width;
  _height = height;
  _idx = idx;
  _window = SDL_CreateWindow("", SDL_WINDOWPOS_UNDEFINED,
                             SDL_WINDOWPOS_UNDEFINED, _width, _height, 0);
  _screen = SDL_GetWindowSurface(_window); // SDL_SetVideoMode(_width,_height,
                                           // 32, SDL_HWSURFACE |
                                           // SDL_DOUBLEBUF);
}
//------------------------------------------------------------------------------
void Renderer::SetSlice(const index_t &xdim, const index_t &ydim,
                        const multi_real_t &origin) {
  _x = xdim;
  _y = ydim;
  for (uint32_t i = 0; i < DIM; ++i)
    _orig[i] = origin[i];
}
//------------------------------------------------------------------------------
int Renderer::Check() {
  SDL_Event event;
  while (SDL_PollEvent(&event)) {
    switch (event.type) {
    case SDL_QUIT:
      return -1;
    case SDL_WINDOWEVENT:
      if (event.window.event == SDL_WINDOWEVENT_CLOSE) {
        _state = -1;
        SDL_DestroyWindow(_window);
        _window = NULL;
      }
      break;
    case SDL_MOUSEBUTTONDOWN:
      _click_x = event.button.x;
      _click_y = event.button.y;
      break;
    case SDL_KEYDOWN:
      switch (event.key.keysym.sym) {
      case SDLK_ESCAPE:
        return -1;
      case SDLK_0:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 0;
        break;
      case SDLK_1:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 1;
        break;
      case SDLK_2:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 2;
        break;
      case SDLK_3:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 3;
        break;
      case SDLK_4:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 4;
        break;
      case SDLK_5:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 5;
        break;
      case SDLK_6:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 6;
        break;
      case SDLK_7:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 7;
        break;
      case SDLK_8:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 8;
        break;
      case SDLK_9:
        _min = std::numeric_limits<real_t>::max();
        _max = std::numeric_limits<real_t>::min();
        _state = 9;
        break;
      case SDLK_RETURN:
        return 10;
      default:
        break;
      };
      break;
    default:
      break;
    };
  }
  return _state;
}
//------------------------------------------------------------------------------
int Renderer::Render(const Grid *grid) { return Render(grid, _min, _max); }
//------------------------------------------------------------------------------
int Renderer::Render(const Grid *grid, const real_t &min, const real_t &max) {
  if (Check() < 0)
    return -1;
  real_t treshold[2] = {0, 0};
  real_t value;
  if (SDL_MUSTLOCK(_screen))
    if (SDL_LockSurface(_screen) < 0)
      return -1;

  _orig[_x] = _length[_x] * _click_x / _width;
  _orig[_y] = _length[_y] * (_height - _click_y - 1) / _height;
  char title[200];
  if (_idx != 0)
    sprintf(title, "Grid Renderer %s %i - grid(%lf", States[_state % 10], _idx,
            _orig[0]);
  else
    sprintf(title, "Grid Renderer %s - grid(%lf", States[_state % 10],
            _orig[0]);
  for (uint32_t d = 1; d < DIM; ++d)
    sprintf(title, "%s,%lf", title, _orig[d]);
  sprintf(title, "%s) = %le", title, grid->Interpolate(_orig));
  SDL_SetWindowTitle(_window, title);

  for (uint32_t x = 0; x < _width; ++x) {
    treshold[1] = _length[_y];
    _orig[_x] = _length[_x] * x / _width;
    for (uint32_t y = 0; y < _height; ++y) {
      _orig[_y] = _length[_y] * (_height - y - 1) / _height;
      if (x == _click_x || y == _click_y) {
        setpixelrgb(_screen, x, y, 255, 255, 255);
        continue;
      }
      if (_grid && _orig[_x] >= treshold[0]) {
        setpixelrgb(_screen, x, y, 20, 20, 20);
      } else if (_grid && _orig[_y] < treshold[1]) {
        setpixelrgb(_screen, x, y, 20, 20, 20);
        treshold[1] -= _h[_y];
      } else {
        value = grid->Interpolate(_orig);
        setpixelhue(_screen, x, y, value, min, max);
        if (value > _max)
          _max = value;
        if (value < _min)
          _min = value;
      }
    }
    if (_grid && _orig[_x] >= treshold[0]) {
      treshold[0] += _h[_x];
    }
  }
  if (SDL_MUSTLOCK(_screen))
    SDL_UnlockSurface(_screen);
  SDL_UpdateWindowSurface(_window);
  // SDL_Flip(_screen);
  return _state;
}
//------------------------------------------------------------------------------
void Renderer::ShowGrid(bool grid) { _grid = grid; }
//------------------------------------------------------------------------------
