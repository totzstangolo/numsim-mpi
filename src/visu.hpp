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

#include "typedef.hpp"
#include "grid.hpp"
#include "SDL2/SDL.h"
//------------------------------------------------------------------------------
#ifndef __VISU_HPP
#define __VISU_HPP
//------------------------------------------------------------------------------
const char States[][42] = {
    "Velocity", // 0
    "U",        // 1
    "V",        // 2
    "P",        // 3
    "",         // 4
    "",         // 5
    "",         // 6
    "",         // 7
    "",         // 8
    ""          // 9
};
//------------------------------------------------------------------------------
/// A grid rendering class using SDL2
class Renderer {
public:
  /// Constructs a Renderer
  Renderer(const multi_real_t &length, const multi_real_t &h);
  ~Renderer();

  /// Initializes the Renderer
  void Init(const index_t &width, const index_t &height, const int &idx = 0);

  /// Defines visible dimensions
  void SetSlice(const index_t &xdim, const index_t &ydim,
                const multi_real_t &origin);

  /// Checks and returns the window status
  int Check();

  /// Updates the window with a given grid (auto-scale)
  int Render(const Grid *grid);
  /// Updates the window with a given grid
  int Render(const Grid *grid, const real_t &min, const real_t &max);

  /// Sets the visability of the grid cells
  void ShowGrid(bool grid);

private:
  index_t _x;
  index_t _y;
  multi_real_t _orig;
  real_t _min;
  real_t _max;
  int _state;
  index_t _width;
  index_t _height;
  SDL_Window *_window;
  SDL_Surface *_screen;
  const multi_real_t &_length;
  const multi_real_t &_h;
  bool _grid;
  index_t _click_x;
  index_t _click_y;
  int _idx;

  static uint32_t _count;
};
//------------------------------------------------------------------------------
/*! \class	Renderer
 *	This class renders grids on a SDL2 GUI using the grids evaluate
 *function.
 *	It is visualizing 2D slices through the grid data with colors.
 */
/*! \fn	Renderer::Renderer (const multi_real_t& h, const multi_real_t& length)
 *	Constructs a Renderer with a given mesh width \p h and total size \p
 *length
 *	of the grid.
 *
 *	\param	h		mesh width to use
 *	\param	length	total size of the area
 */
/*!	\fn	void Renderer::Init (const index_t& width, const index_t&
 *height, const int& idx = 0)
 *	Initializes the SDL window and sets its size.
 *
 *	\param	width	The width of the window in pixels
 *	\param	height	The height of the window in pixels
 *	\param	idx		The Index of the window. If zero, no index is
 *printed in the title.
 */
/*!	\fn void Renderer::SetSlice (const index_t& xdim, const index_t& ydim,
 *const multi_real_t& origin)
 *	Defines witch dimensions to display on the window with \p xdim and \p
 *ydim.
 *	All other dimensions are evaluated at the positions defined in \p
 *origin.
 *
 *	\param	xdim	The dimensions index to display along the x-axis
 *	\param	ydim	The dimensions index to display along the y-axis
 *	\param	origin	Offset to all other dimensions
 */
/*! \fn	int Renderer::Check ()
 *	Returns the status of the SDL window encoded as integer.
 *
 *	\return		The encoded status value of the window.
 *		\arg	-1: The window was closed or the \e ESC key was pressed.
 *		\arg	0~9: The current status of the window. It can be changed
 *by pressing
 *					a number key.
 *		\arg	10: The enter key was pressed. This doesn't effect the
 *status.
 */
/*!	\fn int Renderer::Render (const Grid *grid, const multi_real_t& offset)
 *	Updates the SDL window and checks for pressed keys. The visualization is
 *	auto scaled based on the maximal and minal value in the grid.
 *
 *	\param grid		Pointer to the grid that should be visualized.
 *	\return			The status of the window.
 */
/*!	\fn int Renderer::Render (const Grid *grid, const multi_real_t& offset,
 *const real_t& min, const real_t& max)
 *	Updates the SDL window and checks for pressed keys. The visualization is
 *	scaled based on \p min and \p max. Grid values are cropped by these
 *limits.
 *
 *	\param grid		Pointer to the grid that should be visualized.
 *	\param min		The minimal value.
 *	\param max		The maximal value.
 *	\return			The status of the window.
 */

/*!	\fn void Renderer::ShowGrid (bool grid)
 *	Turns the visability of the grid cells with meshwidth \p _h on or off.
 *
 *	\param grid		Visability of the cells.
 */
//------------------------------------------------------------------------------
#endif // __VISU_HPP
