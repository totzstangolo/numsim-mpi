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

//------------------------------------------------------------------------------
#include "grid.hpp"
#include "typedef.hpp"
//------------------------------------------------------------------------------
#ifndef __COMMUNICATOR_HPP
#define __COMMUNICATOR_HPP
//------------------------------------------------------------------------------
class Communicator {
public:
  /** Communicator constructor; initializes MPI Environment
   *
   * \param [in] argc Number of arguments program was started with
   * \param [in] argv Arguments passed to the program on start
   */
  Communicator(int *argc, char ***argv);

  /** Communicator destructor; finalizes MPI Environment
   */
  ~Communicator();

  /** Returns the position of the current process with respect to the
   *  fields lower left corner
   */
  const multi_index_t &ThreadIdx() const;

  /** Returns the way the domain is partitioned among all processes
   */
  const multi_index_t &ThreadDim() const;

  /** Returns whether this process is a red or a black field
   */
  const bool &EvenOdd() const;

  /** Gets the sum of all values and distributes the result among all
   *  processes
   *
   * \param [in] val The data over which the sum is to be calculated
   */
  real_t gatherSum(const real_t &val) const;

  /** Finds the minimum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the minimum
   */
  real_t gatherMin(const real_t &val) const;

  /** Finds the maximum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the maximum
   */
  real_t gatherMax(const real_t &val) const;

  /** Synchronizes ghost layer
   *
   * \param [in] grid  The values to sync
   */
  void copyBoundary(Grid *grid) const;

  /** Decide whether our left boundary is a domain boundary
   */
  const bool isLeft() const;

  /** Decide whether our right boundary is a domain boundary
   */
  const bool isRight() const;

  /** Decide whether our top boundary is a domain boundary
   */
  const bool isTop() const;

  /** Decide whether our bottom boundary is a domain boundary
   */
  const bool isBottom() const;

  /** Get MPI rank of current process
   */
  const int &getRank() const;

  /** Get number of MPI processes
   */
  const int &getSize() const;

private:
  multi_index_t _tidx;
  multi_index_t _tdim;
  int _rank;
  int _size;
  bool _evenodd;

  /** Function to sync ghost layer on left boundary:
   *  send values of own left boundary to left neighbor and
   *  and receive values from his right boundary
   *
   *   ------------ ------------
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *  |           x|y           |
   *   ------------ ------------
   *
   *   y: values that are sent
   *   x: values that are received
   *
   * \param [in] grid  values whose boundary shall be synced
   */
  bool copyLeftBoundary(Grid *grid) const;

  /** Function to sync ghost layer on right boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
  bool copyRightBoundary(Grid *grid) const;

  /** Function to sync ghost layer on top boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
  bool copyTopBoundary(Grid *grid) const;

  /** Function to sync ghost layer on bottom boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
  bool copyBottomBoundary(Grid *grid) const;
};
//------------------------------------------------------------------------------
#endif // __COMMUNICATOR_HPP
//------------------------------------------------------------------------------
