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
#include "communicator.hpp"

//////////////////
#include <cstdlib>
#include <iostream>
#include <cmath>
//////////////////

//------------------------------------------------------------------------------
Communicator::Communicator(int *argc, char ***argv){
    MPI_Init(argc,argv);
    _tdim[0]=0;
    _tdim[1]=0;
    MPI_Comm_size(MPI_COMM_WORLD,&_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
    MPI_Request reqs[2*(_size-1)];
    MPI_Status stats[2*(_size-1)];
    int recvBuf[2];
    if(_rank!=0){
        MPI_Recv(&_tdim,2,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }else{
        int xDim=1;
        int yDim=1;
        int l_div=1;
        while(yDim<sqrt(_size)+1){
            if(_size % yDim == 0) l_div = yDim;
            yDim++;
        }
        yDim = l_div;
        xDim = _size/yDim;
        if(yDim>xDim){
            l_div = yDim;
            yDim = xDim;
            xDim = l_div;
        }
        _tdim[0]=xDim;
        _tdim[1]=yDim;
        for(int iCount=1;iCount<_size;iCount++){
            MPI_Send(&_tdim,2,MPI_INT,iCount,0,MPI_COMM_WORLD);
        }
    }
    std::cout << "I know about decomposition ("<<_tdim[0]<<","<<_tdim[1]<< ")!"
    << std::endl;
    MPI_Finalize();
}

  /** Communicator destructor; finalizes MPI Environment
   */
Communicator::~Communicator(){

}

  /** Returns the position of the current process with respect to the
   *  fields lower left corner
   */
const multi_index_t &Communicator::ThreadIdx() const{
    return _tidx;
}

  /** Returns the way the domain is partitioned among all processes
   */
const multi_index_t &Communicator::ThreadDim() const{
    return _tdim;
}

  /** Returns whether this process is a red or a black field
   */
const bool &Communicator::EvenOdd() const{
    return true;
}

  /** Gets the sum of all values and distributes the result among all
   *  processes
   *
   * \param [in] val The data over which the sum is to be calculated
   */
real_t Communicator::gatherSum(const real_t &val) const{
    real_t dummy = 1.0;
    return dummy;
}

  /** Finds the minimum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the minimum
   */
real_t Communicator::gatherMin(const real_t &val) const{
    real_t dummy = 1.0;
    return dummy;
}

  /** Finds the maximum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the maximum
   */
real_t Communicator::gatherMax(const real_t &val) const{
    real_t dummy = 1.0;
    return dummy;
}

  /** Synchronizes ghost layer
   *
   * \param [in] grid  The values to sync
   */
void Communicator::copyBoundary(Grid *grid) const{

}

  /** Decide whether our left boundary is a domain boundary
   */
const bool Communicator::isLeft() const{
    return true;
}

  /** Decide whether our right boundary is a domain boundary
   */
const bool Communicator::isRight() const{
    return true;
}

  /** Decide whether our top boundary is a domain boundary
   */
const bool Communicator::isTop() const{
    return true;
}

  /** Decide whether our bottom boundary is a domain boundary
   */
const bool Communicator::isBottom() const{
    return true;
}

  /** Get MPI rank of current process
   */
const int &Communicator::getRank() const{
    return _rank;
}

  /** Get number of MPI processes
   */
const int &Communicator::getSize() const{
    return _size;
}

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
bool Communicator::copyLeftBoundary(Grid *grid) const{
    return true;
}

  /** Function to sync ghost layer on right boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
bool Communicator::copyRightBoundary(Grid *grid) const{
    return true;
}

  /** Function to sync ghost layer on top boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
bool Communicator::copyTopBoundary(Grid *grid) const{
    return true;
}

  /** Function to sync ghost layer on bottom boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
bool Communicator::copyBottomBoundary(Grid *grid) const{
    return true;
}
