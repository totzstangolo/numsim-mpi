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
#include "iterator.hpp"
#include "geometry.hpp"

//////////////////
#include <cstdlib>
#include <iostream>
#include <cmath>
//////////////////

//------------------------------------------------------------------------------
/** Determines _tdim by determining largest divisor of N_processes, while
 *  being smaller than sqrt(N_processes).
 *  _tdim is then {N_processes/divisor, divisor} (except for N_processes = 2).
 *  Rank 0 broadcasts _tdim to the remaining processes.
 */
Communicator::Communicator(int *argc, char ***argv){
    MPI_Init(argc,argv);
    _tdim[0]=0;
    _tdim[1]=0;
    MPI_Comm_size(MPI_COMM_WORLD,&_size);
     MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
    int dim[2] = {1,1};
    int l_div=1;
    while(dim[1]<sqrt(_size)+1){
        if(_size % dim[1] == 0) l_div = dim[1];
        dim[1]++;
    }
    dim[1] = l_div;
    dim[0] = _size/dim[1];
    if(dim[1]>dim[0]){
        l_div = dim[1];
        dim[1] = dim[0];
        dim[0] = l_div;
    }
    _tdim[0]=dim[0];
    _tdim[1]=dim[1];

    _tidx = {_rank % _tdim[0], _rank / _tdim[0]};
}

  /** Communicator destructor; finalizes MPI Environment
   */
Communicator::~Communicator(){
    MPI_Finalize();
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


  /** Gets the sum of all values and distributes the result among all
   *  processes
   *
   * \param [in] val The data over which the sum is to be calculated
   */
real_t Communicator::gatherSum(const real_t &val) const{
    real_t sum;
    MPI_Allreduce(&val,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return sum;
}

  /** Finds the minimum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the minimum
   */
real_t Communicator::gatherMin(const real_t &val) const{
    real_t min;
    MPI_Allreduce(&val,&min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    return min;
}

  /** Finds the maximum of the values and distributes the result among
   *  all processes
   *
   * \param [in] val The data over which to find the maximum
   */
real_t Communicator::gatherMax(const real_t &val) const{
    real_t max;
    MPI_Allreduce(&val,&max,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    return max;
}

  /** Synchronizes ghost layer
   *
   * \param [in] grid  The values to sync
   */
void Communicator::copyBoundary(Grid *grid) const{
    copyLeftBoundary(grid);
    copyRightBoundary(grid);
    copyTopBoundary(grid);
    copyBottomBoundary(grid);
}

  /** Decide whether our left boundary is a domain boundary
   */
const bool Communicator::isLeft() const{
    bool isBound=false;
    if(ThreadIdx()[0]==0) isBound=true;
    return isBound;
}

  /** Decide whether our right boundary is a domain boundary
   */
const bool Communicator::isRight() const{
    bool isBound=false;
    if(ThreadIdx()[0]==ThreadDim()[0]-1) isBound=true;
    return isBound;
}

  /** Decide whether our top boundary is a domain boundary
   */
const bool Communicator::isTop() const{
    bool isBound=false;
    if(ThreadIdx()[1]==ThreadDim()[1]-1) isBound=true;
    return isBound;
}

  /** Decide whether our bottom boundary is a domain boundary
   */
const bool Communicator::isBottom() const{
    bool isBound=false;
    if(ThreadIdx()[1]==0) isBound=true;
    return isBound;
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
    index_t size = grid->getGeometry()->Size()[1];
    double buff[size];
    int counter = 0;
    if(!isRight() && !isLeft()){
	MPI_Status status;
        BoundaryIterator iterL(grid->getGeometry());
        iterL.SetBoundary(1);
        for(iterL.First();iterL.Valid();iterL.Next()){
            buff[counter] = grid->Cell(iterL.Right());
            counter++;
	}
	MPI_Sendrecv_replace(&buff, size, MPI_DOUBLE, _rank-1, 0, _rank+1, 0, MPI_COMM_WORLD,
                       &status);
	BoundaryIterator iterR(grid->getGeometry());
	iterR.SetBoundary(3);
    	counter=0;
        for(iterR.First();iterR.Valid();iterR.Next()){
           grid->Cell(iterR) = buff[counter];
           counter++;
        }
    } else if(!isRight()) {
        MPI_Recv(buff,size,MPI_DOUBLE,_rank+1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	BoundaryIterator iterR(grid->getGeometry());
	iterR.SetBoundary(3);
    	counter=0;
    	for(iterR.First();iterR.Valid();iterR.Next()){
           grid->Cell(iterR) = buff[counter];
           counter++;
    	}
    } else {
        BoundaryIterator iterL(grid->getGeometry());
        iterL.SetBoundary(1);
        for(iterL.First();iterL.Valid();iterL.Next()){
            buff[counter] = grid->Cell(iterL.Right());
            counter++;
        }
        MPI_Send(&buff,size,MPI_DOUBLE,_rank-1,0,MPI_COMM_WORLD);
    }
    return true;
}

  /** Function to sync ghost layer on right boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
bool Communicator::copyRightBoundary(Grid *grid) const{
    index_t size = grid->getGeometry()->Size()[1];
    double buff[size];
    int counter = 0;
    if(!isRight() && !isLeft()){
	MPI_Status status;
        BoundaryIterator iterR(grid->getGeometry());
        iterR.SetBoundary(3);
        for(iterR.First();iterR.Valid();iterR.Next()){
            buff[counter] = grid->Cell(iterR.Left());
            counter++;
        }
	MPI_Sendrecv_replace(&buff, size, MPI_DOUBLE, _rank+1, 1, _rank-1, 1, MPI_COMM_WORLD,
                       &status);
        BoundaryIterator iterL(grid->getGeometry());
        iterL.SetBoundary(1);
        counter=0;
        for(iterL.First();iterL.Valid();iterL.Next()){
            grid->Cell(iterL) = buff[counter];
            counter++;
        }
    } else if(!isLeft()) {
        MPI_Recv(buff,size,MPI_DOUBLE,_rank-1,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	BoundaryIterator iterL(grid->getGeometry());
	iterL.SetBoundary(1);
    	counter=0;
    	for(iterL.First();iterL.Valid();iterL.Next()){
           grid->Cell(iterL) = buff[counter];
           counter++;
    	}
    } else {
        BoundaryIterator iterL(grid->getGeometry());
        iterL.SetBoundary(3);
        for(iterL.First();iterL.Valid();iterL.Next()){
            buff[counter] = grid->Cell(iterL.Left());
            counter++;
        }
        MPI_Send(&buff,size,MPI_DOUBLE,_rank+1,1,MPI_COMM_WORLD);
    }
    return true;
}

  /** Function to sync ghost layer on top boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
bool Communicator::copyTopBoundary(Grid *grid) const{
    index_t size = grid->getGeometry()->Size()[0];
    double buff[size];
    int counter = 0;
   if(!isTop() && !isBottom()){
	MPI_Status status;
        BoundaryIterator iter(grid->getGeometry());
        iter.SetBoundary(2);
        for(iter.First();iter.Valid();iter.Next()){
            buff[counter] = grid->Cell(iter.Down());
            counter++;
        }
	MPI_Sendrecv_replace(&buff, size, MPI_DOUBLE,_rank + ThreadDim()[0], 2, _rank - ThreadDim()[0], 2, MPI_COMM_WORLD,
                       &status);
        iter.SetBoundary(0);
        counter=0;
        for(iter.First();iter.Valid();iter.Next()){
            grid->Cell(iter) = buff[counter];
            counter++;
        }
    } else if(!isBottom()) {
        MPI_Recv(&buff,size,MPI_DOUBLE,_rank-ThreadDim()[0],2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        BoundaryIterator iter(grid->getGeometry());
        iter.SetBoundary(0);
        counter=0;
        for(iter.First();iter.Valid();iter.Next()){
            grid->Cell(iter) = buff[counter];
            counter++;
        }
    } else if(!isTop()){
        BoundaryIterator iter(grid->getGeometry());
        iter.SetBoundary(2);
        for(iter.First();iter.Valid();iter.Next()){
            buff[counter] = grid->Cell(iter.Down());
            counter++;
        }
        MPI_Send(&buff,size,MPI_DOUBLE,_rank+ThreadDim()[0],2,MPI_COMM_WORLD);
    }
    return true;
}

  /** Function to sync ghost layer on bottom boundary
   *  Details analog to left boundary
   *
   * \param [in] grid  values whose boundary shall be synced
   */
bool Communicator::copyBottomBoundary(Grid *grid) const{
    index_t size = grid->getGeometry()->Size()[0];
    double buff[size];
    int counter = 0;
   if(!isTop() && !isBottom()){
	MPI_Status status;
        BoundaryIterator iter(grid->getGeometry());
        iter.SetBoundary(0);
        for(iter.First();iter.Valid();iter.Next()){
            buff[counter] = grid->Cell(iter.Top());
            counter++;
        }
	MPI_Sendrecv_replace(&buff, size, MPI_DOUBLE,_rank - ThreadDim()[0], 3, _rank + ThreadDim()[0], 3, MPI_COMM_WORLD,
                       &status);
        iter.SetBoundary(2);
        counter=0;
        for(iter.First();iter.Valid();iter.Next()){
            grid->Cell(iter) = buff[counter];
            counter++;
        }
    } else if(!isTop()) {
        MPI_Recv(&buff,size,MPI_DOUBLE,_rank+ThreadDim()[0],3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        BoundaryIterator iter(grid->getGeometry());
        iter.SetBoundary(2);
        counter=0;
        for(iter.First();iter.Valid();iter.Next()){
            grid->Cell(iter) = buff[counter];
            counter++;
        }
    } else if(!isBottom()){
        BoundaryIterator iter(grid->getGeometry());
        iter.SetBoundary(0);
        for(iter.First();iter.Valid();iter.Next()){
            buff[counter] = grid->Cell(iter.Top());
            counter++;
        }
        MPI_Send(&buff,size,MPI_DOUBLE,_rank-ThreadDim()[0],3,MPI_COMM_WORLD);
    }
    return true;
}
