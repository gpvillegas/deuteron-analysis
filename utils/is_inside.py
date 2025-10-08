#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:05:24 2025

@author: gvill

  !**********************************************************************
  !                                                               
  !                                                               
  !     SUBPROGRAM DESCRIPTION:                                   
  !          is_inside determines whether point (x,y) is  
  !          inside or outside of the boundary given by the nlim
  !          number of points xlim, ylim                          
  !                                                                  
  !                                                                  
  !     CALLING ARGUMENTS:                                           
  !       x...............x-coordinate of point (input)                            
  !       y...............y-coordinate of point (input)                            
  !       nlim............number of boundary points (input)           
  !       xlim............x coordinates of boundary (input)           
  !       ylim............y coordinates of boundar  (input)           
  !                                                                  
  !     REFERENCES:                                                  
  !          (1)  uses ray casting
  !                                                                  
  !     RECORD OF MODIFICATION:                                      
  !          created 06/03/2020... W. Boeglin                  
  !                                                                  
  !                                                                  
  !                                                                  
  !**********************************************************************
  
  This is a python version of Werner's subroutine: is_inside.
  
  Function is_inside:
      inputs:
          n: int, size of x/y arrays
          x: array, x-coordinates
          y: array, y-coordinates
          nlim: int, number of boundary points
          xlim: array, x-coordinates of boundary
          ylim: array, y-coordinates of boundary
          
      returns:
          is_point_inside: bool, result of evaluating if point (x,y) is inside
                              the boundary determined by the points (xlim,ylim)
  
"""
import numpy as np
import LT.box as B

def is_inside(n, x, y, nlim, xlim, ylim):
    is_point_inside = []
    for i in range(n):
        # trace a horizontal ray from the point (x,y) and count the number
        # of times it intersects with a boundary edge. If odd, point is 
        # inside, if even, point is outside
        ncross = 0
        for j in range(nlim-1):
            if (ylim[j] < y[i]) and (ylim[j+1] < y[i]): continue
            # eliminates points above and below the polygon   
            # if (y[i] > max(ylim)) or (y[i] < min(ylim)): continue
            
            # this condition works by eliminating points along vertical lines
            # on the x boundary points
            if (x[i] == xlim[j]): continue 
        
            t = x[i] - xlim[j]
            s = xlim[j+1] - x[i]
            if (t*s < 0.): continue
        
            di = (ylim[j+1] - ylim[j])/(xlim[j+1] - xlim[j])
            f = ylim[j] + di*(x[i] - xlim[j])
            if (f < y[i]): continue
        
            ncross += 1
    
        is_point_inside.append((ncross % 2) == 1)    
    return np.array(is_point_inside)

#%% example
# array off polygon vertex points
# pa = np.array( [[0,0], [.4,2], [.2, 6.], [0,2.], [0,0]] )
tri = np.array([(0,0),(0.4,2),(0.5,0),(0,0)])
sq = np.array([(0,0),(0,2),(0.4,2),(0.4,0),(0,0)])
hexa = np.array([(0,0),(-0.05,1),(-0.05,2),(0,3),(0.2,3),(0.25,2),(0.25,1),
                 (0.2,0),(0,0)])

# create som psuedo points
xp = np.random.normal(loc = 0.1, size = 10000, scale = 0.15)
yp = np.random.normal(loc = 1.5, size = 10000, scale = 2.)

# plot the points and the polygon
B.pl.figure()
B.pl.plot(xp, yp, '.')
# B.pl.plot(pa[:,0], pa[:,1])
# B.pl.plot(tri[:,0], tri[:,1])
# B.pl.plot(sq[:,0], sq[:,1])
B.pl.plot(hexa[:,0], hexa[:,1])

# check which ponts are inside
# res = is_inside(xp.size, xp, yp, pa[:,0].size, pa[:,0], pa[:, 1])
# res = is_inside(xp.size, xp, yp, tri[:,0].size, tri[:,0], tri[:, 1])
# res = is_inside(xp.size, xp, yp, sq[:,0].size, sq[:,0], sq[:, 1])
res = is_inside(xp.size, xp, yp, hexa[:,0].size, hexa[:,0], hexa[:, 1])

# plot the points inside again 
B.pl.plot(xp[res], yp[res], '.')
