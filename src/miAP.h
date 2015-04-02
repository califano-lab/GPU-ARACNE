//
//  miAP.h
//  aracneGPU
//
//  Created by jing on 3/24/15.
//  Copyright (c) 2015 jing. All rights reserved.
//

#ifndef __aracneGPU__miAP__
#define __aracneGPU__miAP__

#include <stdio.h>

class GenePair
{
  double x;
  double y;
  int xi;
  int yi;
  int maId;
public :
   inline double Get_X() const { return x; }
   inline double Get_Y() const { return y; }
   inline int    Get_XI() const { return xi; } 
   inline int    Get_YI() const { return yi; } 
   inline int  Get_MaID() const { return maId; } 
   inline void Set_X( double X ) { x = X; } 
   inline void Set_Y( double Y ) { y = Y; } 
   inline void Set_XI( int XI ) { xi = XI; } 
   inline void Set_YI( int YI ) { yi = YI; } 
   inline void Set_MaID( int MaID ) { maId = MaID; }

};


class MutualInfo 
{
  
};

#endif /* defined(__aracneGPU__miAP__) */
