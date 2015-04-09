//
//  miAP.cpp
//  aracneGPU
//
//  Created by jing on 3/24/15.
//  Copyright (c) 2015 jing. All rights reserved.
//

#include "miAP.h"
using namespace std;

// function 1 sort GenePair 

typedef std::
int sortX < GenePair, GenePair, bool >
{
  // get X rank
  bool operator() ( const GenePair & lhs, const GenePair & rhs ) const
	  {
	    if ( lhs.Get_X() != rhs.Get_X() ) {
	      return ( lhs.Get_X() < rhs.Get_X() );
	    } else {
	      return ( lhs.Get_MaID() < rhs.Get_MaID() );
	    } 
	  }
  // not yet finish
}

// function 2 compute pairwise MIv
double calMI( GenePairVector &)
{
  // rank
  //
  // init 
   
}


float miAP(){
    // given X, Y as vectors, calculate Mutual information using adaptive partition
    // sort and get index
    // sort ( GenePair) in CPU
}
