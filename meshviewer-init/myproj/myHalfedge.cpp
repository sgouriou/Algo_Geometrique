#include "myHalfedge.h"

int myHalfedge::indexHE = 0;

myHalfedge::myHalfedge(void)
{
	source = NULL; 
	adjacent_face = NULL; 
	next = NULL;  
	prev = NULL;  
	twin = NULL;  
	index = myHalfedge::indexHE++;
}

void myHalfedge::copy(myHalfedge *ie)
{
/**** TODO ****/
	//myHalfedge toReturn;
	//toReturn.source = source;

	//return toReturn;

}

myHalfedge::~myHalfedge(void)
{
}
