#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

#include <iostream>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal()
{
/**** TODO ****/

	// the use of SetNormal is way better than to set it manually.
	myPoint3D* a = adjacent_halfedge->source->point;
	myPoint3D* b = adjacent_halfedge->next->source->point;
	myPoint3D* c = adjacent_halfedge->next->next->source->point;

	//myVector3D vA = myVector3D(b->source->point->X - a->source->point->X, b->source->point->Y - a->source->point->Y, b->source->point->Z - a->source->point->Z);
	//myVector3D vB = myVector3D(c->source->point->X - b->source->point->X, c->source->point->Y - b->source->point->Y, c->source->point->Z - b->source->point->Z);

	//myVector3D toCompute;
	//toCompute = vA.crossproduct(vB);
	//std::cout << " Normals  X: " << toCompute.dX << " Y : " << toCompute.dY << " Z: " << toCompute.dZ << std::endl;
	//normal = &toCompute;

	normal->setNormal(a,b,c);
}

bool myFace::isTriangle()
{

	return (adjacent_halfedge != NULL) && (adjacent_halfedge->next->next->next == adjacent_halfedge);
}

bool myFace::isQuad()
{

	return (adjacent_halfedge != NULL) && (adjacent_halfedge->next->next->next->next == adjacent_halfedge);
}

int myFace::getSize()
{
	int i = 0; 
	myHalfedge* e = adjacent_halfedge;
	do
	{
		i++;
		e = e->next;
	} while (e != adjacent_halfedge);

		return i;
}