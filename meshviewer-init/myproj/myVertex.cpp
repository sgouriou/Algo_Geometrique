#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"

#include <iostream>

myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,1.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal()
{
	/**** TODO ****/



	myVector3D* tempNormal = new myVector3D(0,0,0);

	myHalfedge* tempHE = originof;

	if (tempHE == NULL)
		return;

	do {
		*tempNormal += *(tempHE->adjacent_face->normal);


		tempHE = tempHE->twin->next;



	} while(tempHE != originof);

	tempNormal->normalize();

	*normal = *tempNormal;


}
