#pragma once

class myHalfedge;
class myVector3D;
class myPoint3D;

class myFace
{
public:
	myHalfedge *adjacent_halfedge;

	myVector3D *normal;

	myPoint3D *centroid;

	int index; //use this variable as you wish.
	static int indexF;

	void computeNormal();
	myFace(void);
	~myFace(void);

	bool isTriangle();
	bool isQuad();
	int getSize();
	

};