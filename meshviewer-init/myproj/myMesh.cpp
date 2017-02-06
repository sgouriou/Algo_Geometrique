#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"
#include <algorithm>



using namespace std;


void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}


myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	bool error;
	unsigned int count;

	cout << "Checking mesh for errors...\n";
	cout << "\tNumber of vertices: " << vertices.size() << endl;
	cout << "\tNumber of halfedges: " << halfedges.size() << endl;
	cout << "\tNumber of faces: " << faces.size() << endl << endl;

	cout << "\tChecking for NULL vertices. " << endl;
	for (unsigned int i = 0; i < vertices.size(); i++)
	{
		if (vertices[i] == NULL)
			cout << "\t\tError: vertex " << i << " is NULL.\n";
		else if (vertices[i]->originof == NULL)
			cout << "\t\tError: the originof halfedge of vertex " << i << " is NULL.\n";
	}
	cout << "\t  Ended check.\n\n";

	cout << "\tChecking the halfedges. " << endl;
	error = false;
	count = 0;
	for (unsigned int i = 0; i < halfedges.size(); i++)
	{
		if (halfedges[i]->source == NULL)
			cout << "\tError: Source is NULL for halfedge " << i << endl;
		if (halfedges[i]->twin == NULL) count++;
		if (halfedges[i]->next == NULL || halfedges[i]->prev == NULL)
			cout << "\tError: Next/prev NULL for halfedge " << i << endl;
		if (halfedges[i]->next->prev != halfedges[i] || halfedges[i]->prev->next != halfedges[i])
			cout << "\tError: Next/prev not set properly for halfedge " << i << endl;
		if (halfedges[i]->twin != NULL && halfedges[i] != halfedges[i]->twin->twin)
			cout << "\tError: Twin pair not set properly for halfedge " << i << endl;
	}
	if (count > 0) cout << "\tThis mesh has boundary edges.\n";
	cout << "\t  Ended check.\n\n";

	cout << "\tChecking fans of each vertex.\n";
	for (unsigned int i = 0; i<vertices.size(); i++) {
		myVertex *v = vertices[i];

		myHalfedge *e1 = v->originof;

		unsigned int k = 0;
		do {
			k++;
			e1 = e1->prev->twin;
			if (k > 100000) cout << "\t\tError: Infinite loop when checking adjacent edges for vertex " << i << endl;
		} while (e1 != NULL && e1 != v->originof);
	}
	cout << "\t  Ended check.\n\n";

	cout << "\tChecking edges of each face.\n";
	unsigned int num_incidentedgesoverallfaces = 0;
	bool istriangular = true;
	for (unsigned i = 0; i<faces.size(); i++) {
		myHalfedge *e1 = faces[i]->adjacent_halfedge;
		unsigned int k = 0;
		do {
			k++;
			if (e1 == NULL) cout << "\t\tError: Found NULL edge on boundary of face " << i << endl;
			e1 = e1->next;
			if (k > 100000) cout << "\t\tError: Infinite loop when checking adjacent edges for face " << i << endl;
		} while (e1 != faces[i]->adjacent_halfedge);
		num_incidentedgesoverallfaces += k;
		if (k>3) istriangular = false;
	}
	if (istriangular) cout << "\t\tThe mesh is triangular.\n";
	else cout << "\t\tThe mesh is not triangular.\n";
	if (num_incidentedgesoverallfaces != halfedges.size())
		cout << "\t\tSuspicious: the total number of halfedges is not equal to the sum at each face.\n";
	cout << "\t  Ended check.\n\n";

	cout << "  Ended check of mesh.\n";
}


bool myMesh::readFile(std::string filename)
{

	cout << "hello there \n";
 
	string s, t, u;
	vector<int> faceids;
	myHalfedge **hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

	while (getline(fin, s))
	{
		cout << "got a line here\n";
		stringstream myline(s);
		myline >> t;
		cout << "[" << t << "]" << "\n";
		if (t == "g") {}

		// vertex Case
		else if (t == "v")
		{

			cout << "got a vertex here\n";
			double x = 0.0;
			myline >> u;
			x = stof(u);

			double y = 0.0;
			myline >> u;
			y = stof(u);

			double z = 0.0;
			myline >> u;
			z = stof(u);


			myVertex *v = new myVertex();
			v->point = new myPoint3D(x,y,z);
			cout << "added x : " << x << " y:" << y << " z : " << z << "\n";
			vertices.push_back(v);

		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}

		//Faces 
		else if (t == "f")
		{
			//cout << "Face indices: "; 
			//while (myline >> u) cout << atoi((u.substr(0, u.find("/"))).c_str()) << " ";
			//cout << endl;

			myFace* currentFace = new myFace();
			
			
			vector<int> vertexIndicesofFace;
			vector<myHalfedge*> tempHalfedges;
			vector<string> separated;
			
			string tempStr;
			while(myline >> tempStr)
			{
				separated = split(tempStr, '/');
				cout << "face sep : " <<  separated[0] << endl;
				vertexIndicesofFace.push_back(stoi(separated[0]) - 1);
				tempHalfedges.push_back(new myHalfedge());
			}
			
			for (int i = 0; i < (tempHalfedges.size()); i++)
			{
				//myHalfedge* a = new myHalfedge();

				int next = (i + 1) % tempHalfedges.size();
				int prev = (i - 1 + tempHalfedges.size()) % tempHalfedges.size();


				//Initializing the connection between Vertex and Hedge.
				tempHalfedges[i]->source = vertices[vertexIndicesofFace[i]];
				vertices[vertexIndicesofFace[i]]->originof = tempHalfedges[i];

				tempHalfedges[i]->next = tempHalfedges[next];
				tempHalfedges[i]->prev = tempHalfedges[prev];
				tempHalfedges[i]->adjacent_face = currentFace;


				// pushing new halfEgde
				halfedges.push_back(tempHalfedges[i]);

				// now getting to the mapping part
				twin_map[make_pair(vertexIndicesofFace[i], vertexIndicesofFace[next])] = tempHalfedges[i];
				it = twin_map.find(make_pair(vertexIndicesofFace[next], vertexIndicesofFace[i]));

				if (it != twin_map.end())//if we found a twin candidate
				{
					it->second->twin = tempHalfedges[i];
					tempHalfedges[i]->twin = it->second;

				}
			}


				currentFace->adjacent_halfedge = tempHalfedges[0];
				faces.push_back(currentFace);


			


		}


	}

	checkMesh();
	normalize();

	return true;
}


void myMesh::computeNormals()
{
	/**** TODO ****/

	for (vector<myFace *>::iterator it = faces.begin(); it != faces.end(); it++)
		(*it)->computeNormal();
	for (vector<myVertex *>::iterator it = vertices.begin(); it != vertices.end(); it++)
		(*it)->computeNormal();

	std::cout << "end Compute Normals :" << std::endl;
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}


void myMesh::splitFaceTRIS(myFace *f, myPoint3D *p)
{
	checkFacesEdgeSize();


	myFace* newFace = new myFace();
	myVertex* vp = new myVertex();
	vp->point = p;

	myHalfedge* adj = new myHalfedge();
	adj = f->adjacent_halfedge;

	myHalfedge* triBefore = new myHalfedge();
	myHalfedge* triBeforeT = new myHalfedge();

	myHalfedge* triAfter = new myHalfedge();
	myHalfedge* triAfterT = new myHalfedge();

	triBefore->adjacent_face = newFace;
	triBeforeT->adjacent_face = f;


	triBefore->twin = triBeforeT;
	triBeforeT->twin = triBefore;

	triBefore->source = vp;
	triBeforeT->source = adj->source;

	triBefore->next = adj;
	adj->prev->next = triBeforeT;
	

	triBeforeT->prev = adj->prev;
	adj->prev = triBefore;
	//triBeforeT->prev->next = triBeforeT;
	
	// TriBefore DONE

	triAfter->adjacent_face = newFace;
	triAfterT->adjacent_face = f;

	triAfter->twin = triAfterT;
	triAfterT->twin = triAfter;

	triAfter->source = adj->next->source;
	triAfterT->source = vp;

	triAfter->next = triBefore;
	triBefore->prev = triAfter;

	triAfterT->next = adj->next;
	triAfterT->next->prev = triAfterT;

	triAfter->prev = adj;
	adj->next = triAfter;

	triAfterT->prev = triBeforeT;
	triBeforeT->next = triAfterT;



	adj->adjacent_face = newFace;
	f->adjacent_halfedge = triAfterT;

	newFace->adjacent_halfedge = adj;
	
	/// VERTICES

	vp->originof = triBefore;


	
	// Tri After DONE
	newFace->computeNormal();
	
	
	vertices.push_back(vp);
	halfedges.push_back(triBefore);
	halfedges.push_back(triBeforeT);

	halfedges.push_back(triAfter);
	halfedges.push_back(triAfterT);

	faces.push_back(newFace);

	checkFacesEdgeSize();

	cout << "SPLIT FACES : " << faces.size() << endl;
	cout << "STARTING TRIANGULATE" << endl;
	triangulate(f);
		
	cout << "SPLIT FACES : " << faces.size() << endl;







}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{

	myHalfedge* e2 = new myHalfedge();
	myHalfedge* g2 = new myHalfedge();

	myVertex* vp = new myVertex();
	vp->point = p;
	


	e2->source = vp;
	e2->adjacent_face = (e1->adjacent_face);


	e2->twin = (e1->twin);

	e2->next = e1->next;
	e2->prev = e1;

	e2->next->prev = e2;



	e1->next = e2;

	g2->source = vp;
	g2->twin = e1;
	
	g2->adjacent_face = (e1->twin->adjacent_face);
	g2->next = (e1->twin->next);
	g2->prev = (e1->twin);
	g2->next->prev = g2;
	e1->twin->next = g2;


	e1->twin = g2;

	g2->prev->twin = e2;
	vp->originof = e2;
	vp->originof = g2;

	vertices.push_back(vp);
	halfedges.push_back(e2);
	halfedges.push_back(g2);

	cout << "Got HERE\n";

	//triangulate(e1->adjacent_face);
	//triangulate(e1->twin->adjacent_face);
	


	//myVertex* intermediate = new myVertex();
	//intermediate->point = p;



	//index = 49.3;






	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	

	/**** TODO ****/
	if (f->isTriangle() || /*f->isQuad() || */(f->getSize() % 2 == 1))
	{
		// face is not compatible with splitFaceQuads
		return;
	}

	


	else
	{
		int nbEdges = f->getSize()/2;
		int nbFacesToCreate = nbEdges - 1;

		std::cout << "nbEdges : " << nbEdges << std::endl;
		std::cout << "nbfacesToCreate : " << nbFacesToCreate << std::endl;

		myVertex* vp = new myVertex();
		vp->point = p;

		vector<myHalfedge*> in = vector<myHalfedge*>();
		vector<myHalfedge*> out = vector<myHalfedge*>();


		// inS and outS are the opposite edge of in and out HE on the same face
		// IE : inS->outS->in->out
		vector<myHalfedge*> inS = vector<myHalfedge*>();
		vector<myHalfedge*> outS = vector<myHalfedge*>();

		vector<myFace*> nf = vector<myFace*>();

		for (int i = 0; i < nbEdges; i++)
		{
			in.push_back(new myHalfedge());
			out.push_back(new myHalfedge());

		}

		myHalfedge* temp = f->adjacent_halfedge;
		for (int i = 0; i < nbEdges; i++)
		{


			inS.push_back(temp);
			temp = temp->next;
			outS.push_back(temp);
			temp = temp->next;
		}


		// dont forget to start at 1 when we push into the real one
		nf.push_back(f);
		for (int i = 0; i < nbFacesToCreate; i++)
		{
			nf.push_back(new myFace());
		}

		for (int i = 0; i < nbEdges; i++)
		{

			inS[i]->adjacent_face = nf[i];
			inS[i]->prev = out[i];

			outS[i]->adjacent_face = nf[i];
			//Replaced furtuer
			//outS[i]->next = in[i];
			


			out[i]->source = vp;
			out[i]->next = inS[i];
			out[i]->prev = in[i];
			out[i]->adjacent_face = nf[i];

			//out[i]->next->prev = out[i];
			//out[i]->prev->next = out[i];

			in[i]->source = outS[i]->next->source;

			//replaced here
			 outS[i]->next = in[i];


			in[i]->next = out[i];
			in[i]->prev = outS[i];
			in[i]->adjacent_face = nf[i];


			//added recently
			nf[i]->adjacent_halfedge = out[i];
			vp->originof = out[i];
			//in[i]->next->prev = in[i];
			//in[i]->prev->next = in[i];


			//this offset will complete itself at the last call don't worry
			in[i]->twin = out[(i + 1) % nbEdges];
			out[(i + 1) % nbEdges]->twin = in[i];

		}

		vertices.push_back(vp);
		for (int i = 0; i < nbFacesToCreate; i++)
		{
			// 1+i car nf[0] représente la face initiale
			faces.push_back(nf[1 + i]);
		}

		for (int i = 0; i < nbEdges; i++)
		{
			halfedges.push_back(in[i]);
			halfedges.push_back(out[i]);

			

		}



	}


}


void myMesh::computeCentroids()
{
	//vector<myVertex*> myCentroids;
	myHalfedge* e;
	int size = 0;

	for (int i = 0; i < faces.size(); i++)
	{
		//myCentroids.push_back(new myVertex());
		e = faces[i]->adjacent_halfedge;
		size = 0;
		do
		{
			size++;
			*(faces[i]->centroid) += *(e->source->point);
			//*(myCentroids[i]->point) += *(e->source->point);
			e = e->next;
		} while (e != faces[i]->adjacent_halfedge);

		//*(myCentroids[i]->point) /= size;
		*(faces[i]->centroid) /= size;

	}


	//return myCentroids;
}




void myMesh::subdivisionCatmullClark()
{

	for (int i = 0; i < faces.size(); i++)
	{


		if (faces[i]->adjacent_halfedge == NULL)
		{
			string a;
			std::cout << "First nope : " << i << "\n";
			std::cin >> a;
		}

	}



	/**** TODO ****/
	
	computeCentroids();
	vector<myVertex*> newCentroids;

	for (int i = 0; i < faces.size(); i++)
	{
		newCentroids.push_back(new myVertex());
		newCentroids[i]->point = faces[i]->centroid;
	}
	//maintenant on a les facepoints.


	if (halfedges.size() % 2 != 0)
	{
		cout << "error number of halfedges incorrect" << endl;
	}

	int nbHalfedges = this->halfedges.size() / 2;
	std::map<std::pair<myHalfedge*, myHalfedge*>, myPoint3D*> halfEdgePoint = std::map<std::pair<myHalfedge*, myHalfedge*>, myPoint3D*>();
	//std::map<std::pair<int, int>, myPoint3D*> halfEdgePoint = std::map<std::pair<int, int>, myPoint3D*>();
	std::vector<myPoint3D*> tempVertexPoint = std::vector<myPoint3D*>();

	for (int i = 0; i < halfedges.size(); i++)
	{
		myHalfedge *HEa = halfedges[i];
		myHalfedge *HEb = (HEa->twin);

		std::pair<myHalfedge*, myHalfedge*> temp = std::pair<myHalfedge*, myHalfedge*>(HEa, HEb);
		std::pair<myHalfedge*, myHalfedge*> temp2 = std::pair<myHalfedge*, myHalfedge*>(HEb, HEa);

		std::map<std::pair<myHalfedge*, myHalfedge*>, myPoint3D*>::iterator it;
		std::map<std::pair<myHalfedge*, myHalfedge*>, myPoint3D*>::iterator it2;
		it = halfEdgePoint.find(temp);
		it2 = halfEdgePoint.find(temp2);


		int a = halfedges[i]->index;
		int b = halfedges[i]->twin->index;



		//std::map<std::pair<int, int>, myPoint3D*>::iterator it;
		//std::pair<int, int> temp = std::pair<int, int>(std::min(a,b), std::max(a,b));
		//it = halfEdgePoint.find(temp);



		// pour avoir les edgepoints.
		if (it == halfEdgePoint.end() && it2 == halfEdgePoint.end())
		{

			halfEdgePoint[temp] = new myPoint3D(0,0,0);


			
				*halfEdgePoint[temp] = (*(HEa->adjacent_face->centroid) + *(HEb->adjacent_face->centroid) + *(HEa->source->point) + *(HEb->source->point)) / 4;

			//*halfEdgePoint[temp] = (*(a->adjacent_face->centroid) + *(b->adjacent_face->centroid) + *(a->source->point) + *(b->source->point)) / 4;
			//*halfEdgePoint[temp] = (*(HEa->adjacent_face->centroid) + *(HEb->adjacent_face->centroid) + *(HEa->source->point) + *(HEb->source->point)) / 4;


		}

	}

	for (int i = 0; i < faces.size(); i++)
	{


		if (faces[i]->adjacent_halfedge == NULL)
		{
			string a;
			std::cout << "Second nope : " << i << "\n";
			std::cin >> a;
		}

	}


		//pour avoir les vertexpoints
		

	for(int i = 0; i < vertices.size();i++)
	{
		myPoint3D* q = new myPoint3D(0,0,0);
		myPoint3D* r = new myPoint3D(0, 0, 0);
		std::cout << "gor here\n";
		myHalfedge* e = vertices[i]->originof;
		std::cout << "and here\n";
		int counter = 0;
		do{
			*q += *( e->adjacent_face->centroid);
			*r += (*(e->next->source->point) + *(vertices[i]->point))/2;

			e = e->twin->next;
			counter++;

		} while (e != vertices[i]->originof);

		tempVertexPoint.push_back(new myPoint3D(0,0,0));
		*q /= counter;
		*r = ((*r * 2) / counter);

		*(tempVertexPoint[i]) = *q + *r +  *(vertices[i]->point) * (double)((counter - 3) / counter);

	}




	//set all the new vertices positions
	for (int i = 0; i < vertices.size(); i++)
	{
		*(vertices[i]->point) = *(tempVertexPoint[i]);
	}
	

	// set all the split edge points
	
	for (auto edgesToSplit : halfEdgePoint)
	{

		splitEdge(edgesToSplit.first.first, edgesToSplit.second);
	}

	// to validate the new vertex points



	int facesTemp = faces.size();
	for (int i = 0; i < facesTemp; i++)
	{
			if (faces[i]->adjacent_halfedge == NULL)
		{
			string a;
			std::cout << "Final nope : " << i << "\n";
			std::cin >> a;
		}
		

		else
			this->splitFaceQUADS(faces[i], faces[i]->centroid);
	}



	tempVertexPoint.clear();
	//delete tempVertexPoint;







}


void myMesh::triangulate()
{
	/**** TODO ****/
	//unsigned int size = this->faces.size();
	int size = this->faces.size();
	cout << "Number of faces : " << size  << "\n";
	for (int i =0; i < size;i++)
	{
		//cout << faces[i] << "\n";
		triangulate(faces[i]);
	}

	

}

void myMesh::checkFacesEdgeSize()
{

	for (int i = 0; i < faces.size(); i++)
	{
		cout << "Face " << i  << " : "<< faces[i]->getSize() << endl;
	}

}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{

	if (f == NULL)
	{
		cout << "face NULL\n";
		return false;
	}
		
	if (f->isTriangle())
	{
		cout << "already a triangle\n";
		return false;
	}
		

	//else if (f->isQuad())
	//{
		//myHalfedge* in = new myHalfedge();
		//myHalfedge* out = new myHalfedge();

		//myFace* a = new myFace();
		//myFace* b = new myFace();


		//in->next = f->adjacent_halfedge->next->next;
		//in->prev = f->adjacent_halfedge->next->next->next;
		//in->source = f->adjacent_halfedge->source;
		//in->twin = out;

		//out->next = f->adjacent_halfedge;
		//out->prev = f->adjacent_halfedge->next;
		//out->source = f->adjacent_halfedge->next->next->source;
		//out->twin = in;

		//a->adjacent_halfedge = f->adjacent_halfedge;
		//b->adjacent_halfedge = out;

		//f->adjacent_halfedge->next->next = out;
		//f->adjacent_halfedge = out->prev;
		//f->adjacent_halfedge->next->adjacent_face = a;



	//}

	/**** TODO ****/
	vector<myHalfedge*> in = vector<myHalfedge*>();
	vector<myHalfedge*> out = vector<myHalfedge*>();
	vector<myFace*> nf = vector<myFace*>();

	myHalfedge* e = new myHalfedge();
	e = f->adjacent_halfedge->next;
	int n = f->getSize();
	cout << " size of the face : " << n << endl;

	for (int i = 0; i < n - 2; i++)
	{
		in.push_back(new myHalfedge());
		out.push_back(new myHalfedge());
		nf.push_back(new myFace());
	}
	//because we needone additionnal face
	
	in[0] = f->adjacent_halfedge;

	out[n - 3] = in[0]->prev;


	myHalfedge* nextE = new myHalfedge();
	for (int i = 0; i < n - 2; i++)
	{
		nextE = e->next;

		nf[i]->adjacent_halfedge = e;

		in[i]->next = e;
		in[i]->prev = out[i];
		in[i]->source = f->adjacent_halfedge->source;
		in[i]->adjacent_face = nf[i];

		out[i]->next = in[i];
		out[i]->prev = e;
		out[i]->source = e->twin->source;
		out[i]->adjacent_face = nf[i];

		e->next = out[i];
		e->prev = in[i];
		//e->source;
		e->adjacent_face = nf[i];

		if (i != 0 && i != (n - 3))
		{
			in[i]->twin = out[i - 1];
			out[i]->twin = in[i + 1];
		}

		else if (i == 0)
		{
			//in[i]->twin = out[i - 1];
			out[i]->twin = in[i + 1];
		}

		else if (i == (n - 3))
		{
			in[i]->twin = out[i - 1];
			//out[i]->twin = in[i + 1];
		}

		else
			cout << "error error \n";

		e = nextE;

	}

	
	
	halfedges.push_back(out[0]);
	for (int i = 1; i < n - 3; i++)
	{

			halfedges.push_back(in[i]);
			halfedges.push_back(out[i]);

	}
	halfedges.push_back(in[n-3]);


	*f = *nf[0];
	for (int i = 1; i < n - 2; i++)
	{
		cout << i << " face added \n";
		faces.push_back(nf[i]);

	}

		

	return true;
}

void myMesh::inflateMesh(double dist)
{
	for(myVertex* v : vertices)
	{
		*(v->point) = *(v->point) + (*(v->normal))*dist;

	}


}