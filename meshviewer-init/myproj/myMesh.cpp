#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"




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
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{

	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}


void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
}


void myMesh::triangulate()
{
	/**** TODO ****/
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{
	/**** TODO ****/
	return false;
}

