#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
//#include <OpenGL/gl3.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[2];
int gWidth, gHeight;

GLint modelingMatrixLoc[2];
GLint viewingMatrixLoc[2];
GLint projectionMatrixLoc[2];
GLint eyePosLoc[2];

GLint kdLoc[2];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::vec3 eyePos(0, 0, 0);
glm::mat4 cubeModelingMatrix;
glm::mat4 octaModelingMatrix;
glm::mat4 tetraModelingMatrix;

int activeProgramIndex = 0;

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        return h1 ^ h2;
    }
};

struct vec3_hash
{
    size_t operator()(const glm::vec3& v) const {
        size_t h1 = std::hash<float>{}(v.x);
        size_t h2 = std::hash<float>{}(v.y);
        size_t h3 = std::hash<float>{}(v.z);

        // Combine the hashes using a simple hash function
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct Vertex
{
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Texture
{
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
    GLfloat u, v;
};

struct Normal
{
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int n[], int size) {
        if (size == 3) {
            vIndex.resize(3);
            nIndex.resize(3);
            vIndex[0] = v[0];
            vIndex[1] = v[1];
            vIndex[2] = v[2];
            nIndex[0] = n[0];
            nIndex[1] = n[1];
            nIndex[2] = n[2];
        }
        else if (size == 4) {
			vIndex.resize(4);
			nIndex.resize(4);
			vIndex[0] = v[0];
			vIndex[1] = v[1];
			vIndex[2] = v[2];
			vIndex[3] = v[3];
			nIndex[0] = n[0];
			nIndex[1] = n[1];
			nIndex[2] = n[2];
			nIndex[3] = n[3];
		}
	}
    Face(int v[], int size) {
        if (size == 3) {
            vIndex.resize(3);
            vIndex[0] = v[0];
            vIndex[1] = v[1];
            vIndex[2] = v[2];
        }
        else if (size == 4) {
            vIndex.resize(4);
            vIndex[0] = v[0];
            vIndex[1] = v[1];
            vIndex[2] = v[2];
            vIndex[3] = v[3];
        }
    }
    vector<GLuint> vIndex, nIndex;
};

bool pause = false;

enum DrawingMode { SOLID, LINE, WIREFRAME };
DrawingMode mode = SOLID;

vector<vector<Vertex>> cubeVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> cubeNormals = vector<vector<Normal>>(5);
vector<vector<Face>> cubeQuadFaces = vector<vector<Face>>(5);
vector<vector<Face>> cubeFaces = vector<vector<Face>>(5);

vector<vector<Vertex>> cubeLineFlatVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> cubeLineFlatNormals = vector<vector<Normal>>(5);
vector<vector<Face>> cubeLineFlatFaces = vector<vector<Face>>(5);
vector<GLuint> cubeLineVertexAttribBuffer = vector<GLuint>(5);
vector<GLuint> cubeLineIndexBuffer = vector<GLuint>(5);
vector<int> cubeLineVertexDataSizeInBytes = vector<int>(5);
vector<int> cubeLineNormalDataSizeInBytes = vector<int>(5);
vector<GLuint> cubeLineVAOs = vector<GLuint>(5);

vector<vector<Vertex>> cubeFlatVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> cubeFlatNormals = vector<vector<Normal>>(5);
vector<vector<Face>> cubeFlatFaces = vector<vector<Face>>(5);
vector<GLuint> cubeVertexAttribBuffer = vector<GLuint>(5);
vector<GLuint> cubeIndexBuffer = vector<GLuint>(5);
vector<int> cubeVertexDataSizeInBytes = vector<int>(5);
vector<int> cubeNormalDataSizeInBytes = vector<int>(5);
vector<GLuint> cubeVAOs = vector<GLuint>(5);

glm::vec3 cubeColor(1, 0.2, 0.2);
int cubeLvl = 1;
float cubeAngle = 0.0f;


vector<vector<Vertex>> octaVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> octaNormals = vector<vector<Normal>>(5);
vector<vector<Face>> octaQuadFaces = vector<vector<Face>>(5);
vector<vector<Face>> octaFaces = vector<vector<Face>>(5);

vector<vector<Vertex>> octaLineFlatVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> octaLineFlatNormals = vector<vector<Normal>>(5);
vector<vector<Face>> octaLineFlatFaces = vector<vector<Face>>(5);
vector<GLuint> octaLineVertexAttribBuffer = vector<GLuint>(5);
vector<GLuint> octaLineIndexBuffer = vector<GLuint>(5);
vector<int> octaLineVertexDataSizeInBytes = vector<int>(5);
vector<int> octaLineNormalDataSizeInBytes = vector<int>(5);
vector<GLuint> octaLineVAOs = vector<GLuint>(5);

vector<vector<Vertex>> octaFlatVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> octaFlatNormals = vector<vector<Normal>>(5);
vector<vector<Face>> octaFlatFaces = vector<vector<Face>>(5);
vector<GLuint> octaVertexAttribBuffer = vector<GLuint>(5);
vector<GLuint> octaIndexBuffer = vector<GLuint>(5);
vector<int> octaVertexDataSizeInBytes = vector<int>(5);
vector<int> octaNormalDataSizeInBytes = vector<int>(5);
vector<GLuint> octaVAOs = vector<GLuint>(5);

glm::vec3 octaColor(0.2, 1, 0.2);
int octaLvl = 0;
float octaAngle = 0.0f;

vector<vector<Vertex>> tetraVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> tetraNormals = vector<vector<Normal>>(5);
vector<vector<Face>> tetraQuadFaces = vector<vector<Face>>(5);
vector<vector<Face>> tetraFaces = vector<vector<Face>>(5);

vector<vector<Vertex>> tetraLineFlatVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> tetraLineFlatNormals = vector<vector<Normal>>(5);
vector<vector<Face>> tetraLineFlatFaces = vector<vector<Face>>(5);
vector<GLuint> tetraLineVertexAttribBuffer = vector<GLuint>(5);
vector<GLuint> tetraLineIndexBuffer = vector<GLuint>(5);
vector<int> tetraLineVertexDataSizeInBytes = vector<int>(5);
vector<int> tetraLineNormalDataSizeInBytes = vector<int>(5);
vector<GLuint> tetraLineVAOs = vector<GLuint>(5);

vector<vector<Vertex>> tetraFlatVertices = vector<vector<Vertex>>(5);
vector<vector<Normal>> tetraFlatNormals = vector<vector<Normal>>(5);
vector<vector<Face>> tetraFlatFaces = vector<vector<Face>>(5);
vector<GLuint> tetraVertexAttribBuffer = vector<GLuint>(5);
vector<GLuint> tetraIndexBuffer = vector<GLuint>(5);
vector<int> tetraVertexDataSizeInBytes = vector<int>(5);
vector<int> tetraNormalDataSizeInBytes = vector<int>(5);
vector<GLuint> tetraVAOs = vector<GLuint>(5);

glm::vec3 tetraColor(0.2, 0.2, 1);
int tetraLvl = 0;
float tetraAngle = 0.0f;

bool writeOBJFile(const std::string& filename, const std::vector<Vertex>& vertices, const std::vector<Face>& faces) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to create file " << filename << std::endl;
        return false;
    }

    // Write vertices
    for (const auto& vertex : vertices) {
        file << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }

    // Write faces
    for (const auto& face : faces) {
        file << "f ";
        for (int i = 0; i < 4; ++i)
            file << face.vIndex[i] + 1 << " "; // OBJ format starts indexing from 1
        file << "\n";
    }

    file.close();
    return true;
}

//GLuint gVertexAttribBuffer, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
//int gVertexDataSizeInBytes, gNormalDataSizeInBytes;

bool ParseQuadObj(const string& fileName, vector<Vertex>& gVertices, vector<Normal>& gNormals, vector<Face>& gFaces)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == 'v')
                {
                    if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
                    char c;
                    int vIndex[4], nIndex[4];
                    str >> vIndex[0]; str >> c >> c; // consume "//"
                    str >> nIndex[0];
                    str >> vIndex[1]; str >> c >> c; // consume "//"
                    str >> nIndex[1];
                    str >> vIndex[2]; str >> c >> c; // consume "//"
                    str >> nIndex[2];
                    str >> vIndex[3]; str >> c >> c; // consume "//"
                    str >> nIndex[3];

                    assert(vIndex[0] == nIndex[0] &&
                        vIndex[1] == nIndex[1] &&
                        vIndex[2] == nIndex[2]); // a limitation for now

                    // make indices start from 0
                    for (int c = 0; c < 4; ++c)
                    {
                        vIndex[c] -= 1;
                        nIndex[c] -= 1;
                    }

                    gFaces.push_back(Face(vIndex, nIndex, 4));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    /*
    for (int i = 0; i < gVertices.size(); ++i)
    {
        Vector3 n;

        for (int j = 0; j < gFaces.size(); ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (gFaces[j].vIndex[k] == i)
                {
                    // face j contains vertex i
                    Vector3 a(gVertices[gFaces[j].vIndex[0]].x,
                              gVertices[gFaces[j].vIndex[0]].y,
                              gVertices[gFaces[j].vIndex[0]].z);

                    Vector3 b(gVertices[gFaces[j].vIndex[1]].x,
                              gVertices[gFaces[j].vIndex[1]].y,
                              gVertices[gFaces[j].vIndex[1]].z);

                    Vector3 c(gVertices[gFaces[j].vIndex[2]].x,
                              gVertices[gFaces[j].vIndex[2]].y,
                              gVertices[gFaces[j].vIndex[2]].z);

                    Vector3 ab = b - a;
                    Vector3 ac = c - a;
                    Vector3 normalFromThisFace = (ab.cross(ac)).getNormalized();
                    n += normalFromThisFace;
                }

            }
        }

        n.normalize();

        gNormals.push_back(Normal(n.x, n.y, n.z));
    }
    */

    assert(gVertices.size() == gNormals.size());

    return true;
}

bool ParseObj(const string& fileName, vector<Vertex>& gVertices, vector<Normal>& gNormals, vector<Face>& gFaces)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == 'v')
                {
                    if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
					char c;
					int vIndex[3],  nIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0]; 
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1]; 
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2]; 

					assert(vIndex[0] == nIndex[0] &&
						   vIndex[1] == nIndex[1] &&
						   vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
					}

                    gFaces.push_back(Face(vIndex, nIndex, 3));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

	/*
	for (int i = 0; i < gVertices.size(); ++i)
	{
		Vector3 n;

		for (int j = 0; j < gFaces.size(); ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				if (gFaces[j].vIndex[k] == i)
				{
					// face j contains vertex i
					Vector3 a(gVertices[gFaces[j].vIndex[0]].x, 
							  gVertices[gFaces[j].vIndex[0]].y,
							  gVertices[gFaces[j].vIndex[0]].z);

					Vector3 b(gVertices[gFaces[j].vIndex[1]].x, 
							  gVertices[gFaces[j].vIndex[1]].y,
							  gVertices[gFaces[j].vIndex[1]].z);

					Vector3 c(gVertices[gFaces[j].vIndex[2]].x, 
							  gVertices[gFaces[j].vIndex[2]].y,
							  gVertices[gFaces[j].vIndex[2]].z);

					Vector3 ab = b - a;
					Vector3 ac = c - a;
					Vector3 normalFromThisFace = (ab.cross(ac)).getNormalized();
					n += normalFromThisFace;
				}

			}
		}

		n.normalize();

		gNormals.push_back(Normal(n.x, n.y, n.z));
	}
	*/

	assert(gVertices.size() == gNormals.size());

    return true;
}


void catmullClark(vector<Vertex>& oldVertices, vector<Face>& oldFaces, vector<Vertex>& newVertices, vector<Face>& newFaces)
{
    newVertices.clear();
    newFaces.clear();

    struct PointTracker;
    struct FaceTracker;
    struct EdgeTracker;
    
    struct PointTracker {
        glm::vec3 position;
        vector<FaceTracker*> faces;
        vector<EdgeTracker*> edges;
        glm::vec3 newPosition;
    };

    struct FaceTracker {
        vector<EdgeTracker*> edges;
        vector<PointTracker*> points;
        glm::vec3 facePoint;
    };

    struct EdgeTracker {
		vector<FaceTracker*> faces;
		vector<PointTracker*> points;
		glm::vec3 edgePoint;
        glm::vec3 midPoint;
	};
    

    unordered_map<int, PointTracker*> pointMap;
    unordered_map<int, FaceTracker*> faceMap;
    unordered_map<pair<int, int>, EdgeTracker*, pair_hash> edgeMap;
    
    for (int faceIndex = 0; faceIndex < oldFaces.size(); faceIndex++)
    {
        Face& curFace = oldFaces[faceIndex];
        vector<PointTracker*> faceCorners;
        
        faceMap[faceIndex] = new FaceTracker();

        for (int vertexIndex = 0; vertexIndex < curFace.vIndex.size(); vertexIndex++)
        {
            int positionIndex = curFace.vIndex[vertexIndex];

            PointTracker* pointObject = new PointTracker();
            
            if (pointMap.find(positionIndex) == pointMap.end())
            {
				pointObject->position = glm::vec3(oldVertices[positionIndex].x, oldVertices[positionIndex].y, oldVertices[positionIndex].z);
				pointMap[positionIndex] = pointObject;
                //pointObject = &pointMap[positionIndex];
			}
            else
            {
				pointObject = pointMap[positionIndex];
			}

            pointObject->faces.push_back(faceMap[faceIndex]);

            faceCorners.push_back(pointObject);
		}

        faceMap[faceIndex]->points = faceCorners;

        glm::vec3 newFacePoint(0, 0, 0);

        for (int i = 0; i < faceMap[faceIndex]->points.size(); i++)
        {
			newFacePoint += faceMap[faceIndex]->points[i]->position;
		}
        newFacePoint /= faceMap[faceIndex]->points.size();
        faceMap[faceIndex]->facePoint = newFacePoint;

        vector<EdgeTracker*> IncidentEdges;

        for (int edgeIndex = 0; edgeIndex < curFace.vIndex.size(); edgeIndex++)
        {
            pair<int, int> edge;

            if (curFace.vIndex.size() == 3)
            {
                if (edgeIndex == 0)
                {
					edge = make_pair(curFace.vIndex[0], curFace.vIndex[1]);
				}
                else if (edgeIndex == 1)
                {
					edge = make_pair(curFace.vIndex[1], curFace.vIndex[2]);
				}
                else if (edgeIndex == 2)
                {
					edge = make_pair(curFace.vIndex[2], curFace.vIndex[0]);
				}
            }
            else if (curFace.vIndex.size() == 4)
            {
                if (edgeIndex == 0)
				{
                    edge = make_pair(curFace.vIndex[0], curFace.vIndex[1]);
                }
                else if (edgeIndex == 1)
                {
					edge = make_pair(curFace.vIndex[1], curFace.vIndex[2]);
				}
                else if (edgeIndex == 2)
                {
					edge = make_pair(curFace.vIndex[2], curFace.vIndex[3]);
				}
                else if (edgeIndex == 3)
                {
					edge = make_pair(curFace.vIndex[3], curFace.vIndex[0]);
				}
            }

            if (edge.first > edge.second)
                swap(edge.first, edge.second);

            EdgeTracker *newEdge = new EdgeTracker();

            if (edgeMap.find(edge) == edgeMap.end())
            {
                newEdge->points.push_back(pointMap[edge.first]);
                newEdge->points.push_back(pointMap[edge.second]);

                newEdge->points[0]->edges.push_back(newEdge);
                newEdge->points[1]->edges.push_back(newEdge);

                edgeMap[edge] = newEdge;
                //newEdge = &edgeMap[edge];
            }
            else
            {
                newEdge = edgeMap[edge];
            }

            newEdge->faces.push_back(faceMap[faceIndex]);

            IncidentEdges.push_back(newEdge);
        }

        faceMap[faceIndex]->edges = IncidentEdges;
    }


    for (auto& kvPair : edgeMap)
    {
        
        EdgeTracker* edge = kvPair.second;

        glm::vec3 newEdgePoint(0, 0, 0);
        int count = 0;

        for (int i = 0; i < edge->faces.size(); i++)
        {
            glm::vec3 facePoint = edge->faces[i]->facePoint;
            newEdgePoint += facePoint;
            count++;
        }

        for (int i = 0; i < edge->points.size(); i++)
        {
            glm::vec3 endPoint = edge->points[i]->position;
            newEdgePoint += endPoint;
            count++;
		}

        newEdgePoint /= count;
        edge->edgePoint = newEdgePoint;

        count = 0;
        glm::vec3 newMidPoint (0, 0, 0);

        for (int i = 0; i < edge->points.size(); i++)
        {
			glm::vec3 endPoint = edge->points[i]->position;
			newMidPoint += endPoint;
			count++;
		}
        newMidPoint /= count;

        edge->midPoint = newMidPoint;
	}
    
    for (int i = 0; i < oldVertices.size(); i++)
    {
        PointTracker* point = pointMap[i];
        int n = point->faces.size();
        
        glm::vec3 F(0, 0, 0);
        glm::vec3 R(0, 0, 0);
        glm::vec3 P = point->position;

        for (int i = 0; i < point->faces.size(); i++)
        {
			F += point->faces[i]->facePoint;
		}
        F /= n;

        for (int i = 0; i < point->edges.size(); i++)
        {
            R += point->edges[i]->midPoint;
        }
        R /= n;

        R *= 2;
        P *= (n - 3);

        glm::vec3 newPoint = (F + R + P);
        newPoint /= n;

        point->newPosition = newPoint;


    }

    unordered_map<glm::vec3, int, vec3_hash> vertexIndices;
    int index = 0;

    for (int faceIndex = 0; faceIndex < faceMap.size(); faceIndex++)
    {

        FaceTracker* face = faceMap[faceIndex];

        for (int pointIndex = 0; pointIndex < face->points.size(); pointIndex++)
        {
			PointTracker* point = face->points[pointIndex];

            glm::vec3 a = point->newPosition;
            glm::vec3 b = face->edges[pointIndex % face->edges.size()]->edgePoint;
            glm::vec3 c = face->facePoint;
            glm::vec3 d = face->edges[(pointIndex - 1 + face->edges.size()) % face->edges.size()]->edgePoint;

            int ia, ib, ic, id;

            if (vertexIndices.find(a) == vertexIndices.end())
            {
				ia = index;
				vertexIndices[a] = index;
				index++;
                newVertices.push_back(Vertex(a.x, a.y, a.z));
			}
            else
            {
				ia = vertexIndices[a];
			}

            if (vertexIndices.find(b) == vertexIndices.end())
            {
				ib = index;
				vertexIndices[b] = index;
				index++;
				newVertices.push_back(Vertex(b.x, b.y, b.z));
			}
            else
            {
				ib = vertexIndices[b];
			}

            if (vertexIndices.find(c) == vertexIndices.end())
            {
				ic = index;
				vertexIndices[c] = index;
				index++;
				newVertices.push_back(Vertex(c.x, c.y, c.z));
			}
            else
            {
				ic = vertexIndices[c];
			}

            if (vertexIndices.find(d) == vertexIndices.end())
            {
				id = index;
				vertexIndices[d] = index;
				index++;
				newVertices.push_back(Vertex(d.x, d.y, d.z));
			}
            else
            {
				id = vertexIndices[d];
			}

            int vIndex[4] = { ia, ib, ic, id };
            int nIndex[4] = { ia, ib, ic, id };
			newFaces.push_back(Face(vIndex, nIndex, 4));
		}
    }
}


void QuadFacetoTriFace(const vector<Face>& quadFaces, vector<Face>& triFaces)
{
    triFaces.clear();
    for (int i = 0; i < quadFaces.size(); ++i)
    {
        int vIndex[3], nIndex[3];
        vIndex[0] = quadFaces[i].vIndex[0];
        vIndex[1] = quadFaces[i].vIndex[1];
        vIndex[2] = quadFaces[i].vIndex[2];
        nIndex[0] = quadFaces[i].nIndex[0];
        nIndex[1] = quadFaces[i].nIndex[1];
        nIndex[2] = quadFaces[i].nIndex[2];
        triFaces.push_back(Face(vIndex, nIndex, 3));

        vIndex[0] = quadFaces[i].vIndex[2];
        vIndex[1] = quadFaces[i].vIndex[3];
        vIndex[2] = quadFaces[i].vIndex[0];
        nIndex[0] = quadFaces[i].nIndex[2];
        nIndex[1] = quadFaces[i].nIndex[3];
        nIndex[2] = quadFaces[i].nIndex[0];
        triFaces.push_back(Face(vIndex, nIndex, 3));
	}
}   

void ComputeNormalsFS(vector<Vertex>& oldVertices, vector<Face>& oldFaces, vector<Vertex>& newVertices, vector<Normal>& newNormals, vector<Face>& newFaces, bool isQuad)
{
    newVertices.clear();
    newNormals.clear();
    newFaces.clear();
    int index = 0;
    for (int i = 0; i < oldFaces.size(); i++) {
        Face &face = oldFaces[i];
        Vertex v1 = oldVertices[face.vIndex[0]];
        Vertex v2 = oldVertices[face.vIndex[1]];
        Vertex v3 = oldVertices[face.vIndex[2]];
        glm::vec3 a(v1.x, v1.y, v1.z);
        glm::vec3 b(v2.x, v2.y, v2.z);
        glm::vec3 c(v3.x, v3.y, v3.z);
        glm::vec3 normal = glm::cross(b - a, c - a);
        normal = glm::normalize(normal);
        newVertices.push_back(v1);
        newVertices.push_back(v2);
        newVertices.push_back(v3);
        newNormals.push_back(Normal(normal.x, normal.y, normal.z));
        newNormals.push_back(Normal(normal.x, normal.y, normal.z));
        newNormals.push_back(Normal(normal.x, normal.y, normal.z));
        int vIndex[3] = { index, index + 1, index + 2 };
        int nIndex[3] = { index, index + 1, index + 2 };
        newFaces.push_back(Face(vIndex, nIndex, 3));
        index += 3;
    }
    if (isQuad) {
        for (int i = 0; i < newFaces.size(); i += 2) {
			Face &face1 = newFaces[i];
			Face &face2 = newFaces[i + 1];
			Normal &n1 = newNormals[face1.nIndex[0]];
            Normal &n2 = newNormals[face2.nIndex[0]];
            Normal avg((n1.x + n2.x) / 2, (n1.y + n2.y) / 2, (n1.z + n2.z) / 2);
            newNormals[face1.nIndex[0]] = avg;
            newNormals[face1.nIndex[1]] = avg;
            newNormals[face1.nIndex[2]] = avg;
            newNormals[face2.nIndex[0]] = avg;
            newNormals[face2.nIndex[1]] = avg;
            newNormals[face2.nIndex[2]] = avg;
		}
	}   
}

void ComputeQuadNormalsFS(vector<Vertex>& oldVertices, vector<Face>& oldFaces, vector<Vertex>& newVertices, vector<Normal>& newNormals, vector<Face>& newFaces)
{
    newVertices.clear();
    newNormals.clear();
    newFaces.clear();
    int index = 0;
    for (int i = 0; i < oldFaces.size(); i++) {
		Face &face = oldFaces[i];
		Vertex v1 = oldVertices[face.vIndex[0]];
		Vertex v2 = oldVertices[face.vIndex[1]];
		Vertex v3 = oldVertices[face.vIndex[2]];
		Vertex v4 = oldVertices[face.vIndex[3]];
        glm::vec3 a(v1.x, v1.y, v1.z);
        glm::vec3 b(v2.x, v2.y, v2.z);
        glm::vec3 c(v3.x, v3.y, v3.z);
        glm::vec3 d(v4.x, v4.y, v4.z);
        glm::vec3 normal1 = glm::cross(b - a, c - a);
        glm::vec3 normal2 = glm::cross(c - a, d - a);
        normal1 = glm::normalize(normal1);
        normal2 = glm::normalize(normal2);
        glm::vec3 avg((normal1.x + normal2.x) / 2, (normal1.y + normal2.y) / 2, (normal1.z + normal2.z) / 2);
        newVertices.push_back(v1);
        newVertices.push_back(v2);
        newVertices.push_back(v3);
        newVertices.push_back(v4);
        newNormals.push_back(Normal(avg.x, avg.y, avg.z));
        newNormals.push_back(Normal(avg.x, avg.y, avg.z));
        newNormals.push_back(Normal(avg.x, avg.y, avg.z));
        newNormals.push_back(Normal(avg.x, avg.y, avg.z));
        int vIndex[4] = { index, index + 1, index + 2, index + 3 };
        int nIndex[4] = { index, index + 1, index + 2, index + 3 };
        newFaces.push_back(Face(vIndex, nIndex, 4));
        index += 4;
    }
}

bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

GLuint createVS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

	return fs;
}

void initShaders()
{
	// Create the programs

    gProgram[0] = glCreateProgram();
	gProgram[1] = glCreateProgram();

	// Create the shaders for both programs

    GLuint vs1 = createVS("vert.glsl");
    GLuint fs1 = createFS("frag.glsl");

	GLuint vs2 = createVS("vert2.glsl");
	GLuint fs2 = createFS("frag2.glsl");

	// Attach the shaders to the programs

	glAttachShader(gProgram[0], vs1);
	glAttachShader(gProgram[0], fs1);

	glAttachShader(gProgram[1], vs2);
	glAttachShader(gProgram[1], fs2);

	// Link the programs

    glLinkProgram(gProgram[0]);
	GLint status;
	glGetProgramiv(gProgram[0], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}

	glLinkProgram(gProgram[1]);
	glGetProgramiv(gProgram[1], GL_LINK_STATUS, &status);

	if (status != GL_TRUE)
	{
		cout << "Program link failed" << endl;
		exit(-1);
	}

	// Get the locations of the uniform variables from both programs

	for (int i = 0; i < 2; ++i)
	{
		modelingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "modelingMatrix");
		viewingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "viewingMatrix");
		projectionMatrixLoc[i] = glGetUniformLocation(gProgram[i], "projectionMatrix");
		eyePosLoc[i] = glGetUniformLocation(gProgram[i], "eyePos");
        kdLoc[i] = glGetUniformLocation(gProgram[i], "kd");
	}
}


void initVBO(GLuint &vao, vector<Vertex>& gVertices, int& gVertexDataSizeInBytes, vector<Normal>& gNormals, int& gNormalDataSizeInBytes, vector<Face>& gFaces, GLuint& gVertexAttribBuffer, GLuint& gIndexBuffer)
{
    glGenVertexArrays(1, &vao);
    assert(vao > 0);
    glBindVertexArray(vao);
    cout << "vao = " << vao << endl;

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer);
	glGenBuffers(1, &gIndexBuffer);

	assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
	gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
	int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);
	GLfloat* vertexData = new GLfloat [gVertices.size() * 3];
	GLfloat* normalData = new GLfloat [gNormals.size() * 3];
	GLuint* indexData = new GLuint [gFaces.size() * 3];

    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

	for (int i = 0; i < gVertices.size(); ++i)
	{
		vertexData[3*i] = gVertices[i].x;
		vertexData[3*i+1] = gVertices[i].y;
		vertexData[3*i+2] = gVertices[i].z;
        /*
        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
        */
	}
    /*
    std::cout << "minX = " << minX << std::endl;
    std::cout << "maxX = " << maxX << std::endl;
    std::cout << "minY = " << minY << std::endl;
    std::cout << "maxY = " << maxY << std::endl;
    std::cout << "minZ = " << minZ << std::endl;
    std::cout << "maxZ = " << maxZ << std::endl;
    */
	for (int i = 0; i < gNormals.size(); ++i)
	{
		normalData[3*i] = gNormals[i].x;
		normalData[3*i+1] = gNormals[i].y;
		normalData[3*i+2] = gNormals[i].z;
	}

	for (int i = 0; i < gFaces.size(); ++i)
	{
		indexData[3*i] = gFaces[i].vIndex[0];
		indexData[3*i+1] = gFaces[i].vIndex[1];
		indexData[3*i+2] = gFaces[i].vIndex[2];
	}


	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

	// done copying; can free now
	delete[] vertexData;
	delete[] normalData;
	delete[] indexData;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
}

void initLineVBO(GLuint& vao, vector<Vertex>& gVertices, int& gVertexDataSizeInBytes, vector<Normal>& gNormals, int& gNormalDataSizeInBytes, vector<Face>& gFaces, GLuint& gVertexAttribBuffer, GLuint& gIndexBuffer)
{
    glGenVertexArrays(1, &vao);
    assert(vao > 0);
    glBindVertexArray(vao);
    cout << "vao = " << vao << endl;

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    assert(glGetError() == GL_NONE);

    glGenBuffers(1, &gVertexAttribBuffer);
    glGenBuffers(1, &gIndexBuffer);

    assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

    glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

    gVertexDataSizeInBytes = gVertices.size() * 3 * sizeof(GLfloat);
    gNormalDataSizeInBytes = gNormals.size() * 3 * sizeof(GLfloat);
    int indexDataSizeInBytes;
    if (gFaces[0].vIndex.size() == 3)
        indexDataSizeInBytes = gFaces.size() * 6 * sizeof(GLuint);
    else// if (gFaces[0].vIndex.size() == 4)
		indexDataSizeInBytes = gFaces.size() * 8 * sizeof(GLuint);
    GLfloat* vertexData = new GLfloat[gVertices.size() * 3];
    GLfloat* normalData = new GLfloat[gNormals.size() * 3];
    GLuint* indexData;
    if (gFaces[0].vIndex.size() == 3)
        indexData = new GLuint[gFaces.size() * 6];
    else// if (gFaces[0].vIndex.size() == 4)
		indexData = new GLuint[gFaces.size() * 8];

    float minX = 1e6, maxX = -1e6;
    float minY = 1e6, maxY = -1e6;
    float minZ = 1e6, maxZ = -1e6;

    for (int i = 0; i < gVertices.size(); ++i)
    {
        vertexData[3 * i] = gVertices[i].x;
        vertexData[3 * i + 1] = gVertices[i].y;
        vertexData[3 * i + 2] = gVertices[i].z;
        /*
        minX = std::min(minX, gVertices[i].x);
        maxX = std::max(maxX, gVertices[i].x);
        minY = std::min(minY, gVertices[i].y);
        maxY = std::max(maxY, gVertices[i].y);
        minZ = std::min(minZ, gVertices[i].z);
        maxZ = std::max(maxZ, gVertices[i].z);
        */
    }
    /*
    std::cout << "minX = " << minX << std::endl;
    std::cout << "maxX = " << maxX << std::endl;
    std::cout << "minY = " << minY << std::endl;
    std::cout << "maxY = " << maxY << std::endl;
    std::cout << "minZ = " << minZ << std::endl;
    std::cout << "maxZ = " << maxZ << std::endl;
    */
    for (int i = 0; i < gNormals.size(); ++i)
    {
        normalData[3 * i] = gNormals[i].x;
        normalData[3 * i + 1] = gNormals[i].y;
        normalData[3 * i + 2] = gNormals[i].z;
    }
    if (gFaces[0].vIndex.size() == 3)
    {
        for (int i = 0; i < gFaces.size(); ++i)
        {
            indexData[6 * i] = gFaces[i].vIndex[0];
            indexData[6 * i + 1] = gFaces[i].vIndex[1];
            indexData[6 * i + 2] = gFaces[i].vIndex[1];
            indexData[6 * i + 3] = gFaces[i].vIndex[2];
            indexData[6 * i + 4] = gFaces[i].vIndex[2];
            indexData[6 * i + 5] = gFaces[i].vIndex[0];
        }
    }
    else if (gFaces[0].vIndex.size() == 4)
    {
        for (int i = 0; i < gFaces.size(); i++)
        {
            indexData[8 * i] = gFaces[i].vIndex[0];
			indexData[8 * i + 1] = gFaces[i].vIndex[1];
			indexData[8 * i + 2] = gFaces[i].vIndex[1];
			indexData[8 * i + 3] = gFaces[i].vIndex[2];
			indexData[8 * i + 4] = gFaces[i].vIndex[2];
			indexData[8 * i + 5] = gFaces[i].vIndex[3];
			indexData[8 * i + 6] = gFaces[i].vIndex[3];
			indexData[8 * i + 7] = gFaces[i].vIndex[0];
        }
    }
    


    glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
    glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    // done copying; can free now
    delete[] vertexData;
    delete[] normalData;
    delete[] indexData;

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
}

void init() 
{
	ParseQuadObj("cube.obj", cubeVertices[0], cubeNormals[0], cubeQuadFaces[0]);

    ParseObj("octahedron.obj", octaVertices[0], octaNormals[0], octaFaces[0]);
    octaQuadFaces[0] = octaFaces[0];

    ParseObj("tetrahedron.obj", tetraVertices[0], tetraNormals[0], tetraFaces[0]);
    tetraQuadFaces[0] = tetraFaces[0];

    glEnable(GL_DEPTH_TEST);
    initShaders();

    for (int i = 0; i < 5; i++)
    {
        ComputeQuadNormalsFS(cubeVertices[i], cubeQuadFaces[i], cubeLineFlatVertices[i], cubeLineFlatNormals[i], cubeLineFlatFaces[i]);
        initLineVBO(cubeLineVAOs[i], cubeLineFlatVertices[i], cubeLineVertexDataSizeInBytes[i], cubeLineFlatNormals[i], cubeLineVertexDataSizeInBytes[i], cubeLineFlatFaces[i], cubeLineVertexAttribBuffer[i], cubeLineIndexBuffer[i]);
        QuadFacetoTriFace(cubeQuadFaces[i], cubeFaces[i]);
        ComputeNormalsFS(cubeVertices[i], cubeFaces[i], cubeFlatVertices[i], cubeFlatNormals[i], cubeFlatFaces[i], true);
        initVBO(cubeVAOs[i], cubeFlatVertices[i], cubeVertexDataSizeInBytes[i], cubeFlatNormals[i], cubeVertexDataSizeInBytes[i], cubeFlatFaces[i], cubeVertexAttribBuffer[i], cubeIndexBuffer[i]);
        if (i != 4)
            catmullClark(cubeVertices[i], cubeQuadFaces[i], cubeVertices[i + 1], cubeQuadFaces[i + 1]);
    }

    for (int i = 0; i < 5; i++)
    {
        if (i > 0)
        {
            ComputeQuadNormalsFS(octaVertices[i], octaQuadFaces[i], octaLineFlatVertices[i], octaLineFlatNormals[i], octaLineFlatFaces[i]);
            initLineVBO(octaLineVAOs[i], octaLineFlatVertices[i], octaLineVertexDataSizeInBytes[i], octaLineFlatNormals[i], octaLineVertexDataSizeInBytes[i], octaLineFlatFaces[i], octaLineVertexAttribBuffer[i], octaLineIndexBuffer[i]);
            QuadFacetoTriFace(octaQuadFaces[i], octaFaces[i]);
        }
		ComputeNormalsFS(octaVertices[i], octaFaces[i], octaFlatVertices[i], octaFlatNormals[i], octaFlatFaces[i], (i > 0 ? true : false));
        if (i == 0)
        {
            octaLineFlatVertices[0] = octaFlatVertices[0];
            octaLineFlatNormals[0] = octaFlatNormals[0];
            octaLineFlatFaces[0] = octaFlatFaces[0];
            initLineVBO(octaLineVAOs[0], octaLineFlatVertices[0], octaLineVertexDataSizeInBytes[0], octaLineFlatNormals[0], octaLineVertexDataSizeInBytes[0], octaLineFlatFaces[0], octaLineVertexAttribBuffer[0], octaLineIndexBuffer[0]);
        }
		initVBO(octaVAOs[i], octaFlatVertices[i], octaVertexDataSizeInBytes[i], octaFlatNormals[i], octaVertexDataSizeInBytes[i], octaFlatFaces[i], octaVertexAttribBuffer[i], octaIndexBuffer[i]);
		if (i != 4)
			catmullClark(octaVertices[i], octaQuadFaces[i], octaVertices[i + 1], octaQuadFaces[i + 1]);
	}
    
    for (int i = 0; i < 5; i++)
    {
        if (i > 0)
        {
            ComputeQuadNormalsFS(tetraVertices[i], tetraQuadFaces[i], tetraLineFlatVertices[i], tetraLineFlatNormals[i], tetraLineFlatFaces[i]);
			initLineVBO(tetraLineVAOs[i], tetraLineFlatVertices[i], tetraLineVertexDataSizeInBytes[i], tetraLineFlatNormals[i], tetraLineVertexDataSizeInBytes[i], tetraLineFlatFaces[i], tetraLineVertexAttribBuffer[i], tetraLineIndexBuffer[i]);
			QuadFacetoTriFace(tetraQuadFaces[i], tetraFaces[i]);
        }
        ComputeNormalsFS(tetraVertices[i], tetraFaces[i], tetraFlatVertices[i], tetraFlatNormals[i], tetraFlatFaces[i], (i > 0 ? true : false));
        if (i == 0)
        {
			tetraLineFlatVertices[0] = tetraFlatVertices[0];
			tetraLineFlatNormals[0] = tetraFlatNormals[0];
			tetraLineFlatFaces[0] = tetraFlatFaces[0];
			initLineVBO(tetraLineVAOs[0], tetraLineFlatVertices[0], tetraLineVertexDataSizeInBytes[0], tetraLineFlatNormals[0], tetraLineVertexDataSizeInBytes[0], tetraLineFlatFaces[0], tetraLineVertexAttribBuffer[0], tetraLineIndexBuffer[0]);
		}
        initVBO(tetraVAOs[i], tetraFlatVertices[i], tetraVertexDataSizeInBytes[i], tetraFlatNormals[i], tetraVertexDataSizeInBytes[i], tetraFlatFaces[i], tetraVertexAttribBuffer[i], tetraIndexBuffer[i]);
		if (i != 4)
			catmullClark(tetraVertices[i], tetraQuadFaces[i], tetraVertices[i + 1], tetraQuadFaces[i + 1]);
    }
}

void drawModel(GLuint& vao, GLuint& gVertexAttribBuffer, GLuint& gIndexBuffer, int& gVertexDataSizeInBytes, vector<Face>& gFaces)
{
    glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
}

void drawLineModel(GLuint& vao, GLuint& gVertexAttribBuffer, GLuint& gIndexBuffer, int& gVertexDataSizeInBytes, vector<Face>& gFaces)
{
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

    if (gFaces[0].vIndex.size() == 3)
		glDrawElements(GL_LINES, gFaces.size() * 6, GL_UNSIGNED_INT, 0);
	else if (gFaces[0].vIndex.size() == 4)
		glDrawElements(GL_LINES, gFaces.size() * 8, GL_UNSIGNED_INT, 0);
}

void display()
{
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    //glDisable(GL_CULL_FACE);

    glUseProgram(gProgram[activeProgramIndex]);
    glUniformMatrix4fv(projectionMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
    glUniformMatrix4fv(viewingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
    glUniform3fv(eyePosLoc[activeProgramIndex], 1, glm::value_ptr(eyePos));

    
    float cubeAngleRad = (float)(cubeAngle / 180.0) * M_PI;
    glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0), glm::vec3(0.0f, 0.0f, -13.0f));
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0), cubeAngleRad, glm::vec3(0.0, 1.0, 0.0));
    glm::mat4 scalingMatrix = glm::scale(glm::mat4(1.0), glm::vec3(1.0, 1.0, 1.0));
    cubeModelingMatrix = translationMatrix * rotationMatrix * scalingMatrix;
        
    
    float octaAngleRad = (float)(octaAngle / 180.0) * M_PI;
    glm::mat4 octaTranslationMatrix = glm::translate(glm::mat4(1.0), glm::vec3(0.0f, 0.0f, -13.0f));
    glm::mat4 octaRotationMatrix = glm::rotate(glm::mat4(1.0), octaAngleRad, glm::vec3(0.0, 1.0, 0.0));
    glm::mat4 octaOffsetMatrix = glm::translate(glm::mat4(1.0), glm::vec3(3.5f, 0.0f, 0.0f));
    glm::mat4 octaScalingMatrix = glm::scale(glm::mat4(1.0), glm::vec3(1.0, 1.0, 1.0));
    octaModelingMatrix = octaTranslationMatrix * octaRotationMatrix * octaOffsetMatrix * octaRotationMatrix * octaScalingMatrix;

    
    float tetraAngleRad = (float)(tetraAngle / 180.0) * M_PI;
    glm::mat4 tetraTranslationMatrix = glm::translate(glm::mat4(1.0), glm::vec3(0.0f, 0.0f, -13.0f));
    glm::mat4 tetraRotationMatrix = glm::rotate(glm::mat4(1.0), tetraAngleRad, glm::vec3(0.0, 1.0, 0.0));
    glm::mat4 tetraOffsetMatrix = glm::translate(glm::mat4(1.0), glm::vec3(1.5f, 0.0f, 0.0f));
    glm::mat4 tetraScalingMatrix = glm::scale(glm::mat4(1.0), glm::vec3(0.4, 0.4, 0.4));
    tetraModelingMatrix = octaModelingMatrix * tetraRotationMatrix * tetraOffsetMatrix * tetraRotationMatrix * tetraScalingMatrix;

    if (mode == SOLID)
    {
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(cubeModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(cubeColor));
        drawModel(cubeVAOs[cubeLvl], cubeVertexAttribBuffer[cubeLvl], cubeIndexBuffer[cubeLvl], cubeVertexDataSizeInBytes[cubeLvl], cubeFlatFaces[cubeLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(octaModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(octaColor));
        drawModel(octaVAOs[octaLvl], octaVertexAttribBuffer[octaLvl], octaIndexBuffer[octaLvl], octaVertexDataSizeInBytes[octaLvl], octaFlatFaces[octaLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(tetraModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(tetraColor));
        drawModel(tetraVAOs[tetraLvl], tetraVertexAttribBuffer[tetraLvl], tetraIndexBuffer[tetraLvl], tetraVertexDataSizeInBytes[tetraLvl], tetraFlatFaces[tetraLvl]);
    }
    else if (mode == LINE)
    {
        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(cubeModelingMatrix));
        drawModel(cubeVAOs[cubeLvl], cubeVertexAttribBuffer[cubeLvl], cubeIndexBuffer[cubeLvl], cubeVertexDataSizeInBytes[cubeLvl], cubeFlatFaces[cubeLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(octaModelingMatrix));
        drawModel(octaVAOs[octaLvl], octaVertexAttribBuffer[octaLvl], octaIndexBuffer[octaLvl], octaVertexDataSizeInBytes[octaLvl], octaFlatFaces[octaLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(tetraModelingMatrix));
        drawModel(tetraVAOs[tetraLvl], tetraVertexAttribBuffer[tetraLvl], tetraIndexBuffer[tetraLvl], tetraVertexDataSizeInBytes[tetraLvl], tetraFlatFaces[tetraLvl]);
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
        glLineWidth(2.0);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(cubeModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(cubeColor));
        drawLineModel(cubeLineVAOs[cubeLvl], cubeLineVertexAttribBuffer[cubeLvl], cubeLineIndexBuffer[cubeLvl], cubeLineVertexDataSizeInBytes[cubeLvl], cubeLineFlatFaces[cubeLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(octaModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(octaColor));
        drawLineModel(octaLineVAOs[octaLvl], octaLineVertexAttribBuffer[octaLvl], octaLineIndexBuffer[octaLvl], octaLineVertexDataSizeInBytes[octaLvl], octaLineFlatFaces[octaLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(tetraModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(tetraColor));
        drawLineModel(tetraLineVAOs[tetraLvl], tetraLineVertexAttribBuffer[tetraLvl], tetraLineIndexBuffer[tetraLvl], tetraLineVertexDataSizeInBytes[tetraLvl], tetraLineFlatFaces[tetraLvl]);
        glLineWidth(1.0);
    }
    else
    {
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(cubeModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(cubeColor));
        drawLineModel(cubeLineVAOs[cubeLvl], cubeLineVertexAttribBuffer[cubeLvl], cubeLineIndexBuffer[cubeLvl], cubeLineVertexDataSizeInBytes[cubeLvl], cubeLineFlatFaces[cubeLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(octaModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(octaColor));
        drawLineModel(octaLineVAOs[octaLvl], octaLineVertexAttribBuffer[octaLvl], octaLineIndexBuffer[octaLvl], octaLineVertexDataSizeInBytes[octaLvl], octaLineFlatFaces[octaLvl]);
        glUniformMatrix4fv(modelingMatrixLoc[activeProgramIndex], 1, GL_FALSE, glm::value_ptr(tetraModelingMatrix));
        glUniform3fv(kdLoc[activeProgramIndex], 1, glm::value_ptr(tetraColor));
        drawLineModel(tetraLineVAOs[tetraLvl], tetraLineVertexAttribBuffer[tetraLvl], tetraLineIndexBuffer[tetraLvl], tetraLineVertexDataSizeInBytes[tetraLvl], tetraLineFlatFaces[tetraLvl]);
    }
        

    if (not pause)
    {
        cubeAngle += 0.25;
        octaAngle += 0.5;
        tetraAngle += 1;
    }
}

void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);

    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glOrtho(-10, 10, -10, 10, -10, 10);
    //gluPerspective(45, 1, 1, 100);

	// Use perspective projection

	float fovyRad = (float) (45.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, 1.0f, 1.0f, 100.0f);

	// Assume default camera position and orientation (camera is at
	// (0, 0, 0) with looking at -z direction and its up vector pointing
	// at +y direction)

	viewingMatrix = glm::mat4(1);

    //glMatrixMode(GL_MODELVIEW);
    //glLoadIdentity();
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_E && action == GLFW_PRESS)
    {
        if (cubeLvl < 4)
			cubeLvl++;
        if (octaLvl < 4)
			octaLvl++;
        if (tetraLvl < 4)
            tetraLvl++;
    }
    else if (key == GLFW_KEY_T && action == GLFW_PRESS)
    {
		if (cubeLvl > 0)
            cubeLvl--;
        if (octaLvl > 0)
            octaLvl--;
        if (tetraLvl > 0)
            tetraLvl--;
    }
    else if (key == GLFW_KEY_M && action == GLFW_PRESS)
    {
        if (mode == SOLID)
            mode = LINE;
        else if (mode == LINE)
            mode = WIREFRAME;
        else
            mode = SOLID;
    }
    else if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        pause = false;
        mode = SOLID;

        cubeLvl = 1;
        cubeAngle = 0;

        octaLvl = 0;
        octaAngle = 0;

        tetraLvl = 0;
        tetraAngle = 0;
    }
    else if (key == GLFW_KEY_S && action == GLFW_PRESS)
    {
        pause = not pause;
    }
}

void mainLoop(GLFWwindow* window)
{
    while (!glfwWindowShouldClose(window))
    {
        display();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
    GLFWwindow* window;
    if (!glfwInit())
    {
        exit(-1);
    }

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    int width = 640, height = 480;
    window = glfwCreateWindow(width, height, "Simple `Example", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, width, height); // need to call this once ourselves
    mainLoop(window); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
