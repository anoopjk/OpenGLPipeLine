/**
* minigl.cpp
* -------------------------------
* Implement miniGL here.
* Do not use any additional files
*/
// #includes of all the necessary libraries
#include <iostream>
#include <algorithm>
#include <vector>
#include <stack>
#include <cstdio>
#include <cstdlib>


// including the minigl header file
#include "minigl.h"


//declaring global variables for matrix row and column size
const unsigned int  mrow_size = 4;
const unsigned int  mcol_size = 4;
#define PI 3.14159265
//namespace of std
using namespace std;

///////////////////////////////////////////////////////////
//declarations of global variables
MGLpixel pixel_RGB[3];

//declaring the classes required for the operations
class mMatrix;
class mVertex;
class mPixel;

int mgl_shape;
int mgl_matrix;




// class definitions
class mMatrix
{

public:
	MGLfloat M[4][4];
	
	////////////////////////////////////////////////////////
	// matrix constructors
	mMatrix()
	{
		zero_matrix();
		initialize_matrix(0, 0, 0);
	}
	/////////////////////////////////////////////////////////
	mMatrix(MGLfloat X, MGLfloat Y, MGLfloat Z)
	{
		zero_matrix();
		initialize_matrix(X, Y, Z);
	}
	////////////////////////////////////////////
	// constructing a zero matrix 
	void zero_matrix()
	{
		for (unsigned int i = 0; i < mrow_size; ++i)
		{
			for (unsigned int j = 0; j < mcol_size; ++j)
			{
				M[i][j] = 0.0;
			}
		}
	}
	void initialize_matrix(MGLfloat X, MGLfloat Y, MGLfloat Z)
	{
		
		M[3][3] = 1; //w

		// set the x, y, z coordinate
		M[0][3] = X;
		M[1][3] = Y;
		M[2][3] = Z;
	}
	////////////////////////////////////////////////////////
	//Assignement operator overloading
	mMatrix& operator=(const mMatrix& rhs)
	{
		if (this != &rhs)
		{
			for (unsigned int i = 0; i < mrow_size; ++i)
			{
				for (unsigned int j = 0; j < mcol_size; ++j)
				{
					M[i][j] = rhs.M[i][j];
				}
			}
		}

		return *this;
	}
	/////////////////////////////////////////////////////////////
	// Currently unused. mglMultMatrix is used instead.
	//matrix multiplication overloaded*
	mMatrix operator*(const mMatrix& rhs)
	{
		mMatrix result;
		result.zero_matrix();

		for (unsigned int i = 0; i < mrow_size; ++i)
		{
			for (unsigned int j = 0; j < mcol_size; ++j)
			{
				for (unsigned int k = 0; k < mcol_size; ++k)
				{
					result.M[i][j] += M[i][k] * rhs.M[k][j];
				}

			}
		}

		return result;
	}
	/////////////////////////////////////////////////////////////////////////
	// outputting a matrix, overloading << operator
	friend ostream& operator<<(ostream& out, const mMatrix& matrix)
	{
		for (unsigned int i = 0; i < mrow_size; ++i)
		{
			for (unsigned int j = 0; j < mcol_size; ++j)
			{
				out << matrix.M[i][j] << " ";
			}
			out << "\n";
		}
		return out;

	}
	///////////////////////////////////////////////////////////////////////////
	//setup a scaling matrix

	void scale_matrix(float X, float Y, float Z)
	{
		zero_matrix();
		//setting the main diagonal part
		M[0][0] = X;
		M[1][1] = Y;
		M[2][2] = Z;
		M[3][3] = 1;
	}
	////////////////////////////////////////////////////////////////////////////
	//translation operation
	void translate_matrix(MGLfloat X, MGLfloat Y, MGLfloat Z)
	{
		zero_matrix();
		
		M[0][0] = 1;
		M[1][1] = 1;
		M[2][2] = 1;
		M[3][3] = 1;
		// setting the translation part
		M[0][3] = X;
		M[1][3] = Y;
		M[2][3] = Z;
	}
	////////////////////////////////////////////////////////////////////////////

		
};

///matrix stacks declaration
stack<mMatrix> Model_stack;
stack<mMatrix> Projection_stack;


///////////////////////////////////////////////////////////////
class mPixel
{
public:
	int x, y;
	MGLpixel pixel_color;
	MGLfloat z;

	mPixel(int X, int Y, MGLpixel c, MGLfloat Z)
	{
		x = X;
		y = Y;
		pixel_color = c;
		z = Z;
	}
	
};
///////////////////////////////////////
vector<mPixel> frame_buffer;
vector<mPixel> z_buffer;
////////////////////////////////////////////////////////////////////
//class mVertex
class mVertex
{
public:
	MGLfloat x, y, z, w;
	MGLfloat x_image, y_image, z_image, w_image;
	MGLpixel vertex_color[3];

	mVertex()
	{
		x = 0;
		y = 0;
		z = 0;
		w = 1;

		for (unsigned int i = 0; i < 3; ++i)
		{
			vertex_color[i] = pixel_RGB[i];
		}
	}

	mVertex(MGLfloat X, MGLfloat Y, MGLfloat Z, MGLfloat W)
	{
		x = X;
		y = Y;
		z = Z;
		w = W;

		for (unsigned int i = 0; i < 3; ++i)
		{
			vertex_color[i] = pixel_RGB[i];
		}

	}

	////////////////////////////////////////////////
	//overloaded assignement operator
	mVertex& operator=(const  mVertex& rhs)
	{
		if (this != &rhs)
		{
			x = rhs.x;
			y = rhs.y;
			z = rhs.z;
			w = rhs.w;

			x_image = rhs.x_image;
			y_image = rhs.y_image;
			z_image = rhs.z_image;
			w_image = rhs.w_image;


			for (unsigned int i = 0; i < 3; ++i)
			{
				vertex_color[i] = rhs.vertex_color[i];
			}

		}
		return *this;
	}

	///////////////////////////////////////////////////////
	//vector-matrix multiplication (v*M)
	mVertex operator*(const mMatrix& rhs)
	{
		MGLfloat v[4] = { x, y, z, w };
		MGLfloat result_vector[4] = { 0,0,0,0 };

		for (unsigned int i = 0; i < mrow_size; ++i)
		{
			for (unsigned int j = 0; j < mcol_size; ++j)
			{
				result_vector[i] += v[j] * rhs.M[i][j];
			}
		}
		
		mVertex result(result_vector[0], result_vector[1], result_vector[2], result_vector[3]);
		return result;
	}

	void world2screen(MGLsize width, MGLsize height)
	{
		MGLfloat xtilde = (x*width) / 2;
		MGLfloat ytilde = (y*height) / 2;

	//	perspective division

			x = xtilde / w;
			y = ytilde / w;
			z = z / w;
			w = 1;

			if (x > width)
			{
				x = width;
			}

			x_image = x;
			y_image = y;
			z_image = z;
			w_image = w;
	}

	////////////////////////////////////////////////

	void transform()
	{
		mMatrix model = Model_stack.top();
		mMatrix project = Projection_stack.top();

		mMatrix Mtrans;
		Mtrans.translate_matrix(1, 1, 1);

		mVertex V(x, y, z, w);
		// conversion is in the order ,
		// Modelview -> projection -> translation -> updating
		//operations represented in vertex*matrix multiplication
		V = V*model;
		V = V*project;
		V = V*Mtrans;

		x = V.x;
		y = V.y;
		z = V.z;
		w = V.w;

	}
};
//////////////////////////////////////////////
vector< vector<mVertex> > geometry_list;
vector<mVertex> vertex_list;

////////////////////////////////////////////////////////////////
void set_pixel(int x, int y, MGLpixel c, MGLfloat z)
{
	mPixel P(x, y, c, z);
	frame_buffer.push_back(P);
	z_buffer.push_back(P);
}

//////////////////////////////////////////////////////////
vector<float>  BB_cord(const vector<mVertex>& v)
{
	vector<float> result;
	unsigned int num_vertices = v.size();
	float xmin = v[0].x_image;
	float ymin = v[0].y_image;
	float xmax = v[0].x_image;
	float ymax = v[0].y_image;
	
	for (unsigned int i = 1; i < num_vertices; ++i)
	{
		if (v[i].x_image < xmin)
			xmin = v[i].x_image;

		if (v[i].y_image < ymin)
			ymin = v[i].y_image;

		if (v[i].x_image > xmax)
			xmax = v[i].x_image;

		if (v[i].y_image > ymax)
			ymax = v[i].y_image;
	}
	result.push_back(xmin);
	result.push_back(ymin);
	result.push_back(xmax);
	result.push_back(ymax);
	
	return result;
}


///////////////////////////////////////////////////
void world2screen_list(MGLsize width, MGLsize height, vector<mVertex>& vl)
{
	for (unsigned int i = 0; i < vl.size(); ++i)
	{
		vl[i].world2screen(width, height);
	}
}

///////////////////////////////////////////////////
//line equation
float line(float x, float y,float x0, float y0, float x1,float y1)
{
	float A = y0 - y1;
	float B = x1 - x0;
	float C = (x0*y1) - (x1*y0);

	//cout << "line" << A*x + B*y + C << endl;
	return A*x + B*y + C;
}

///////////////////////////////////////////
//gouraud shading
MGLpixel gouraud(float alpha, float beta, float gamma, const MGLpixel* c0, const MGLpixel* c1, const MGLpixel* c2)
{
	MGLfloat cr = alpha*c0[0] + beta*c1[0] + gamma*c2[0];
	MGLfloat cg = alpha*c0[1] + beta*c1[1] + gamma*c2[1];
	MGLfloat cb = alpha*c0[2] + beta*c1[2] + gamma*c2[2];

	MGLpixel result = 0;
	MGL_SET_RED(result, (MGLpixel)cr);
	MGL_SET_GREEN(result, (MGLpixel)cg);
	MGL_SET_BLUE(result, (MGLpixel)cb);

	return result;

}

//Drawing a triangle based on barycentric coordinates
void draw_triangle(float x, float y, const mVertex& a, const mVertex& b, const mVertex& c)
{
	float xa = a.x_image;
	float ya = a.y_image;

	float xb = b.x_image;
	float yb = b.y_image;

	float xc = c.x_image;
	float yc = c.y_image;

	float alpha = line(x, y, xb, yb, xc, yc) / line(xa, ya, xb, yb, xc, yc);
	float beta = line(x, y, xc, yc, xa, ya) / line(xb, yb, xc, yc, xa , ya);
	float gamma = line(x, y, xa, ya, xb, yb) / line(xc, yc, xa, ya, xb, yb);

	if ((alpha >= 0 && alpha < 1) && (beta >= 0 && beta < 1) && (gamma >= 0 && gamma < 1))
	{

		MGLpixel color = gouraud(alpha, beta, gamma, a.vertex_color, b.vertex_color, c.vertex_color);
		MGLfloat Z = alpha*a.z_image + beta*b.z_image + gamma*c.z_image;

		set_pixel((int)x, (int)y, color, Z);
	}
	
}
//////////////////////////////////////////////////////
void rasterize_triangle(MGLsize width, MGLsize height, const vector<mVertex>& v)
{
	vector<float> BBox = BB_cord(v);
	float xmin = BBox.at(0);
	float ymin = BBox.at(1);
	float xmax = BBox.at(2);
	float ymax = BBox.at(3);

	cout << "bbox" << xmax << ' ' << ymax << endl;
	for (float x = xmin; x <= xmax; ++x)
	{
		for (float y = ymin; y <= ymax; ++y)
		{
			draw_triangle(x, y, v.at(0), v.at(1), v.at(2));
		}
	}
}

//////////////////////////////////////////////////////////
void rasterize_quad(MGLsize width, MGLsize height, const vector<mVertex>& v)
{
	vector<float> BBox = BB_cord(v);
	float xmin = BBox.at(0);
	float ymin = BBox.at(1);
	float xmax = BBox.at(2);
	float ymax = BBox.at(3);


	for (float x = xmin; x <= xmax; ++x)
	{
		for (float y = ymin; y <= ymax; ++y)
		{
			draw_triangle(x, y, v.at(0), v.at(1), v.at(2));
			draw_triangle(x, y, v.at(0), v.at(2), v.at(3));
		}
	}
}


/**
* Standard macro to report errors
*/
inline void MGL_ERROR(const char* description) {
	printf("%s\n", description);
	exit(1);
}

bool z_sort(mPixel a, mPixel b)
{
	return a.z > b.z; //check condition for sorting the pixels based on z
}
/**
* Read pixel data starting with the pixel at coordinates
* (0, 0), up to (width,  height), into the array
* pointed to by data.  The boundaries are lower-inclusive,
* that is, a call with width = height = 1 would just read
* the pixel at (0, 0).
*
* Rasterization and z-buffering should be performed when
* this function is called, so that the data array is filled
* with the actual pixel values that should be displayed on
* the two-dimensional screen.
*/
void mglReadPixels(MGLsize width,
	MGLsize height,
	MGLpixel *data)
{
	for (unsigned int i = 0; i < geometry_list.size(); ++i)
	{
		world2screen_list(width, height, geometry_list.at(i));
		
		MGLsize vertex_count = geometry_list.at(i).size();

		if (vertex_count == 3)
		{
			rasterize_triangle(width, height, geometry_list.at(i));
		}
		else if (vertex_count == 4)
		{
			rasterize_quad(width, height, geometry_list.at(i));
		}
	}

	//sort function http://www.cplusplus.com/reference/algorithm/sort/
	sort(z_buffer.begin(), z_buffer.end(), z_sort);//sorting based on the z element of each pixel

	for (unsigned int i = 0; i < z_buffer.size(); ++i)
	{
		int x = z_buffer.at(i).x;
		int y = z_buffer.at(i).y;

		MGLpixel color = z_buffer.at(i).pixel_color;
		data[y*width + x] = color; // all the pixel coordinates are represented by this data array 
	}
	geometry_list.clear();
}

/**
* Start specifying the vertices for a group of primitives,
* whose type is specified by the given mode.
*/
void mglBegin(MGLpoly_mode mode)
{
	if (mode == MGL_TRIANGLES || mode == MGL_QUADS)
		mgl_shape = mode;

	else
		mgl_shape = -1;
}

/**
* Stop specifying the vertices for a group of primitives.
*/
void mglEnd()
{
	if (!vertex_list.empty())
	{
		geometry_list.push_back(vertex_list);
		vertex_list.clear();
	}
}

/**
* Specify a two-dimensional vertex; the x- and y-coordinates
* are explicitly specified, while the z-coordinate is assumed
* to be zero.  Must appear between calls to mglBegin() and
* mglEnd().
*/
void mglVertex2(MGLfloat x,
	MGLfloat y)
{
	if (mgl_shape == -1)
	{
		MGL_ERROR("Error: Missing mglBegin! Aborting mission.\n");
		exit(1);
	}
	if (mgl_shape == MGL_TRIANGLES && vertex_list.size() == 3)
	{
		geometry_list.push_back(vertex_list);
		vertex_list.clear();
	}
	if (mgl_shape == MGL_QUADS && vertex_list.size() == 4)
	{
		geometry_list.push_back(vertex_list);
		vertex_list.clear();
	}

	mVertex v(x, y, 0, 1);

	v.transform();
	vertex_list.push_back(v);
}

/**
* Specify a three-dimensional vertex.  Must appear between
* calls to mglBegin() and mglEnd().
*/
void mglVertex3(MGLfloat x,
	MGLfloat y,
	MGLfloat z)
{
	if (mgl_shape == -1)
	{
		MGL_ERROR("Error: Missing mglBegin! Aborting mission.\n");
		exit(1);
	}
	if (mgl_shape == MGL_TRIANGLES && vertex_list.size() == 3)
	{
		geometry_list.push_back(vertex_list);
		vertex_list.clear();
	}
	if (mgl_shape == MGL_QUADS && vertex_list.size() == 4)
	{
		geometry_list.push_back(vertex_list);
		vertex_list.clear();
	}

	mVertex v(x, y, z, 1);

	v.transform();
	vertex_list.push_back(v);
}

/**
* Set the current matrix mode (modelview or projection).
*/
void mglMatrixMode(MGLmatrix_mode mode)
{
	if (mode == MGL_MODELVIEW || mode == MGL_PROJECTION)
		mgl_matrix = mode;
	else
		mgl_matrix = -1;
}

/**
* Push a copy of the current matrix onto the stack for the
* current matrix mode.
*/
void mglPushMatrix()
{
	if (mgl_matrix == MGL_MODELVIEW && !Model_stack.empty())
		Model_stack.push(Model_stack.top());

	else if (mgl_matrix == MGL_PROJECTION && !Projection_stack.empty())
		Projection_stack.push(Projection_stack.top());

}

/**
* Pop the top matrix from the stack for the current matrix
* mode.
*/
void mglPopMatrix()
{
	if (mgl_matrix == MGL_MODELVIEW && !Model_stack.empty())
	{
		if (!Model_stack.empty())
			Model_stack.pop();
	}

	else if (mgl_matrix == MGL_PROJECTION && !Projection_stack.empty())
	{
		if (!Projection_stack.empty())
			Projection_stack.pop();
	}
}

/**
* Replace the current matrix with the identity.
*/
void mglLoadIdentity()
{
	mMatrix Identity;
	Identity.zero_matrix();
	//cout << Identity.M;
	Identity.M[0][0] = 1;
	
	Identity.M[1][1] = 1;
	Identity.M[2][2] = 1;
	Identity.M[3][3] = 1;
	//cout << Identity;
	if (mgl_matrix == MGL_PROJECTION)
	{

		if (!Projection_stack.empty())
		{
			Projection_stack.pop();
		}
		Projection_stack.push(Identity);

	}
	else if (mgl_matrix == MGL_MODELVIEW)
	{
		if (!Model_stack.empty())
		{
			Model_stack.pop();
		}
		Model_stack.push(Identity);
	}

}

/**
* Replace the current matrix with an arbitrary 4x4 matrix,
* specified in column-major order.  That is, the matrix
* is stored as:
*
*   ( a0  a4  a8  a12 )
*   ( a1  a5  a9  a13 )
*   ( a2  a6  a10 a14 )
*   ( a3  a7  a11 a15 )
*
* where ai is the i'th entry of the array.
*/
void mglLoadMatrix(const MGLfloat *matrix)
{
	mMatrix m;
	for (MGLsize k = 0; k < 16; ++k)
	{
		for (MGLsize i = 0; i < mrow_size; ++i)
		{
			for (MGLsize j = 0; j < mcol_size; ++j)
			{
				m.M[i][j] = matrix[k];
			}
		}

	}

	if (mgl_matrix == MGL_MODELVIEW && !Model_stack.empty())
		Model_stack.push(m);

	else if (mgl_matrix == MGL_PROJECTION && !Projection_stack.empty())
		Projection_stack.push(m);
}

/**
* Multiply the current matrix by an arbitrary 4x4 matrix,
* specified in column-major order.  That is, the matrix
* is stored as:
*
*   ( a0  a4  a8  a12 )
*   ( a1  a5  a9  a13 )
*   ( a2  a6  a10 a14 )
*   ( a3  a7  a11 a15 )
*
* where ai is the i'th entry of the array.
*/
void mglMultMatrix(const MGLfloat *matrix)
{
	mMatrix m;
	for (MGLsize k = 0; k < 16; ++k)
	{
		for (MGLsize i = 0; i < mrow_size; ++i)
		{
			for (MGLsize j = 0; j < mcol_size; ++j)
			{
				m.M[i][j] = matrix[k];
			}
		}
	}

	if (mgl_matrix == MGL_MODELVIEW && !Model_stack.empty())
		Model_stack.top() = Model_stack.top()*m;

	else if (mgl_matrix == MGL_PROJECTION && !Projection_stack.empty())
		Projection_stack.top() = Projection_stack.top()*m;

}

/**
* Multiply the current matrix by the translation matrix
* for the translation vector given by (x, y, z).
*/
void mglTranslate(MGLfloat x,
	MGLfloat y,
	MGLfloat z)
{
	mMatrix Mtrans;
	Mtrans.translate_matrix(x, y, z);

	if (mgl_matrix == MGL_PROJECTION){
		Projection_stack.top() = Projection_stack.top()*Mtrans;
		
	}
	else if (mgl_matrix == MGL_MODELVIEW){
		Model_stack.top() = Model_stack.top()*Mtrans;
		
	}
}

/**
* Multiply the current matrix by the rotation matrix
* for a rotation of (angle) degrees about the vector
* from the origin to the point (x, y, z).
*/
//the equations are given in the following link
//reference https://www.opengl.org/sdk/docs/man2/xhtml/glFrustum.xml

void mglRotate(MGLfloat angle,
	MGLfloat x,
	MGLfloat y,
	MGLfloat z)
{

	MGLfloat s = sin(angle * PI / 180);//converting into radians
	MGLfloat c = cos(angle * PI / 180);

	MGLfloat magnitude = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	if (magnitude != 0)
	{
		x = x / magnitude;
		y = y / magnitude;
		z = z / magnitude;
	}
	mMatrix R(0, 0, 0);

	R.M[0][0] = x*x * (1 - c) + c;
	R.M[0][1] = x*y * (1 - c) - z*s;
	R.M[0][2] = x*z * (1 - c) + y*s;

	R.M[1][0] = y*x * (1 - c) + z*s;
	R.M[1][1] = y*y * (1 - c) + c;
	R.M[1][2] = y*z * (1 - c) - x*s;

	R.M[2][0] = x*z * (1 - c) - y*s;
	R.M[2][1] = y*z * (1 - c) + x*s;
	R.M[2][2] = z*z * (1 - c) + c;

	if (mgl_matrix == MGL_PROJECTION)
		Projection_stack.top() = Projection_stack.top()*R;
	
	else if (mgl_matrix == MGL_MODELVIEW)
		Model_stack.top() = Model_stack.top()*R;
	
}

/**
* Multiply the current matrix by the scale matrix
* for the given scale factors.
*/
void mglScale(MGLfloat x,
	MGLfloat y,
	MGLfloat z)
{
	mMatrix Mscale;
	Mscale.scale_matrix(x, y, z);

	if (mgl_matrix == MGL_PROJECTION)
		Projection_stack.top() = Projection_stack.top()*Mscale;
	
	else if (mgl_matrix == MGL_MODELVIEW)
		Model_stack.top() = Model_stack.top()*Mscale;
	

}

/**
* Multiply the current matrix by the perspective matrix
* with the given clipping plane coordinates.
*/
//http://www.songho.ca/opengl/gl_transform.html
void mglFrustum(MGLfloat left,
	MGLfloat right,
	MGLfloat bottom,
	MGLfloat top,
	MGLfloat near,
	MGLfloat far)
{
	
	mMatrix Mfrustum;

	Mfrustum.M[0][0] = (2 * near) / (right - left);
	Mfrustum.M[1][1] = (2 * near) / (top - bottom);


	Mfrustum.M[0][2] = (right + left)/(right-left);
	Mfrustum.M[1][2] = (top + bottom) / (top - bottom);
	Mfrustum.M[2][2] = -(far + near) / (far - near);
	Mfrustum.M[2][3] = -(2 * far*near) / (far - near);
	Mfrustum.M[3][2] = -1;
	Mfrustum.M[3][3] = 0;

	if (mgl_matrix == MGL_PROJECTION)
		Projection_stack.top() = Projection_stack.top()*Mfrustum;
	
	else if (mgl_matrix == MGL_MODELVIEW)
		Model_stack.top() = Model_stack.top()*Mfrustum;
	

}

/**
* Multiply the current matrix by the orthographic matrix
* with the given clipping plane coordinates.
*///http://www.songho.ca/opengl/gl_transform.html

void mglOrtho(MGLfloat left,
	MGLfloat right,
	MGLfloat bottom,
	MGLfloat top,
	MGLfloat near,
	MGLfloat far)
{
	
	mMatrix Mortho;
	Mortho.M[0][0] = 2 / (right - left);
	Mortho.M[0][3] = -(right + left) / (right - left);
	
	Mortho.M[1][1] = 2 / (top - bottom);
	Mortho.M[1][3] = -(top + bottom) / (top - bottom);
	Mortho.M[2][2] = -2 / (far - near);

	Mortho.M[0][3] = -(right + left) / (right - left);
	
	Mortho.M[2][3] = -(far + near) / (far - near);

	if (mgl_matrix == MGL_PROJECTION)
		Projection_stack.top() = Projection_stack.top()*Mortho;
	
	else if (mgl_matrix == MGL_MODELVIEW)
		Model_stack.top() = Model_stack.top()*Mortho;
	
}

/**
* Set the current color for drawn shapes.
*/
void mglColor(MGLbyte red,
	MGLbyte green,
	MGLbyte blue)
{
	//I am assigning this colours in the class mVertex
	pixel_RGB[0] = red;
	pixel_RGB[1] = green;
	pixel_RGB[2] = blue;
}
