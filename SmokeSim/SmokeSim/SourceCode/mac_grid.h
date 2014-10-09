// Modified by Peter Kutz, 2011 and 2012.

#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)

#include "open_gl_headers.h"
#include "vec.h"
#include "grid_data.h"
#include "grid_data_matrix.h"

class Camera;

class MACGrid
{

public:
	MACGrid();
	~MACGrid();
	MACGrid(const MACGrid& orig);
	MACGrid& operator=(const MACGrid& orig);

	void reset();

	void draw(const Camera& c);
	void updateSources();
	void advectVelocity(double dt);
	void addExternalForces(double dt);
	void project(double dt);
	void advectTemperature(double dt);
	void advectDensity(double dt);

protected:

	// Setup:
	void initialize();

	// Simulation:
	void computeBouyancy(double dt);
	void computeVorticityConfinement(double dt);

	// Rendering:
	struct Cube { vec3 pos; vec4 color; double dist; };
	void drawWireGrid();
	void drawSmokeCubes(const Camera& c);
	void drawSmoke(const Camera& c);
	void drawCube(const MACGrid::Cube& c);
	void drawFace(const MACGrid::Cube& c);
	void drawVelocities();
	vec4 getRenderColor(int i, int j, int k);
	vec4 getRenderColor(const vec3& pt);
	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);

	// GridData accessors:
	enum Direction { X, Y, Z };
	vec3 getVelocity(const vec3& pt);
	double getVelocityX(const vec3& pt);
	double getVelocityY(const vec3& pt);
	double getVelocityZ(const vec3& pt);
	double getTemperature(const vec3& pt);
	double getDensity(const vec3& pt);
	vec3 getCenter(int i, int j, int k);

	// Sets up the A matrix:
	void setUpAMatrix();

	// Conjugate gradient stuff:
	bool conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);
	bool isValidCell(int i, int j, int k);
	bool checkDivergence(double threshold);

	// Fluid grid cell properties:
	GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
	GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
	GridDataZ mW; // W component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
	GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
	GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
	GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ

	// The A matrix:
	GridDataMatrix AMatrix;

public:

	enum RenderMode { CUBES, SHEETS };
	static RenderMode theRenderMode;
	static bool theDisplayVel;
	
	// Saves smoke in CIS 460 volumetric format:
	void saveSmoke(const char* fileName);

};

#endif