// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>


// Globals:
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_X \
	for(int k = 0; k < theDim[MACGrid::Z]; k++) \
		for(int j = 0; j < theDim[MACGrid::Y]; j++) \
			for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_Y \
   for(int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Z \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 



MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   setUpAMatrix();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO: Set initial values for density, temperature, and velocity.

	// create a source at cell
	mD(20,1,15) = 10.0;
	mT(20,1,15) = 150.0;
	mV(20,1,15) = 20.0;
}

void MACGrid::advectVelocity(double dt)
{
	// TODO: Calculate new velocities and store in target.
	FOR_EACH_FACE_X {
		vec3 uPt = getCenter(i,j,k) - 0.5 * vec3(theCellSize,0,0);
		vec3 uVel = getVelocity(uPt);
		vec3 uSrc = uPt - uVel*dt;
		float uAdvected;
		int a, b, c;
		mU.getCell(uSrc, a, b, c);
		if (!isValidCell(a, b, c)) {
			if (a < 0) a = 0;
			if (a > theDim[X]) a = theDim[X];
			if (b < 0) b = 0;
			if (b > theDim[Y]) b = theDim[Y];
			if (c < 0) c = 0;
			if (c > theDim[Z]) c = theDim[Z];
		} uAdvected = mU.interpolate(uSrc);
		target.mU(i,j,k) = uAdvected;
	}

	FOR_EACH_FACE_Y {
		vec3 vPt = getCenter(i,j,k) - 0.5 * vec3(0,theCellSize,0);
		vec3 vVel = getVelocity(vPt);
		vec3 vSrc = vPt - vVel*dt;
		float vAdvected;
		int a, b, c;
		mV.getCell(vSrc, a, b, c);
		if (!isValidCell(a, b, c)) {
			if (a < 0) a = 0;
			if (a > theDim[X]) a = theDim[X];
			if (b < 0) b = 0;
			if (b > theDim[Y]) b = theDim[Y];
			if (c < 0) c = 0;
			if (c > theDim[Z]) c = theDim[Z];
		} vAdvected = mV.interpolate(vSrc);
		target.mV(i,j,k) = vAdvected;
	}

	FOR_EACH_FACE_Z {
		vec3 wPt = getCenter(i,j,k) - 0.5 * vec3(0,0,theCellSize);
		vec3 wVel = getVelocity(wPt);
		vec3 wSrc = wPt - wVel*dt;
		float wAdvected;
		int a, b, c;
		mW.getCell(wSrc, a, b, c);
		if (!isValidCell(a, b, c)) {
			if (a < 0) a = 0;
			if (a > theDim[X]) a = theDim[X];
			if (b < 0) b = 0;
			if (b > theDim[Y]) b = theDim[Y];
			if (c < 0) c = 0;
			if (c > theDim[Z]) c = theDim[Z];
		} wAdvected = mW.interpolate(wSrc);
		target.mW(i,j,k) = wAdvected;
	}
 
    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
	FOR_EACH_CELL {

		// get cell point at cell center
		vec3 cellPt = getCenter(i,j,k);

		// find velocity at cell point
		vec3 cellVel = getVelocity(cellPt);

		// track back from cell point using velocity
		vec3 cellSrc = cellPt - cellVel*dt;

		// find velocity at source
		double tempAdvected;
		int a, b, c;
		mT.getCell(cellSrc, a, b, c);
		if (!isValidCell(a, b, c)) tempAdvected = mT(i,j,k);
		else					   tempAdvected = mT.interpolate(cellSrc);

		// set advected velocities
		target.mT(i,j,k) = tempAdvected;
	}

    // Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new density and store in target.
	FOR_EACH_CELL {

		// get cell point at cell center
		vec3 cellPt = target.getCenter(i,j,k);

		// find velocity at cell point
		vec3 cellVel = getVelocity(cellPt);

		// track back from cell point using velocity
		vec3 cellSrc = cellPt - cellVel*dt;

		// find velocity at source
		double densityAdvected;
		int a, b, c;
		mD.getCell(cellSrc, a, b, c);
		if (!isValidCell(a, b, c)) densityAdvected = mD(i,j,k);
		else					   densityAdvected = mD.interpolate(cellSrc);

		// set advected velocities
		target.mD(i,j,k) = densityAdvected;
	}

    // Then save the result to our object.
    mD = target.mD;
}

void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target.
	float alpha = 0.1;
	float beta = 0.3;
	FOR_EACH_CELL {
		if (i == 0 || j == 0 || k == 0) continue;
		float avgTemp = (mT(i,j,k) + mT(i,j-1,k)) / 2.0f;
		float avgDens = (mD(i,j,k) + mD(i,j-1,k)) / 2.0f;
		target.mV(i,j,k) = mV(i,j,k) - alpha*avgDens + beta*avgTemp;
	}
	// Then save the result to our object.
	mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	// Apply the forces to the current velocity and store the result in target.

	// compute vorticity at each cell center
	GridData vortX = GridData(); vortX.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	GridData vortY = GridData(); vortY.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	GridData vortZ = GridData();vortZ.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	FOR_EACH_CELL {
		vortX(i,j,k) = (mW(i,j+1,k)-mW(i,j,k)) / theCellSize - (mV(i,j,k+1)-mV(i,j,k)) / theCellSize;
		vortY(i,j,k) = (mU(i,j,k+1)-mU(i,j,k)) / theCellSize - (mW(i+1,j,k)-mW(i,j,k)) / theCellSize;
		vortZ(i,j,k) = (mV(i+1,j,k)-mV(i,j,k)) / theCellSize - (mU(i,j+1,k)-mU(i,j,k)) / theCellSize;
	}

	// compute vorticity confinement force at each center based on vorticity
	GridData fconf_centerX = GridData(); fconf_centerX.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	GridData fconf_centerY = GridData(); fconf_centerY.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	GridData fconf_centerZ = GridData(); fconf_centerZ.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	fconf_centerX.initialize();
	fconf_centerY.initialize();
	fconf_centerZ.initialize();
	FOR_EACH_CELL {

		// ignore boundary cells
		if (i == 0 || j == 0 || k == 0) continue;
		if (i == theDim[X]-1 || j == theDim[Y]-1 || k == theDim[Z]-1) continue;

		// compute gradient of vorticity using central differences
		float p, q;
		p = vec3(vortX(i+1,j,k), vortY(i+1,j,k), vortZ(i+1,j,k)).Length();
		q = vec3(vortX(i-1,j,k), vortY(i-1,j,k), vortZ(i-1,j,k)).Length(); 
		float grad1 = (p - q) / (2.0*theCellSize);

		p = vec3(vortX(i,j+1,k), vortY(i,j+1,k), vortZ(i,j+1,k)).Length();
		q = vec3(vortX(i,j-1,k), vortY(i,j-1,k), vortZ(i,j-1,k)).Length(); 
		float grad2 = (p - q) / (2.0*theCellSize);

		p = vec3(vortX(i,j,k+1), vortY(i+1,j,k+1), vortZ(i+1,j,k+1)).Length();
		q = vec3(vortX(i,j,k-1), vortY(i-1,j,k-1), vortZ(i-1,j,k-1)).Length(); 
		float grad3 = (p - q) / (2.0*theCellSize);

		vec3 gradVort(grad1, grad2, grad3);

		// compute "N" vector
		vec3 N_ijk = gradVort / (gradVort.Length() + 10e-20);

		// add this to the velocities
		float epsilon = 0.55;
		vec3 vorticity = vec3(vortX(i,j,k), vortY(i,j,k), vortZ(i,j,k));
		vec3 fconf_center = epsilon * theCellSize * vorticity.Cross(N_ijk);
		fconf_centerX(i,j,k) = fconf_center[X];
		fconf_centerY(i,j,k) = fconf_center[Y];
		fconf_centerZ(i,j,k) = fconf_center[Z];
	}

	// find appropriate averages at each face
	FOR_EACH_FACE_X{
		if (i == 0 || i == theDim[X]) continue;
		target.mU(i,j,k) += (fconf_centerX(i,j,k) + fconf_centerX(i-1,j,k)) / 2.0f;
	}
	FOR_EACH_FACE_Y{
		if (j == 0 || j == theDim[Y]) continue;
		target.mV(i,j,k) += (fconf_centerY(i,j,k) + fconf_centerY(i,j-1,k)) / 2.0f;
	}
	FOR_EACH_FACE_Z{
		if (k == 0 || k == theDim[Z]) continue;
		target.mW(i,j,k) += (fconf_centerZ(i,j,k) + fconf_centerZ(i,j,k-1)) / 2.0f;
	}

	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
	// boundary constraints
	FOR_EACH_FACE_X {
		if (i == 0 || i == theDim[X]+1)
			mU(i,j,k) = 0;}
	FOR_EACH_FACE_Y {
		if (j == 0 || j == theDim[Y]+1)
			mV(i,j,k) = 0;}
	FOR_EACH_FACE_Z {
		if (k == 0 || k == theDim[Z]+1)
			mW(i,j,k) = 0;}

	// TODO: Solve Ap = d for pressure.
	// 1. Construct d
	float rho = 1.0;
	float scale = -theCellSize*theCellSize*rho/dt;
	GridData d = GridData();
	d.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	FOR_EACH_CELL {
		float dU = (mU(i+1,j,k) - mU(i,j,k)) / theCellSize;
		float dV = (mV(i,j+1,k) - mV(i,j,k)) / theCellSize;
		float dW = (mW(i,j,k+1) - mW(i,j,k)) / theCellSize;
		float dValue = scale * (dU + dV + dW);
		d(i,j,k) = dValue;
	}

	// 2. Construct A
	setUpAMatrix();

	// 3. Solve for p
	GridData p = GridData();
	p.data().resize(theDim[X]*theDim[Y]*theDim[Z]);
	bool result = conjugateGradient(AMatrix, p, d, 1000, 0.01);

	// Subtract pressure from our velocity and save in target.
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	FOR_EACH_CELL {
		if (i > 0) target.mU(i,j,k) -= (dt / rho) * (p(i,j,k) - p(i-1,j,k)) / theCellSize;
		if (j > 0) target.mV(i,j,k) -= (dt / rho) * (p(i,j,k) - p(i,j-1,k)) / theCellSize;
		if (k > 0) target.mW(i,j,k) -= (dt / rho) * (p(i,j,k) - p(i,j,k-1)) / theCellSize;
		target.mP(i,j,k) = p(i,j,k);
	}

	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
	
	// boundary constraints (again)
	FOR_EACH_FACE_X {
		if (i == 0 || i == theDim[X]+1)
			mU(i,j,k) = 0;}
	FOR_EACH_FACE_Y {
		if (j == 0 || j == theDim[Y]+1)
			mV(i,j,k) = 0;}
	FOR_EACH_FACE_Z {
		if (k == 0 || k == theDim[Z]+1)
			mW(i,j,k) = 0;}

	// IMPLEMENT THIS AS A SANITY CHECK:
	//assert (checkDivergence(0.1));
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}







/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:
	z = r;

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		// For now, just bypass the preconditioner:
		z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}


bool MACGrid::checkDivergence(double threshold) {

	FOR_EACH_CELL {
		float dU = (mU(i+1,j,k) - mU(i,j,k)) / theCellSize;
		float dV = (mV(i,j+1,k) - mV(i,j,k)) / theCellSize;
		float dW = (mW(i,j,k+1) - mW(i,j,k)) / theCellSize;
		float divergence = dU + dV + dW;
		if (divergence > threshold)
			return false;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}





/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
    return vec4(1.0, 1.0, 1.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
