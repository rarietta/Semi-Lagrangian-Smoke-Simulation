#include "grid_data.h"


GridData::GridData() :
   mDfltValue(0.0), mMax(0.0,0.0,0.0)
{
}

GridData::GridData(const GridData& orig) :
   mDfltValue(orig.mDfltValue)
{
   mData = orig.mData;
   mMax = orig.mMax;
}

GridData::~GridData() 
{
}

std::vector<double>& GridData::data()
{
   return mData;
}

GridData& GridData::operator=(const GridData& orig)
{
   if (this == &orig)
   {
      return *this;
   }
   mDfltValue = orig.mDfltValue;
   mData = orig.mData;
   mMax = orig.mMax;
   return *this;
}

void GridData::initialize(double dfltValue)
{
   mDfltValue = dfltValue;
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridData::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

const double GridData::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

void GridData::getCell(const vec3& pt, int& i, int& j, int& k)
{
   vec3 pos = worldToSelf(pt); 
   i = (int) (pos[0]/theCellSize);
   j = (int) (pos[1]/theCellSize);
   k = (int) (pos[2]/theCellSize);   
}

double GridData::interpolate(const vec3& pt)
{

	// TODO: Implement sharper cubic interpolation here.


	// LINEAR INTERPOLATION:
	vec3 pos = worldToSelf(pt);

	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);
	int k = (int) (pos[2]/theCellSize);

	double scale = 1.0/theCellSize;  
	double fractx = scale*(pos[0] - i*theCellSize);
	double fracty = scale*(pos[1] - j*theCellSize);
	double fractz = scale*(pos[2] - k*theCellSize);

	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);

	// Y @ low X, low Z:
	double tmp1 = (*this)(i,j,k);
	double tmp2 = (*this)(i,j+1,k);
	// Y @ high X, low Z:
	double tmp3 = (*this)(i+1,j,k);
	double tmp4 = (*this)(i+1,j+1,k);

	// Y @ low X, high Z:
	double tmp5 = (*this)(i,j,k+1);
	double tmp6 = (*this)(i,j+1,k+1);
	// Y @ high X, high Z:
	double tmp7 = (*this)(i+1,j,k+1);
	double tmp8 = (*this)(i+1,j+1,k+1);

	// Y @ low X, low Z
	double tmp12 = LERP(tmp1, tmp2, fracty);
	// Y @ high X, low Z
	double tmp34 = LERP(tmp3, tmp4, fracty);

	// Y @ low X, high Z
	double tmp56 = LERP(tmp5, tmp6, fracty);
	// Y @ high X, high Z
	double tmp78 = LERP(tmp7, tmp8, fracty);

	// X @ low Z
	double tmp1234 = LERP (tmp12, tmp34, fractx);
	// X @ high Z
	double tmp5678 = LERP (tmp56, tmp78, fractx);

	// Z
	double tmp = LERP(tmp1234, tmp5678, fractz);
	return tmp;

}

vec3 GridData::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
   return out;
}

GridDataX::GridDataX() : GridData()
{
}

GridDataX::~GridDataX()
{
}

void GridDataX::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*(theDim[0]+1);
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize((theDim[0]+1)*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataX::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

const double GridDataX::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataX::worldToSelf(const vec3& pt) const
{   
   vec3 out;
   out[0] = min(max(0.0, pt[0]), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataY::GridDataY() : GridData()
{
}

GridDataY::~GridDataY()
{
}

void GridDataY::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*(theDim[1]+1);
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*(theDim[1]+1)*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataY::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

const double GridDataY::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataY::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataZ::GridDataZ() : GridData()
{
}

GridDataZ::~GridDataZ()
{
}

void GridDataZ::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*(theDim[2]+1);
   mData.resize(theDim[0]*theDim[1]*(theDim[2]+1), false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataZ::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

const double GridDataZ::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

vec3 GridDataZ::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]), mMax[2]);
   return out;
}
