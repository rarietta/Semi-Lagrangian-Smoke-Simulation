// ==========================================================================
// Copyright (C) 2008 Aline Normoyle
// ==========================================================================

#ifndef camera_H_
#define camera_H_

#include "open_gl_headers.h"
#include "vec.h"
#include "matrix.h"

class Camera
{
public:
   Camera();
   virtual ~Camera();

   // Draw projection and eyepoint
   virtual void draw();

   // Print eyepoint position and basis
   virtual void print();

   // Initialize the camera with glyLookAt parameters
   virtual void set(const vec3& eyepos, const vec3& look, const vec3& up);

   // Get camera state
   virtual void setPosition(const vec3& pos);
   virtual const vec3& getPosition() const;
   virtual const vec3& getUp() const;
   virtual const vec3& getBackward() const;
   virtual const vec3& getRight() const;
   virtual vec3 getRelativePosition(float left, float up, float forward);
   virtual float heading() const;
   virtual float pitch() const;

   // Camera frustrum managements
   virtual void setProjection(
      float vfov, float aspect, float zNear, float zFar);
   virtual void getProjection(
      float* vfov, float* aspect, float* zNear, float* zFar);
   virtual void getViewport(int& x, int& y, int& w, int& h);

   // Relative movement commands
   virtual void moveLeft(int scale = 1.0);
   virtual void moveRight(int scale = 1.0);
   virtual void moveUp(int scale = 1.0);
   virtual void moveDown(int scale = 1.0);
   virtual void moveForward(int scale = 1.0);
   virtual void moveBack(int scale = 1.0);

   virtual void turnLeft(int scale = 1.0);
   virtual void turnRight(int scale = 1.0);
   virtual void turnUp(int scale = 1.0);
   virtual void turnDown(int scale = 1.0);

   virtual void orbitLeft(int scale = 1.0);
   virtual void orbitRight(int scale = 1.0);
   virtual void orbitUp(int scale = 1.0);
   virtual void orbitDown(int scale = 1.0);

   // Reset to original state
   virtual void reset();

   // Conversion utilities between screen and world coordinates
   virtual bool screenToWorld(int screenX, int screenY, vec3& worldCoords);
   virtual bool worldToScreen(const vec3& worldCoords, int& screenX, int& screenY);

   // Get camera to world matrix
   virtual math::matrix<double> cameraToWorldMatrix();

protected:
   enum Dir { NONE, F, B, L, R, U, D, TL, TR, TU, TD} myDir, myTurnDir;
   virtual void turn(vec3& v, vec3& n, float amount);
   virtual void move(float dU, float dV, float dN);
   virtual void orbit(float h, float p);

protected:
   float mSpeed, mTurnRate;

   vec3 eye; // camera position
   float mHeading, mPitch, mRadius;
   float mVfov, mAspect, mNear, mFar; // projection parameters
   
   // Basis of camera local coord system
   vec3 u; // up
   vec3 v; // v points right
   vec3 n; // -n points forward

   // Cache useful values
   GLdouble myProjMatrix[16];
   GLdouble myModelMatrix[16];
   GLint myViewport[4];

public:

   // Defaults
   static vec3 dfltEye, dfltUp, dfltLook;
   static float dfltVfov, dfltAspect, dfltNear, dfltFar; 
   static float dfltSpeed, dfltTurnRate;
};

#endif
