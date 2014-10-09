#include "smoke_sim.h"
#include "camera.h"
#include "fps.h"
#include "constants.h"
#include <stdio.h>
#include <cmath>
#include "open_gl_headers.h"
#include "basic_math.h"

// Geometry and whatnot
SmokeSim theSmokeSim;
Camera theCamera;
mmc::FpsTracker theFpsTracker;

// UI Helpers
int lastX = 0, lastY = 0;
int theMenu = 0;
int theButtonState = 0;
int theModifierState = 0;
bool isRunning = true;

int savedWidth = 0;
int savedHeight = 0;

void initCamera()
{
   double w = theDim[0]*theCellSize;   
   double h = theDim[1]*theCellSize;   
   double d = theDim[2]*theCellSize;   
   double angle = 0.5*theCamera.dfltVfov*BasicMath::PI/180.0;
   double dist;
   if (w > h) dist = w*0.5/std::tan(angle);  // aspect is 1, so i can do this
   else dist = h*0.5/std::tan(angle);
   theCamera.dfltEye.set(w*0.5, h*0.5, -(dist+d*0.5));
   theCamera.dfltLook.set(w*0.5, h*0.5, 0.0);
   theCamera.reset();
}

void onMouseMotionCb(int x, int y)
{
   int deltaX = lastX - x;
   int deltaY = lastY - y;
   bool moveLeftRight = abs(deltaX) > abs(deltaY);
   bool moveUpDown = !moveLeftRight;

   if (theButtonState == GLUT_LEFT_BUTTON)  // Rotate
   {
      if (moveLeftRight && deltaX > 0) theCamera.orbitLeft(deltaX);
      else if (moveLeftRight && deltaX < 0) theCamera.orbitRight(-deltaX);
      else if (moveUpDown && deltaY > 0) theCamera.orbitUp(deltaY);
      else if (moveUpDown && deltaY < 0) theCamera.orbitDown(-deltaY);
   }
   else if (theButtonState == GLUT_MIDDLE_BUTTON) // Zoom
   {
      if (moveUpDown && deltaY > 0) theCamera.moveForward(deltaY);
      else if (moveUpDown && deltaY < 0) theCamera.moveBack(-deltaY);
   }    

   if (theModifierState & GLUT_ACTIVE_ALT) // camera move
   {
      if (theButtonState == GLUT_RIGHT_BUTTON) // Pan
      {
         if (moveLeftRight && deltaX > 0) theCamera.moveLeft(deltaX);
         else if (moveLeftRight && deltaX < 0) theCamera.moveRight(-deltaX);
         else if (moveUpDown && deltaY > 0) theCamera.moveUp(deltaY);
         else if (moveUpDown && deltaY < 0) theCamera.moveDown(-deltaY);
      }   
   }
 
   lastX = x;
   lastY = y;
   glutPostRedisplay();
}

void onMouseCb(int button, int state, int x, int y)
{
   theButtonState = button;
   theModifierState = glutGetModifiers();
   lastX = x;
   lastY = y;

   glutSetMenu(theMenu);
   if (theModifierState & GLUT_ACTIVE_ALT)
   {
      glutDetachMenu(GLUT_RIGHT_BUTTON);
   }
   else
   {
      glutAttachMenu(GLUT_RIGHT_BUTTON);
   }

   onMouseMotionCb(x, y);
}


void onKeyboardCb(unsigned char key, int x, int y)
{
   if (key == ' ') theCamera.reset();
   else if (key == '0') MACGrid::theRenderMode = MACGrid::CUBES;
   else if (key == '1') MACGrid::theRenderMode = MACGrid::SHEETS;
   else if (key == 'v') MACGrid::theDisplayVel = !MACGrid::theDisplayVel;
   else if (key == 'r') theSmokeSim.setRecording(!theSmokeSim.isRecording(), savedWidth, savedHeight);
   else if (key == '>') isRunning = true;
   else if (key == '=') isRunning = false;
   else if (key == '<') theSmokeSim.reset();
   else if (key == 27) exit(0); // ESC Key
   glutPostRedisplay();
}

void onMenuCb(int value)
{
   switch (value)
   {
   case -1: exit(0);
   case -6: theSmokeSim.reset(); break;
   default: onKeyboardCb(value, 0, 0); break;
   }
}

void onKeyboardSpecialCb(int key, int x, int y)
{
}

void onTimerCb(int value)
{
   if (isRunning) theSmokeSim.step();
   glutTimerFunc(theMillisecondsPerFrame, onTimerCb, 0);
   glutPostRedisplay();
}

void onResizeCb(int width, int height)
{
	// Save the width and height:
	savedWidth = width;
	savedHeight = height;
	
   // Update viewport
   glViewport(0, 0, width, height);

   // Update camera projection's aspect ratio
   float vfov, aspect, zNear, zFar;
   theCamera.getProjection(&vfov, &aspect, &zNear, &zFar);
   theCamera.setProjection(vfov, ((GLfloat) width)/height, zNear, zFar);
}

void drawOverlay()
{
  // Draw Overlay
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glPushAttrib(GL_LIGHTING_BIT);
     glDisable(GL_LIGHTING);

     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluOrtho2D(0.0, 1.0, 0.0, 1.0);

     glMatrixMode(GL_MODELVIEW);
     glLoadIdentity();
     glRasterPos2f(0.01, 0.01);
     
     char info[1024];
     sprintf(info, "Framerate: %3.1f  |  Frame: %u  |  %s", 
         theFpsTracker.fpsAverage(), theSmokeSim.getTotalFrames(),
         theSmokeSim.isRecording()? "Recording..." : "");
 
     for (unsigned int i = 0; i < strlen(info); i++)
     {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, info[i]);
     }
  glPopAttrib();
}

void onDrawCb()
{
	// Keep track of time
	theFpsTracker.timestamp();

	// Draw Scene and overlay
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	theCamera.draw();
	theSmokeSim.draw(theCamera);
	drawOverlay();
	glutSwapBuffers();
}

void init(void)
{
    initCamera();
    glClearColor(0.1, 0.1, 0.1, 1.0);

    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_NORMALIZE);
    glDisable(GL_LIGHTING);
    glCullFace(GL_BACK);
}


int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Fluid Simulation - CIS563");
    glutDisplayFunc(onDrawCb);
    glutKeyboardFunc(onKeyboardCb);
    glutSpecialFunc(onKeyboardSpecialCb);
    glutMouseFunc(onMouseCb);
    glutMotionFunc(onMouseMotionCb); 
    glutTimerFunc(theMillisecondsPerFrame, onTimerCb, 0);
    glutReshapeFunc(onResizeCb);

    int viewMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Toggle velocities\t'v'", 'v');
    glutAddMenuEntry("Render density as cubes\t'0'", '0');
    glutAddMenuEntry("Render density as sheets\t'1'", '1');

    theMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Start\t'>'", '>');
    glutAddMenuEntry("Pause\t'='", '=');
    glutAddMenuEntry("Reset\t'<'", '<');
    glutAddMenuEntry("Reset camera\t' '", ' ');
    glutAddMenuEntry("Record\t'r'", 'r');
    glutAddSubMenu("Display", viewMenu);
    glutAddMenuEntry("_________________", -1);
    glutAddMenuEntry("Exit", 27);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    init();

    glutMainLoop();
    return 0;             
}

