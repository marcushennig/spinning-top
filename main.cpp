//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Date:		 10/21/2004  
// Name:          Marcus Hennig, University of Chapel Hill
// File:		 main.cpp
// Description:  symmetric spinning top, project phys61
//			
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// INCLUDE FILES
// the usual stuff
#include <iostream>
#include <fstream>
#include <string.h>
//#include <stdio.h>
#include <math.h>
//#include <stdlib.h>
#include <memory.h>
#include <time.h>

// opengl to visualize the simulation
//#include <windows.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

// my own stuff
#include "Geometry.h"
#include "vector4d.h"
#include "trackball.h"

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// NAMESPACES
// I will never use, but a old bad habit
using namespace std;

double norm(double x)
{
	if (x < 0)
		return -x;
	else
		return x;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// GLOBAL VARIABLES AND STRUCTS
// althought I hate global parameters, it's necessary to increase speed  

// the parameter structur contains all current viariables, especially the angles values
// are very IMPORTANT for the simulation, a lot of stuff in this struct is time independent
struct parameter
{
	 double theta;
	 double phi;
	 double psi;
	 double p_theta;
	 double p_phi;
	 double p_psi;
	 double Dtheta;
	 double Dphi;
	 double Dpsi;
	 double M;
	 double g;
	 double l;
	 double J1;
	 double J3;
	 double E;
	 int       itermax;
	 double accuracy;
	 double dt;
};

// little buffer to store my datas for the graphs
class dBuffer
{
private:
	 double* data;
	 int size;
	 int length;
public:
	 dBuffer(int n) 
	 { 
		  size = n; 
		  length = 0;
		  data = new double[n]; 
	 }
	 
	 ~dBuffer() 
	 { 
		  if(data) delete [] data; 
	 }
	 
	 double operator[](int j)
	 {
		  return data[j];
	 }
	 
	 int GetSize(void)
	 {
		  return size;
	 }
	
	 int GetLength(void)
	 {
		  return length;
	 }
	 
	 void Clear(void)
	 {
		  length = 0;
	 } 
	 
	 double* GetPointer(void)
	 {
		  return data;
	 }
	 
	 void Push(double x, double y)
	 {
		  if(length==size)
		  {
				memmove(data, data+2, (length-2)*sizeof(double));
				data[length-2] = x;
				data[length-1] = y;
		  }
		  if(length<size)
		  {
				data[length++] = x;
				data[length++] = y;
		  }
	 }
	 
	 void GetRangeY(double& min, double& max)
	 {
		  if(length>0)
		  {
				min = data[1]; max = min;
				double d;
				int j = 1;
				while(j<length)
				{
					 d= data[j];
					 if(d>max)
						  max = d;
					 else if(d<min)
						  min = d;
					 j += 2;
				}
		  }
	 }
	 
	 void GetRangeX(double& min, double& max)
	 {
		  if(length>0)
		  {
				min = data[0]; max = min;
				double d;
				int j = 0;
				while(j<length)
				{
					 d= data[j];
					 if(d>max)
						  max = d;
					 else if(d<min)
						  min = d;
					 j += 2;
				}
		  }
	 }
};

// global visual parameters
typedef struct 
{
	 GLdouble x,y,z;
} recVec;

typedef struct 
{
	 recVec viewPos;		// View position
	 recVec viewDir;		// View direction vector
	 recVec viewUp;		// View up direction
	 recVec rotPoint;		// Point to rotate about
	 GLdouble focalLength;	// Focal Length along view direction
	 GLdouble aperture;	// gCamera aperture
	 GLdouble eyeSep;		// Eye separation
	 GLint screenWidth,screenHeight; // current window/screen height and width
} recCamera;

// sorry, but this small letter have the main role in this play (see declaration of parameter) 
// and it's my favorite letter since ceasars name begins with it :-)
parameter c;

//global time and timestep
double t = 0.0;

// some paramters for opengl and the glut API
// gPoint... etc store a handle for a small ASM-program 
int gMainWindow = 0;
GLuint gTopList   = NULL;
GLuint gPointList = NULL;
GLuint gWireList = NULL;
GLuint gSolidList = NULL;

// global dimension of my drawing window
int width;
int height;

int currentMenuItem = 0;

// some switchers
bool bThetaPotential = TRUE;
bool bPhiTheta	   = TRUE;
bool bDPhiTheta	   = FALSE;
bool bPeriodicKicks    = FALSE;
bool bPhiMotion	   = TRUE;
bool bThetaMotion    = TRUE;
bool bPsiMotion	   = TRUE;
bool bSolid		   = TRUE;
bool bWire		   = FALSE;
bool bPoints		   = FALSE; 
bool bTop		   = FALSE;
bool bAnimation	   = TRUE;
bool bWriteTextFile    = FALSE;

// initialize global buffers
dBuffer bufferThetaPotential(100);
dBuffer bufferPhiTheta(400);
dBuffer bufferDPhiTheta(100);

// file object to record the data
ofstream outFile;

// rotation parameters
const double DTOR = 0.0174532925;

GLboolean gTrackBall = GL_FALSE;

recCamera gCamera;
recVec gOrigin				 = { 0.0, 0.0, 0.0 };
GLfloat gShapeSize			 = 11.0f;
GLfloat gWorldRotation [4]	 = { 155.0, 0.0, -1.0, 0.0 };
GLfloat gTrackBallRotation [4] = { 0.0, 0.0, 0.0, 0.0 };

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// actualize the energy within the global structure c, I know it's better to use c class to keep classiness
//but in order to increase speed I'mforced to introduce some rough ugly styles :(
void UpdateEnergy(void)
{
	c.E = c.p_theta*c.p_theta / (2 * c.J1) + c.M*c.g*c.l*cos(c.theta) + c.p_psi*c.p_psi / (2*c.J3) + pow( (c.p_phi - c.p_psi*cos(c.theta)), 2.0) / (2*c.J1*pow(sin(c.theta), 2.0)); 
}

//update the angelar momentum in case that only the time derivatives of the 
// euler angles are known
void UpdateAngularMomentum(void)
{
	 c.p_theta = c.J1 * c.Dtheta;
	 c.p_phi    = c.J3 * cos(c.theta) * c.Dpsi + (c.J1 * pow(sin(c.theta), 2) + c.J3 * pow(cos(c.theta), 2)) * c.Dphi;
	 c.p_psi    = c.J3 * (cos(c.theta) * c.Dphi + c.Dpsi);
}

void UpdateAngularVelocity(void)
{
	 c.Dtheta = c.p_theta / c.J1;
	 c.Dphi    = (c.p_phi - c.p_psi * cos(c.theta)) / (pow(sin(c.theta), 2) * c.J1);
	 c.Dpsi    = (c.J1 * c.p_psi + c.J3 * c.p_psi * pow(cos(c.theta), 2) - c.J3 * c.p_phi * cos(c.theta)) / (c.J1 * c.J3 * pow(sin(c.theta), 2));
}

// see in the good old GOLDSTEIN, topic THE HAVY SYMMETRICAL TOP WITH ONE FIXED POINT and study 
// the chapter, this function describes the time derivative of phi wihich depends on theta
double DPhi(double theta)
{
	 return (c.p_phi - c.p_psi*cos(theta))/(c.J1*pow(sin(theta), 2));
}

// this small potential function reveals a first impression how the motion could be, also check up
// GOLDSTEIN
double ThetaPotential(double theta)
{ 
	 if(norm(theta)<1E-10 && norm(c.p_phi - c.p_psi)<1E-10)
		  return  c.M*c.g*c.l;
	 else
		  return c.M*c.g*c.l*cos(theta) + c.J1/2.0 * pow((c.p_phi - c.p_psi*cos(theta))/(c.J1*sin(theta)),2);
	 
}

// this is grad H multiplied with the symplectric matrix, Dt(x) = S*grad(H), motion equation provieded by the 
// Hamilton formalism.... or the right hand side of our 1st order differential equation system
vector4d F(vector4d x)
{
	 double theta     = x[1];
	 double phi        = x[2];
	 double psi         = x[3];
	 double p_theta = x[4];
		  
	 double s_theta = sin( theta );
	 double c_theta = cos( theta );
	 double A          = c.p_phi - c.p_psi * c_theta;
	 
	 return vector4d( p_theta / c.J1,
				   A / (c.J1 * s_theta*s_theta),
				   c.p_psi / c.J3 - c_theta *A / (c.J1 * s_theta*s_theta),
				   A * ( -c.p_psi / (c.J1 * s_theta) + A * c_theta / (c.J1 * s_theta*s_theta*s_theta) )  + c.M*c.g*c.l*s_theta );
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// step for the Runge-Kutta method 
vector4d step( vector4d x, double h)
{
	 vector4d k1 = h * F( x );
	 vector4d k2 = h * F( x + 0.5 * k1 );
	 vector4d k3 = h * F( x + 0.5 * k2 );
	 vector4d k4 = h * F( x + 1.0 * k3 );
	 vector4d t =  ( k1 + 2*k2 + 2*k3 + k4 ) / 6.0;
	 return t;
}

// Runge - Kutta method with step size controlling
void NDSolveSpinningTop( double a, double b )
{
	 vector4d x( c.theta, c.phi, c.psi, c.p_theta ); 
	 
	 double t  = a;
	 
	 double h = (b-a) / 10.0; 
	 double hmax = (b-a) / 2.0;
	 
	 double error = 0;
	 vector4d u, v, d;
	 
	 double rho = 0.9;
	 int q = 2;
	 int p = 4;
	 double H;			
	 double eps = 1E-12;
	 int iter = 0;
	 
	 while( (t < b) && (iter < c.itermax))
	 {
		  iter++;
		  h = min( h, b - t ); 
		  H = h; 
		  u = x + step(x, H/2);
		  u = u + step(u, H/2);
		  v = x + step(x,  H );
		  d = u - v;
		  error = sqrt( d*d )/15.0;
		  
		  if (error < c.accuracy)
		  {
				t = t + H;
				x = u;
		  }					 
		  if(error >= eps)
				h = min(q*H, min(rho*H*pow(c.accuracy/error, 0.2), hmax) ); 
		  else
				h = min(q*H, hmax);
				
		  if (h < eps)
				h = 2*eps;
	 }; 
	 
	 c.theta     = x[1];
	 c.phi        = x[2];
	 c.psi         = x[3];
	 c.p_theta = x[4];
}

// in case the user is not sensible with the intial conditions, he or she can reset 
// all intial parameter
void ResetInitialConditions(void)
{
	 c.J1			 = 0.75;
	 c.J3			 = 5.10;
	 c.M			 = 1.65;
	 c.g			 = 9.81;
	 c.l			 = 80;
	 c.accuracy	 = 1E-4;
	 c.itermax	 = 10000;
	 c.theta		 = 0.5*M_PI;
	 c.phi		 = 0;
	 c.psi		 = 0;
	 c.Dtheta		 = 0.0;
	 c.Dphi		 = 0.0;
	 c.Dpsi		 = 56;
	 c.dt			 = 0.001;
	 UpdateAngularMomentum();
	 UpdateEnergy();
	 
	 t   = 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// the same as before
void gCameraReset(void)
{
	 gCamera.aperture = 40;
	 gCamera.focalLength = 15;
	 gCamera.rotPoint = gOrigin;
	 
	 gCamera.viewPos.x = 0.0;
	 gCamera.viewPos.y = 0.0;
	 gCamera.viewPos.z = -gCamera.focalLength;
	 gCamera.viewDir.x = -gCamera.viewPos.x; 
	 gCamera.viewDir.y = -gCamera.viewPos.y; 
	 gCamera.viewDir.z = -gCamera.viewPos.z;
	 
	 gCamera.viewUp.x = 0;  
	 gCamera.viewUp.y = 1; 
	 gCamera.viewUp.z = 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CREATING OF THE GEOMETRY DATA
//
float f(float x)
{
	 return sinf(x);
}

float Df(float x)
{
	 return cosf(x);
}

void Surface( float x[], float n[], float u, float v)
{
	 x[0] = cosf(u) * f(v);
	 x[1] = v;
	 x[2] = sinf(u)  * f(v);
	 n[0] = -cosf(u) * f(v);
	 n[1] = f(v) * Df(v);
	 n[2] =-sinf(u) * f(v);
}

void CreateTop(void)
{
	gTopList = glGenLists(1);
	 glNewList(gTopList, GL_COMPILE); 
	 {
		  float v[3];
		  float n[3];
		  float r = 0;
		  float s = 0;
		  float dr  = 2 * M_PI / 50;
		  float ds = M_PI / 50;
	 
	 
		  glBegin(GL_QUADS);
		  while(s < M_PI)
		  {
				r = 0;
				while(r < 2*M_PI)
				{
					 Surface( v, n, r, s );
					 glNormal3fv( n );
					 glVertex3fv ( v );
				
					 Surface( v, n, r+dr, s );
					 glNormal3fv( n );
					 glVertex3fv ( v );
				
					 Surface( v, n, r+dr, s+ds );
					 glNormal3fv( n );
					 glVertex3fv ( v );
				
					 Surface( v, n, r, s+ds );
					 glNormal3fv( n );
					 glVertex3fv ( v );
				
					 r += dr;
				};
				s += ds;
		  };
		  glEnd();
	 }
	 glEndList();
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// LIGHTING

void SetLighting(unsigned int mode)
{
	 GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
	 GLfloat mat_shininess[] = {90.0};
	 
	 GLfloat position[4] = {7.0, 0.0, 12.0,0.0};
	 GLfloat ambient[4]  = {0.2,0.2,0.2,1.0};
	 GLfloat diffuse[4]  = {1.0,1.0,1.0,1.0};
	 GLfloat specular[4] = {1.0,1.0,1.0,1.0};
	 
	 glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	 glMaterialfv (GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	 
	 glEnable(GL_COLOR_MATERIAL);
	 glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
	 
	 switch (mode) {
		  case 0:
				break;
		  case 1:
				glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
				glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_FALSE);
				break;
		  case 2:
				glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
				glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
				break;
		  case 3:
				glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
				glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_FALSE);
				break;
		  case 4:
				glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
				glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
				break;
	 }
	 
	 glLightfv(GL_LIGHT0,GL_POSITION,position);
	 glLightfv(GL_LIGHT0,GL_AMBIENT,ambient);
	 glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);
	 glLightfv(GL_LIGHT0,GL_SPECULAR,specular);
	glEnable(GL_LIGHT0);
}

void InitGL(void)
{
	 glShadeModel(GL_SMOOTH);
	 glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	 SetLighting(4);
	 gCameraReset ();
	 CreateTop();
	 BuildGeometry(kTranguloidTrefoil, 4, 50, 3, &gSolidList, &gWireList, &gPointList);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------
// draw some strings on the screen
void DrawGLString(GLfloat x, GLfloat y, char *string, int size = 10)
{
	 int length;
	 
	 glRasterPos2f(x, y);
	 length = (int) strlen(string);
	 
	 for(int i = 0; i < length; i++) 
	 {
		  switch( size )
		  {
				case 10:
					 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, string[i]);
					 break;
				case 12:
					 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, string[i]); 
					 break;
				case 18:
					 glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]); 				  
		  };
	 }
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void DrawGLText (void)
{
	 char outString [256] = "";
	 GLint lineSpacing = 20;
	 GLint line = 0;
	 GLint startOffest = -7;
	 	 
	 //mark current menu
	 double y[12];
	 line = 3;
	 for(int i=0; i<12; i++)
	 {
		  y[i] = lineSpacing * line++ ;
	 }
	 
	 glBegin(GL_QUADS);
		  glColor3f( 1.0f, 0.4f, 1.0f );	
		  glVertex2d( 0, y[currentMenuItem] );
		  glVertex2d( 0, y[currentMenuItem] - lineSpacing );
		  glColor3f( 0.3f, 0.3f, 0.4f );
		  glVertex2d( 200, y[currentMenuItem] - lineSpacing );
		  glVertex2d( 200, y[currentMenuItem] );
	 glEnd();
	 
	 glColor3f (1.0, 1.0, 1.0);
	 
	 sprintf (outString, "Inital Conditions", c.M);
	 DrawGLString (10, (lineSpacing * 2) + startOffest, outString,18);
	 
	 sprintf (outString, "mass m   = %5.2f kg", c.M);
	 DrawGLString (10, y[0]+ startOffest, outString);
	 sprintf (outString, "gravity g = %5.2f m/s^2", c.g);
	 DrawGLString (10, y[1]+ startOffest, outString);
	 sprintf (outString, "length l   = %5.2f m", c.l);
	 DrawGLString (10, y[2]+ startOffest, outString);
	 sprintf (outString, "moment of intertia J1  = %5.2f kg*m^2", c.J1);
	 DrawGLString (10, y[3]+ startOffest, outString);
	 sprintf (outString, "moment of intertia J3  = %5.2f kg*m^2", c.J3);
	 DrawGLString (10, y[4]+ startOffest, outString);
	 sprintf (outString, "theta  = %5.2f Pi", c.theta / M_PI);
	 DrawGLString (10, y[5]+ startOffest, outString);
	 sprintf (outString, "phi  = %5.2f Pi", c.phi / M_PI);
	 DrawGLString (10, y[6]+ startOffest, outString);
	 sprintf (outString, "psi  = %5.2f Pi", c.psi / M_PI);
	 DrawGLString (10, y[7]+ startOffest, outString);
	 sprintf (outString, "Dtheta  = %5.2f Pi (p_theta  = %5.2f)", c.Dtheta / M_PI, c.p_theta );
	 DrawGLString (10, y[8]+ startOffest, outString);
	 sprintf (outString, "Dphi  = %5.2f Pi (p_phi  = %5.2f)", c.Dphi / M_PI, c.p_phi  );
	 DrawGLString (10, y[9]+ startOffest, outString);
	 sprintf (outString, "Dpsi  = %5.2f Pi (p_psi  = %5.2f)", c.Dpsi / M_PI, c.p_psi  );	 
	 DrawGLString (10, y[10]+ startOffest, outString);
	 sprintf (outString, "dt  = %2.5f s", c.dt );	 
	 DrawGLString (11, y[11]+ startOffest, outString);
 
	 line = 15;
	 glColor3f(1, 0.6, 0);
	 sprintf (outString, "time t  = %8.3f s", t);
	 DrawGLString (10, lineSpacing * (line++)+ startOffest, outString);
	 sprintf (outString, "Energy E = %5.4f J", c.E);
	 DrawGLString (10, lineSpacing * (line++)+ startOffest, outString);
	 
	 sprintf (outString, "motion theta: %d phi: %d psi: %d", bThetaMotion, bPhiMotion, bPsiMotion );
	 DrawGLString (10, lineSpacing * (line++)+ startOffest, outString);
	 
	 if(bWriteTextFile) DrawGLString (10, lineSpacing * (line++)+ startOffest, "store data");
	 
	 glColor3f(1, 1, 1);
	 if(bThetaPotential) DrawGLString(10, 2*height/3, "theta potential", 12);
	 if(bPhiTheta) DrawGLString(10+width/4, 2*height/3, "parameter function (theta(t), phi(t))", 12);
	 if(bDPhiTheta) DrawGLString(10+5*width/8, 2*height/3, "parameter function (theta(t), Dphi(t))", 12);
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// AXES
void DrawGLAxes(float sx=1, float sy=1, float sz=1 )
{
	 glLineWidth(2.0f);
	 glBegin(GL_LINES);
		  glColor3f( 0.0f, 0.0f, 1.0f );	  
		  glVertex3f( 0, 0, 0 );
		  glVertex3f( sx, 0, 0 );
		  glColor3f( 0.0f, 1.0f, 0.0f );	
		  glVertex3f( 0, 0, 0 );
		  glVertex3f( 0, sy, 0 );
		  glColor3f( 1.0f, 0.0f, 0.0f );	
		  glVertex3f( 0, 0, 0 );
		  glVertex3f( 0, 0, sz );
	 glEnd();
	 glPushMatrix();
		  glColor3f( 1.0f, 0.0f, 0.0f );	  
		  glTranslatef( 0, 0, sz );
		  glutSolidCone( 0.2, 0.4, 10, 10 );
	 glPopMatrix();
	 glPushMatrix();
		  glColor3f( 0.0f, 1.0f, 0.0f );	  
		  glTranslatef( 0, sy, 0 );
		  glRotatef( -90, 1, 0, 0 );
		  glutSolidCone( 0.2, 0.4, 10, 10 );
	 glPopMatrix();
	 glPushMatrix();
		  glColor3f( 0.0f, 0.0f, 1.0f );	  
		  glTranslatef( sx, 0., 0 );
		  glRotatef( 90, 0, 1, 0 );
		  glutSolidCone( 0.2, 0.4, 10, 10 );
	 glPopMatrix();

}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// simulates instantaneous kicks, all 20 time steps a constant value is add to g0
double Gravity(double t)
{
	 const double g0 = 9.81;
	 double T   = 20*c.dt;
	 double n = t / T;
	 double d = norm(n - round(n));
	 if(d < 1E-8)
		  return g0+100;
	 else 
		  return g0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// draw a very pretty body, colorful, sophisticated twisted... this was the idea of my girlfriend
void DrawGLBody(void)
{	 
	 DrawGLAxes( 4, 4, 4 );
	 
	 //evaluate small timestep, function acutalize the parameterlist
	 if( bAnimation )
	 {
		  NDSolveSpinningTop(t, t + c.dt);
		  t += c.dt; 
		  UpdateAngularVelocity();
		  
		  // record results 
		  if( bWriteTextFile )
		  {
				outFile<<t<<"\t"<<c.theta<<"\t"<<c.phi<<"\t"<<c.psi<<"\n";
		  }
	 }
	 
	 if( bAnimation && bPhiTheta )  bufferPhiTheta.Push(c.phi, c.theta);
	 if( bAnimation && bDPhiTheta) bufferDPhiTheta.Push(c.theta, DPhi(c.theta));
	 if( bAnimation && bPeriodicKicks ) c.g = Gravity( t );
	 // eulerrotation
	 if( bPhiMotion )     glRotated( 180 * c.phi / M_PI,    0, 1, 0 );
	 if( bThetaMotion ) glRotated( 180 * c.theta / M_PI, 1, 0, 0 );
	 if( bPsiMotion )      glRotated( 180 * c.psi / M_PI,     0, 1, 0 );

	 DrawGLAxes(1, 4, 1 );
	 glColor3f(1, 0.5,0);
	 
	 if(!bTop)
	 {
		  if(bSolid) 
				glCallList (gSolidList);
		  else if(bWire)
				glCallList (gWireList);
		  else if(bPoints)
				glCallList (gPointList);
	 }
	 else
	 {
		  glCallList (gTopList);
	 }
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// colorize the screen 
void DrawGLScene(void)
{
	 glClearColor (0.3f, 0.3f, 0.4f, 1.0f);
	 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	 
	 width  = glutGet(GLUT_WINDOW_WIDTH);
	 height = glutGet(GLUT_WINDOW_HEIGHT);
	 	 
	 for(int screen=0; screen<5; screen++)
	 {
		  if(screen==0)
		  {
				int w = width;
				int h = height;
				
				gCamera.screenWidth = w;
				gCamera.screenHeight = h;
				
				glViewport(0, 0, w, h);
				
				GLdouble xmin, xmax, ymin, ymax;
				// far frustum plane
				GLdouble zFar = -gCamera.viewPos.z + gShapeSize * 0.5;
				// near frustum plane clamped at 1.0
				GLdouble zNear = min(-gCamera.viewPos.z - gShapeSize * 0.5, 1.0);
				// window aspect ratio
				GLdouble aspect = gCamera.screenWidth / (GLdouble)gCamera.screenHeight; 
				
				glMatrixMode (GL_PROJECTION);
				glLoadIdentity ();
				
				if (aspect > 1.0) {
					 ymax = zNear * tan (gCamera.aperture * 0.5 * DTOR);
					 ymin = -ymax;
					 xmin = ymin * aspect;
					 xmax = ymax * aspect;
				} else {
					 xmax = zNear * tan (gCamera.aperture * 0.5 * DTOR);
					 xmin = -xmax;
					 ymin = xmin / aspect;
					 ymax = xmax / aspect;
				}
				glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
				
				glMatrixMode (GL_MODELVIEW);
				glLoadIdentity ();
				gluLookAt (gCamera.viewPos.x, gCamera.viewPos.y, gCamera.viewPos.z,
							  gCamera.viewPos.x + gCamera.viewDir.x,
							  gCamera.viewPos.y + gCamera.viewDir.y,
							  gCamera.viewPos.z + gCamera.viewDir.z,
							  gCamera.viewUp.x, gCamera.viewUp.y ,gCamera.viewUp.z);
				
				// track ball rotation
				glRotatef (gTrackBallRotation[0], gTrackBallRotation[1], gTrackBallRotation[2], gTrackBallRotation[3]);
				glRotatef (gWorldRotation[0], gWorldRotation[1], gWorldRotation[2], gWorldRotation[3]);
				
				glEnable(GL_LIGHTING);
				glEnable(GL_DEPTH_TEST);
		
				DrawGLBody();
		  }
		  if(screen==1)
		  {
				glViewport(0, 0, width, height);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				gluOrtho2D( 0, width, height, 0);
				glDisable(GL_LIGHTING);
				glDisable(GL_DEPTH_TEST);
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				DrawGLText();
		  }
		  if(screen==2 && bThetaPotential)
		  {
				int w = width/4;
				int h = height/3;

				glViewport(0, 0, w, h);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				
				// create vertexdata for the graph
				int N = bufferThetaPotential.GetSize()/2;
			
				double theta = -3.14;
				double dtheta = 2*3.14 / N;
				double min = ThetaPotential(theta);
				double V;
				double *data = bufferThetaPotential.GetPointer();
				int j=0;
				while(j<2*N)				
				 {
					 V = ThetaPotential(theta);
					 if(V>c.E) V = c.E;   
					 if(V<min) min = V;
					 data[j++] = theta;
					 data[j++] = V;
					 theta += dtheta;
				}
				
				double xmin = -M_PI;
				double xmax = M_PI;
				double ymin = min;
				double ymax = c.E;
				double dx  = 0.05*(xmax - xmin);
				double dy  = 0.05*(ymax - ymin);
				
				gluOrtho2D( xmin-dx, xmax+dx, ymin-dy, ymax+dy );
				
				glDisable(GL_LIGHTING);
				glDisable(GL_DEPTH_TEST);
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glColor3f( 1, 1, 1 );
				
				glBegin(GL_LINE_STRIP);
					 glVertex2d(xmin, ymin);
					 glVertex2d(xmin, ymax);
					 glVertex2d(xmax, ymax);
					 glVertex2d(xmax, ymin);
					 glVertex2d(xmin, ymin);
				glEnd();
				
				glEnableClientState ( GL_VERTEX_ARRAY ) ;
					 glVertexPointer(2, GL_DOUBLE, 0, data);
					 glDrawArrays(GL_LINE_STRIP, 0, N);
				glDisableClientState(GL_VERTEX_ARRAY);
				
				glPointSize(8);
				glColor3f(1,0,0);
				glBegin(GL_POINTS);
					 glVertex2d(c.theta,ThetaPotential(c.theta));
				glEnd();
		  }
		  if(screen==3 && bPhiTheta)
		  {
				int w;
				
				if( bDPhiTheta == FALSE )
					 w = 3*width/4;
				else 
					 w = 3*width/8;
				
				int h = height/3;
				
				glViewport(width/4, 0, w, h);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();

				double xmin, xmax, ymin, ymax;
				bufferPhiTheta.GetRangeX( xmin, xmax);
				bufferPhiTheta.GetRangeY( ymin, ymax);
				
				double dx  = 0.05*(xmax - xmin);
				double dy  = 0.05*(ymax - ymin);
				
				gluOrtho2D( xmin-dx, xmax+dx, ymin-dy, ymax+dy );
				
				glDisable(GL_LIGHTING);
				glDisable(GL_DEPTH_TEST);
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
		
				glColor3f( 1, 1, 1 );
				
				glBegin(GL_LINE_STRIP);
					 glVertex2d(xmin, ymin);
					 glVertex2d(xmin, ymax);
					 glVertex2d(xmax, ymax);
					 glVertex2d(xmax, ymin);
					 glVertex2d(xmin, ymin);
				glEnd();
				glBegin(GL_LINE_STRIP);
					 
				glEnd();
				glEnableClientState ( GL_VERTEX_ARRAY ) ;
				glVertexPointer(2, GL_DOUBLE, 0, bufferPhiTheta.GetPointer());
				glDrawArrays(GL_LINE_STRIP, 0, bufferPhiTheta.GetLength()/2);
				glDisableClientState(GL_VERTEX_ARRAY);
		  }
		  if(screen==4 && bDPhiTheta)
		  {
				int w = 3*width/8;
				int h = height/3;
				
				glViewport((5*width)/8, 0, w, h);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				
				// create vertexdata for the graph
				double xmin, xmax, ymin, ymax;
				bufferDPhiTheta.GetRangeX( xmin, xmax);
				bufferDPhiTheta.GetRangeY( ymin, ymax);
				
				double dx  = 0.05*(xmax - xmin);
				double dy  = 0.05*(ymax - ymin);
				
				gluOrtho2D( xmin-dx, xmax+dx, ymin-dy, ymax+dy );
				
				glDisable(GL_LIGHTING);
				glDisable(GL_DEPTH_TEST);
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glColor3f( 1, 1, 1 );
				
				glBegin(GL_LINE_STRIP);
					 glVertex2d(xmin, ymin);
					 glVertex2d(xmin, ymax);
					 glVertex2d(xmax, ymax);
					 glVertex2d(xmax, ymin);
					 glVertex2d(xmin, ymin);
				glEnd();
				if( ymin<=0 && 0<ymax )
					 
				glColor3f( 1, 0.7, 0 );
				glBegin(GL_LINE_STRIP);
					 glVertex2d(xmin, 0);
					 glVertex2d(xmax, 0);
				glEnd();
				
				glColor3f( 1, 1, 1 );
				glEnableClientState ( GL_VERTEX_ARRAY ) ;
					 glVertexPointer(2, GL_DOUBLE, 0, bufferDPhiTheta.GetPointer());
					 glDrawArrays(GL_LINE_STRIP, 0, bufferDPhiTheta.GetLength()/2);
				glDisableClientState(GL_VERTEX_ARRAY);
				
				glPointSize(8);
				glColor3f(1,0,0);
				glBegin(GL_POINTS);
					 glVertex2d(c.theta,DPhi(c.theta));
				glEnd();
				
		  }
	 }
	
	 glFlush();
	 glutSwapBuffers();
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void Keyboard ( unsigned char key, int x, int y )  // Create Keyboard Function
{
	 switch ( key )
	 {
		  case '1':
				bThetaPotential = !bThetaPotential;
				break;
		  case '2':
				bPhiTheta = !bPhiTheta;
				break;
		  case '3':
				bDPhiTheta = !bDPhiTheta;
				break;
		  case 'q':
				bThetaMotion =  !bThetaMotion;
				break;
		  case 'w':
				bPhiMotion = !bPhiMotion;
				break;
		  case 'e':
				bPsiMotion = !bPsiMotion; 
				break;
		  case 'g':
				bPeriodicKicks =!bPeriodicKicks;
				break;
		  case 'r':
				ResetInitialConditions();
				break;
		  case 't':
				bTop = !bTop; 
				break;
		  case 'y':
				if(bSolid)
				{
					 bSolid  = FALSE;
					 bWire   = TRUE;
					 bPoints = FALSE;
				}
				else if(bWire)
				{
					 bSolid  = FALSE;
					 bWire   = FALSE;
					 bPoints = TRUE;
				}
				else if(bPoints)
				{
					 bSolid  = TRUE;
					 bWire   = FALSE;
					 bPoints = FALSE;
				}

				break;
		  case 'u':
				t = 0;
				bAnimation = !bAnimation;
				break;
		  case 'i':
				if(bWriteTextFile)
				{
					 bWriteTextFile = FALSE;
					 outFile.close();
				}
				else
				{
					 bWriteTextFile = TRUE;
					 /*time_t rawtime;
					 struct tm * timeinfo;
					 time ( &rawtime );
					 timeinfo = localtime ( &rawtime );*/
					 outFile.open("solution.txt", ofstream::out | ofstream::app);
				}
				break;
		  case 27:        
				exit ( 0 );   
				break;       
		  default:        				
				break;
	 }
}

void SpecialKeys ( int a_keys, int x, int y )  // Create Special Function (required for arrow keys)
{
	 switch ( a_keys )
	 {
		  case GLUT_KEY_UP:
				currentMenuItem--;
				if(currentMenuItem < 0) currentMenuItem=11;
				break;
		  case GLUT_KEY_DOWN: 
				currentMenuItem++;
				if(currentMenuItem > 11) currentMenuItem=0;
				break;
		  case GLUT_KEY_LEFT:
				switch(currentMenuItem)
				{
					 case 0:
						  c.M -= 0.01;
						  if(c.M<0) c.M=0.0;
						  UpdateEnergy();
						  break;
					 case 1: 
						  c.g -= 0.1;
						  UpdateEnergy();						  
						  break;
					 case 2: 
						  c.l -= 0.1;
						  if(c.l<0) c.l=0.0;
						  UpdateEnergy();
						  break;
					 case 3: 
						  c.J1 -= 0.01;
						  if(c.J1<0) c.J1=0.0;
						  UpdateEnergy();
						  break;
					 case 4: 
						  c.J3 -= 0.01;
						  if(c.J3<0) c.J3=0.0;
						  UpdateEnergy();
						  break;
					 case 5: 
						  c.theta -= 0.01 * M_PI;
						  UpdateEnergy();
						  break;
					 case 6: 
						  c.phi -= 0.01 * M_PI;
						  UpdateEnergy();
						  break;
					 case 7: 
						  c.psi -= 0.01 * M_PI;
						  UpdateEnergy();
						  break;
					 case 8: 
						  c.Dtheta -= 0.1 * M_PI;
						  UpdateEnergy();
						  UpdateAngularMomentum();
						  break;
					 case 9: 
						  c.Dphi -= 0.1 * M_PI;
						  UpdateEnergy();
						  UpdateAngularMomentum();
						  break;
					 case 10: 
						  c.Dpsi -= 0.1 * M_PI;
						  UpdateEnergy();
						  UpdateAngularMomentum();
						  break;
					 case 11:
						  c.dt += 0.0001;
						  break;
				};
				break;
		  case GLUT_KEY_RIGHT:
				switch(currentMenuItem)
				{
					 case 0:
						  c.M += 0.01;
						  UpdateEnergy();
						  break;
					 case 1: 
						  c.g += 0.1;
						  UpdateEnergy();
						  break;
					 case 2: 
						  c.l += 0.1;
						  UpdateEnergy();
						  break;
					 case 3: 
						  c.J1 += 0.01;
						  UpdateEnergy();
						  break;
					 case 4: 
						  c.J3 += 0.01;
						  UpdateEnergy();
						  break;
					 case 5: 
						  c.theta += 0.01 * M_PI;
						  UpdateEnergy();
						  break;
					 case 6: 
						  c.phi += 0.01 * M_PI;
						  UpdateEnergy();
						  break;
					 case 7: 
						  c.psi += 0.01 * M_PI;
						  UpdateEnergy();;
						  break;
					 case 8: 
						  c.Dtheta += 0.1 * M_PI;
						  UpdateAngularMomentum();
						  UpdateEnergy();
						  break;
					 case 9: 
						  c.Dphi += 0.1 * M_PI;
						  UpdateAngularMomentum();
						  UpdateEnergy();
						  break;
					 case 10: 
						  c.Dpsi += 0.1 * M_PI;
						  UpdateAngularMomentum();
						  UpdateEnergy();						  
						  break;				
					 case 11:
						  c.dt -= 0.0001;
						  break;

				};
				break;
		  case GLUT_KEY_F1:     // When Up Arrow Is Pressed...
				glutFullScreen ( ); // Go Into Full Screen Mode
				break;
		  case GLUT_KEY_F2:               // When Down Arrow Is Pressed...
				glutReshapeWindow ( 1024, 768 ); 
				break;
		  default:
				break;
	 }
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void MouseTrackball (int x, int y)
{
	 if (gTrackBall) {
		  rollToTrackball (x, y, gTrackBallRotation);
		  glutPostRedisplay();
	 }
}

void Mouse (int button, int state, int x, int y)
{
	 if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN))
	 {
		  startTrackball (x, y, 0, 0, gCamera.screenWidth, gCamera.screenHeight);
		  glutMotionFunc (MouseTrackball);
		  gTrackBall = GL_TRUE;
	 }
	 else if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_UP)) 
	 {
		  gTrackBall = GL_FALSE;
		  glutMotionFunc (NULL);
		  rollToTrackball (x, y, gTrackBallRotation);
		  if (gTrackBallRotation[0] != 0.0)
				addToRotationTrackball (gTrackBallRotation, gWorldRotation);
		  gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
	 }
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// main function
int main (int argc, const char * argv[])
{	 
	 ResetInitialConditions();
	 
	 glutInit(&argc, (char**)argv);
	 glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // non-stereo for main window
	 glutInitWindowPosition (0, 0);
	 glutInitWindowSize (1024, 768);
	 gMainWindow = glutCreateWindow("Spinning Top");
	 
	 InitGL(); // standard GL init
	 
	 glutDisplayFunc  (DrawGLScene);   
	 glutKeyboardFunc (Keyboard);
	 glutSpecialFunc  (SpecialKeys);
	 glutMouseFunc	  (Mouse);
	 glutIdleFunc	  (DrawGLScene);
	 
	 glutMainLoop();	 
	 return 0;
}