/*
 *  vector3d.cpp
 *  SpinningTop
 *
 *  Created by hennig on 10/21/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "vector4d.h"

vector4d::vector4d()
{
	 x = 0.0;
	 y = 0.0;
	 z = 0.0;
	 w = 0.0;
}

vector4d::vector4d( const vector4d& v)
{
	 x  = v.x;
	 y  = v.y;
	 z  = v.z;
	 w = v.w;
}

vector4d::vector4d( double a, double b, double c, double d )
{
	 x  = a;
	 y  = b;
	 z  = c;
	 w = d;
}

vector4d::vector4d( double* v )
{
	 x  = v[0];
	 y  = v[1];
	 z  = v[2];
	 w = v[3];
}

vector4d::~vector4d()
{}

double vector4d::operator[](int n)
{
	 switch( n )
	 {
		  case  1: return x;
		  case  2: return y;
		  case  3: return z;
		  case  4: return w;
		  default: return 0;
	 };
}

vector4d& vector4d::operator=(const vector4d& v)
{ 
	 x  = v.x;
	 y  = v.y;
	 z  = v.z;
	 w = v.w;
	 
	 return *this;
}

vector4d operator+( vector4d v1, vector4d v2 )
{
	 return vector4d(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w );
}

vector4d operator-( vector4d v1, vector4d v2 )
{
	 return vector4d( v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w );
}

double operator*( vector4d v1, vector4d v2 )
{
	 return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
}

vector4d& vector4d::operator+=( vector4d v )
{
	 x  += v.x;
	 y  += v.y;
	 z  += v.z;
	 w += v.w;
	 return *this;
}

vector4d& vector4d::operator-=( vector4d v )
{
	 x  -= v.x;
	 y  -= v.y;
	 z  -= v.z;
	 w -= v.w;
	 
	 return *this;
}

vector4d& vector4d::operator++()
{
	 x++;
	 y++;
	 z++;
	 w++;
	 
	 return *this;
}

vector4d vector4d::operator++(int)
{
	 vector4d temp = *this;
	 
	 x++;
	 y++;
	 z++;
	 w++;
	 
	 return temp;
}

vector4d operator*( vector4d v, double s )
{
	 return vector4d(v.x*s, v.y*s, v.z*s, v.w*s);
}

vector4d operator*( double s, vector4d v )
{
	 return vector4d(v.x*s, v.y*s, v.z*s, v.w*s);
}

vector4d& vector4d::operator*=( double a )
{
	 x  *= a;
	 y  *= a;
	 z  *= a;
	 w *= a;
	 
	 return *this;
}

vector4d operator/( vector4d v, double s )
{
	 return vector4d(v.x/s, v.y/s, v.z/s, v.w/s);
}
