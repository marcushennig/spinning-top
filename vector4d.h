/*
 *  vector3d.h
 *  SpinningTop
 *
 *  Created by hennig on 10/21/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>

class vector4d
{
	 
private:
	 double x;
	 double y;
	 double z;
	 double w;
	 
public:
	 vector4d ( );
	  vector4d ( const vector4d& );
	 vector4d ( double, double, double, double );
	 vector4d ( double[] );
      ~vector4d( );
		
	 double operator[](int);
		
	 vector4d& operator=(const vector4d&);
	 	 
	 friend vector4d operator+( vector4d, vector4d );
	 friend vector4d operator- ( vector4d, vector4d );
	 friend double    operator*( vector4d, vector4d );
	 
	 vector4d& operator+=( vector4d );
	 vector4d& operator -=( vector4d );
	 
	 vector4d& operator++();
	 vector4d operator++(int);
	 
	 friend vector4d operator*( double, vector4d );
	 friend vector4d operator*( vector4d, double );
	 
	 vector4d& operator*=( double );
	 
	 friend vector4d operator/( vector4d, double );
	 vector4d& operator/=( double );
};

