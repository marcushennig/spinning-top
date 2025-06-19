package com.example;

/**
 * Simple 4D vector utility used by the simulation.
 */
public class Vector4d {
    public double x;
    public double y;
    public double z;
    public double w;

    /** Creates a zero vector. */
    public Vector4d() {
        this(0,0,0,0);
    }

    /**
     * Creates a vector with the supplied components.
     *
     * @param a first component
     * @param b second component
     * @param c third component
     * @param d fourth component
     */
    public Vector4d(double a,double b,double c,double d) {
        x=a; y=b; z=c; w=d;
    }

    /**
     * Copy constructor.
     *
     * @param v vector to copy
     */
    public Vector4d(Vector4d v){
        this(v.x,v.y,v.z,v.w);
    }

    /**
     * Returns the i-th component (1 based).
     *
     * @param i index in range 1–4
     * @return component value or 0 if out of range
     */
    public double get(int i) {
        switch(i) {
            case 1: return x;
            case 2: return y;
            case 3: return z;
            case 4: return w;
            default: return 0;
        }
    }

    /** Adds two vectors. */
    public Vector4d add(Vector4d v) { return new Vector4d(x+v.x,y+v.y,z+v.z,w+v.w); }
    /** Subtracts the supplied vector from this one. */
    public Vector4d sub(Vector4d v) { return new Vector4d(x-v.x,y-v.y,z-v.z,w-v.w); }
    /** Dot product of two vectors. */
    public double dot(Vector4d v) { return x*v.x+y*v.y+z*v.z+w*v.w; }
    /** Scales all components by {@code s}. */
    public Vector4d scale(double s){ return new Vector4d(x*s,y*s,z*s,w*s); }
    /** Divides all components by {@code s}. */
    public Vector4d div(double s){ return new Vector4d(x/s,y/s,z/s,w/s); }

    /** Adds another vector to this one in place. */
    public Vector4d plusAssign(Vector4d v){ x+=v.x; y+=v.y; z+=v.z; w+=v.w; return this; }
    /** Scales this vector in place. */
    public Vector4d timesAssign(double s){ x*=s; y*=s; z*=s; w*=s; return this; }
}
