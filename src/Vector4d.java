public class Vector4d {
    public double x;
    public double y;
    public double z;
    public double w;

    public Vector4d() {
        this(0,0,0,0);
    }
    public Vector4d(double a,double b,double c,double d) {
        x=a; y=b; z=c; w=d;
    }
    public Vector4d(Vector4d v){
        this(v.x,v.y,v.z,v.w);
    }

    public double get(int i) {
        switch(i) {
            case 1: return x;
            case 2: return y;
            case 3: return z;
            case 4: return w;
            default: return 0;
        }
    }

    public Vector4d add(Vector4d v) { return new Vector4d(x+v.x,y+v.y,z+v.z,w+v.w); }
    public Vector4d sub(Vector4d v) { return new Vector4d(x-v.x,y-v.y,z-v.z,w-v.w); }
    public double dot(Vector4d v) { return x*v.x+y*v.y+z*v.z+w*v.w; }
    public Vector4d scale(double s){ return new Vector4d(x*s,y*s,z*s,w*s); }
    public Vector4d div(double s){ return new Vector4d(x/s,y/s,z/s,w/s); }

    public Vector4d plusAssign(Vector4d v){ x+=v.x; y+=v.y; z+=v.z; w+=v.w; return this; }
    public Vector4d timesAssign(double s){ x*=s; y*=s; z*=s; w*=s; return this; }
}
