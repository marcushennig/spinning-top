package com.example;

/** Numerical model for the spinning top dynamics. */
public class SpinningTopSimulation {
    private Parameter c = new Parameter();

    public SpinningTopSimulation() {
        resetInitialConditions();
    }

    private double norm(double x) {
        return x < 0 ? -x : x;
    }

    public void updateEnergy() {
        c.E = c.p_theta*c.p_theta/(2*c.J1) + c.M*c.g*c.l*Math.cos(c.theta)
                + c.p_psi*c.p_psi/(2*c.J3)
                + Math.pow((c.p_phi - c.p_psi*Math.cos(c.theta)),2.0)
                  /(2*c.J1*Math.pow(Math.sin(c.theta),2.0));
    }

    public void updateAngularMomentum() {
        c.p_theta = c.J1 * c.Dtheta;
        c.p_phi = c.J3*Math.cos(c.theta)*c.Dpsi
                + (c.J1*Math.pow(Math.sin(c.theta),2)+c.J3*Math.pow(Math.cos(c.theta),2))*c.Dphi;
        c.p_psi = c.J3*(Math.cos(c.theta)*c.Dphi + c.Dpsi);
    }

    public void updateAngularVelocity() {
        c.Dtheta = c.p_theta/c.J1;
        c.Dphi = (c.p_phi - c.p_psi*Math.cos(c.theta))/(Math.pow(Math.sin(c.theta),2)*c.J1);
        c.Dpsi = (c.J1*c.p_psi + c.J3*c.p_psi*Math.pow(Math.cos(c.theta),2) - c.J3*c.p_phi*Math.cos(c.theta))
                /(c.J1*c.J3*Math.pow(Math.sin(c.theta),2));
    }

    public Vector4d F(Vector4d x){
        double theta = x.get(1);
        double phi = x.get(2);
        double psi = x.get(3);
        double p_theta = x.get(4);
        double s_theta = Math.sin(theta);
        double c_theta = Math.cos(theta);
        double A = c.p_phi - c.p_psi * c_theta;

        return new Vector4d(
            p_theta / c.J1,
            A / (c.J1*s_theta*s_theta),
            c.p_psi / c.J3 - c_theta*A/(c.J1*s_theta*s_theta),
            A * (-c.p_psi/(c.J1*s_theta) + A * c_theta/(c.J1*s_theta*s_theta*s_theta)) + c.M*c.g*c.l*s_theta
        );
    }

    public Vector4d step(Vector4d x, double h){
        Vector4d k1 = F(x).scale(h);
        Vector4d k2 = F(x.add(k1.scale(0.5))).scale(h);
        Vector4d k3 = F(x.add(k2.scale(0.5))).scale(h);
        Vector4d k4 = F(x.add(k3)).scale(h);
        return k1.add(k2.scale(2)).add(k3.scale(2)).add(k4).div(6.0);
    }

    public void NDSolve(double a,double b){
        Vector4d x = new Vector4d(c.theta,c.phi,c.psi,c.p_theta);
        double t=a;
        double h=(b-a)/10.0;
        double hmax=(b-a)/2.0;
        double eps=1E-12;
        int iter=0;
        while(t<b && iter<c.itermax){
            iter++; h=Math.min(h,b-t);
            double H=h;
            Vector4d u=x.add(step(x,H/2));
            u=u.add(step(u,H/2));
            Vector4d v=x.add(step(x,H));
            Vector4d d=u.sub(v);
            double error=Math.sqrt(d.dot(d))/15.0;
            if(error < c.accuracy){ t=t+H; x=u; }
            if(error>=eps) h=Math.min(2*H, Math.min(0.9*H*Math.pow(c.accuracy/error,0.2),hmax));
            else h=Math.min(2*H, hmax);
            if(h<eps) h=2*eps;
        }
        c.theta=x.get(1);
        c.phi=x.get(2);
        c.psi=x.get(3);
        c.p_theta=x.get(4);
    }

    public void resetInitialConditions() {
        c.J1=0.75;
        c.J3=5.10;
        c.M=1.65;
        c.g=9.81;
        c.l=80;
        c.accuracy=1E-4;
        c.itermax=10000;
        c.theta=0.5*Math.PI;
        c.phi=0;
        c.psi=0;
        c.Dtheta=0.0;
        c.Dphi=0.0;
        c.Dpsi=56;
        c.dt=0.001;
        updateAngularMomentum();
        updateEnergy();
    }

    public Parameter getParameter(){ return c; }
}
