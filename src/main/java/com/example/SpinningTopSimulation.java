package com.example;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

/**
 * Numerical model for the spinning top dynamics.
 * This class implements a simple symplectic integrator
 * for a spinning top using canonical coordinates.
 */
public class SpinningTopSimulation {
    /** Simulation parameters and current state. */
    private Parameter c = new Parameter();

    /**
     * Constructs a new simulation instance and initialises the
     * state with default values.
     */
    public SpinningTopSimulation() {
        resetInitialConditions();
    }

    /**
     * Recomputes the Hamiltonian energy from the current state.
     */
    public void updateEnergy() {
        c.E = c.state.p_theta*c.state.p_theta/(2*c.constants.J1)
                + c.constants.M*c.constants.g*c.constants.l*Math.cos(c.state.theta)
                + c.state.p_psi*c.state.p_psi/(2*c.constants.J3)
                + Math.pow((c.state.p_phi - c.state.p_psi*Math.cos(c.state.theta)), 2.0)
                  /(2*c.constants.J1*Math.pow(Math.sin(c.state.theta),2.0));
    }

    /**
     * Updates the canonical momenta based on the current angular velocities.
     */
    public void updateAngularMomentum() {
        c.state.p_theta = c.constants.J1 * c.Dtheta;
        c.state.p_phi = c.constants.J3*Math.cos(c.state.theta)*c.Dpsi
                + (c.constants.J1*Math.pow(Math.sin(c.state.theta),2)
                   + c.constants.J3*Math.pow(Math.cos(c.state.theta),2))*c.Dphi;
        c.state.p_psi = c.constants.J3*(Math.cos(c.state.theta)*c.Dphi + c.Dpsi);
    }

    /**
     * Updates the angular velocities from the current momenta.
     */
    public void updateAngularVelocity() {
        c.Dtheta = c.state.p_theta/c.constants.J1;
        c.Dphi = (c.state.p_phi - c.state.p_psi*Math.cos(c.state.theta))
                /(Math.pow(Math.sin(c.state.theta),2)*c.constants.J1);
        c.Dpsi = (c.constants.J1*c.state.p_psi + c.constants.J3*c.state.p_psi*Math.pow(Math.cos(c.state.theta),2)
                - c.constants.J3*c.state.p_phi*Math.cos(c.state.theta))
                /(c.constants.J1*c.constants.J3*Math.pow(Math.sin(c.state.theta),2));
    }

    /**
     * Right hand side of the Hamiltonian equations.
     *
     * @param x current state vector
     * @return time derivative of {@code x}
     */
    public RealVector F(RealVector x){
        double theta = x.getEntry(0);
        double p_theta = x.getEntry(3);
        double s_theta = Math.sin(theta);
        double c_theta = Math.cos(theta);
        double A = c.state.p_phi - c.state.p_psi * c_theta;

        return new ArrayRealVector(new double[]{
            p_theta / c.constants.J1,
            A / (c.constants.J1*s_theta*s_theta),
            c.state.p_psi / c.constants.J3 - c_theta*A/(c.constants.J1*s_theta*s_theta),
            A * (-c.state.p_psi/(c.constants.J1*s_theta) + A * c_theta/(c.constants.J1*s_theta*s_theta*s_theta))
                    + c.constants.M*c.constants.g*c.constants.l*s_theta
        });
    }

    /**
     * Advances the state vector using a single fourth order Runge–Kutta step.
     *
     * @param x current state vector
     * @param h integration step size
     * @return increment to be added to {@code x}
     */
    public RealVector step(RealVector x, double h){
        RealVector k1 = F(x).mapMultiply(h);
        RealVector k2 = F(x.add(k1.mapMultiply(0.5))).mapMultiply(h);
        RealVector k3 = F(x.add(k2.mapMultiply(0.5))).mapMultiply(h);
        RealVector k4 = F(x.add(k3)).mapMultiply(h);
        return k1.add(k2.mapMultiply(2)).add(k3.mapMultiply(2)).add(k4).mapDivide(6.0);
    }

    /**
     * Integrates the equations of motion between two time points.
     *
     * @param a start time
     * @param b end time
     */
    public void NDSolve(double a,double b){
        RealVector x = new ArrayRealVector(new double[]{c.state.theta,c.state.phi,c.state.psi,c.state.p_theta});
        double t=a;
        double h=(b-a)/10.0;
        double hmax=(b-a)/2.0;
        double eps=1E-12;
        int iter=0;
        while(t<b && iter<c.itermax){
            iter++; h=Math.min(h,b-t);
            double H=h;
            RealVector u=x.add(step(x,H/2));
            u=u.add(step(u,H/2));
            RealVector v=x.add(step(x,H));
            RealVector d=u.subtract(v);
            double error=Math.sqrt(d.dotProduct(d))/15.0;
            if(error < c.accuracy){ t=t+H; x=u; }
            if(error>=eps) h=Math.min(2*H, Math.min(0.9*H*Math.pow(c.accuracy/error,0.2),hmax));
            else h=Math.min(2*H, hmax);
            if(h<eps) h=2*eps;
        }
        c.state.theta=x.getEntry(0);
        c.state.phi=x.getEntry(1);
        c.state.psi=x.getEntry(2);
        c.state.p_theta=x.getEntry(3);
    }

    /**
     * Restores a set of default parameters and state.
     */
    public void resetInitialConditions() {
        c.constants.J1 = 0.75;
        c.constants.J3 = 5.10;
        c.constants.M = 1.65;
        c.constants.g = 9.81;
        c.constants.l = 80;

        c.accuracy = 1E-4;
        c.itermax = 10000;

        c.state.theta = 0.5*Math.PI;
        c.state.phi = 0;
        c.state.psi = 0;

        c.Dtheta = 0.0;
        c.Dphi = 0.0;
        c.Dpsi = 56;

        c.dt = 0.001;

        updateAngularMomentum();
        updateEnergy();
    }

    /**
     * Exposes the current simulation parameters.
     *
     * @return parameter container
     */
    public Parameter getParameter(){ return c; }
}
