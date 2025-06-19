package com.example;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

/**
 * Numerical model for the spinning top dynamics.
 * This class implements a simple symplectic integrator
 * for a spinning top using canonical coordinates.
 */
public class SpinningTopSimulation {
    /** Dynamic state of the spinning top. */
    private SpinningTopState state = new SpinningTopState();
    /** Physical constants used in the model. */
    private PhysicalConstants constants = new PhysicalConstants();
    /** Maximum number of iterations for the integrator. */
    private int itermax;
    /** Desired accuracy for the integrator. */
    private double accuracy;
    /** Integration time step used by the viewer. */
    private double dt;

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
        state.energy = state.pTheta*state.pTheta/(2*constants.J1)
                + constants.M*constants.g*constants.l*Math.cos(state.theta)
                + state.pPsi*state.pPsi/(2*constants.J3)
                + Math.pow((state.pPhi - state.pPsi*Math.cos(state.theta)), 2.0)
                  /(2*constants.J1*Math.pow(Math.sin(state.theta),2.0));
    }

    /**
     * Updates the canonical momenta based on the current angular velocities.
     */
    public void updateAngularMomentum() {
        state.pTheta = constants.J1 * state.thetaDot;
        state.pPhi = constants.J3*Math.cos(state.theta)*state.psiDot
                + (constants.J1*Math.pow(Math.sin(state.theta),2)
                   + constants.J3*Math.pow(Math.cos(state.theta),2))*state.phiDot;
        state.pPsi = constants.J3*(Math.cos(state.theta)*state.phiDot + state.psiDot);
    }

    /**
     * Updates the angular velocities from the current momenta.
     */
    public void updateAngularVelocity() {
        state.thetaDot = state.pTheta/constants.J1;
        state.phiDot = (state.pPhi - state.pPsi*Math.cos(state.theta))
                /(Math.pow(Math.sin(state.theta),2)*constants.J1);
        state.psiDot = (constants.J1*state.pPsi + constants.J3*state.pPsi*Math.pow(Math.cos(state.theta),2)
                - constants.J3*state.pPhi*Math.cos(state.theta))
                /(constants.J1*constants.J3*Math.pow(Math.sin(state.theta),2));
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
        double A = state.pPhi - state.pPsi * c_theta;

        return new ArrayRealVector(new double[]{
            p_theta / constants.J1,
            A / (constants.J1*s_theta*s_theta),
            state.pPsi / constants.J3 - c_theta*A/(constants.J1*s_theta*s_theta),
            A * (-state.pPsi/(constants.J1*s_theta) + A * c_theta/(constants.J1*s_theta*s_theta*s_theta))
                    + constants.M*constants.g*constants.l*s_theta
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
        RealVector x = new ArrayRealVector(new double[]{state.theta,state.phi,state.psi,state.pTheta});
        double t=a;
        double h=(b-a)/10.0;
        double hmax=(b-a)/2.0;
        double eps=1E-12;
        int iter=0;
        while(t<b && iter<itermax){
            iter++; h=Math.min(h,b-t);
            double H=h;
            RealVector u=x.add(step(x,H/2));
            u=u.add(step(u,H/2));
            RealVector v=x.add(step(x,H));
            RealVector d=u.subtract(v);
            double error=Math.sqrt(d.dotProduct(d))/15.0;
            if(error < accuracy){ t=t+H; x=u; }
            if(error>=eps) h=Math.min(2*H, Math.min(0.9*H*Math.pow(accuracy/error,0.2),hmax));
            else h=Math.min(2*H, hmax);
            if(h<eps) h=2*eps;
        }
        state.theta=x.getEntry(0);
        state.phi=x.getEntry(1);
        state.psi=x.getEntry(2);
        state.pTheta=x.getEntry(3);
    }

    /**
     * Restores a set of default parameters and state.
     */
    public void resetInitialConditions() {
        constants.J1 = 0.75;
        constants.J3 = 5.10;
        constants.M = 1.65;
        constants.g = 9.81;
        constants.l = 80;

        accuracy = 1E-4;
        itermax = 10000;

        state.theta = 0.5*Math.PI;
        state.phi = 0;
        state.psi = 0;

        state.thetaDot = 0.0;
        state.phiDot = 0.0;
        state.psiDot = 56;

        dt = 0.001;

        updateAngularMomentum();
        updateEnergy();
    }

    /** Returns the current dynamic state. */
    public SpinningTopState getState() { return state; }

    /** Returns the simulation constants. */
    public PhysicalConstants getConstants() { return constants; }

    /** Returns the integration time step used by the viewer. */
    public double getDt() { return dt; }
}
