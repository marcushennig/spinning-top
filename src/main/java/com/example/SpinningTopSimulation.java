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
    
    private final int itermax = 10000;
    /** Desired accuracy for the integrator. */
    private final double accuracy = 1E-6;
    /** Integration time step used by the viewer. */
    private final double dt = 0.0001;

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
        var sinTheta = Math.sin(state.theta);
        var cosTheta = Math.cos(state.theta);

        var kineticTheta = state.pTheta * state.pTheta / (2 * constants.inertiaPerpendicular);
        var kineticPsi = state.pPsi * state.pPsi / (2 * constants.inertiaAxis);
        var potential = constants.mass * constants.g * constants.l * cosTheta;
        var phiTerm = state.pPhi - state.pPsi * cosTheta;
        var kineticPhi = (phiTerm * phiTerm) / (2 * constants.inertiaPerpendicular * sinTheta * sinTheta);

        state.energy = kineticTheta + kineticPsi + kineticPhi + potential;
    }

    /**
     * Updates the canonical momenta based on the current angular velocities.
     */
    public void updateAngularMomentum() {
        double sinTheta = Math.sin(state.theta);
        double cosTheta = Math.cos(state.theta);
        double sinThetaSq = sinTheta * sinTheta;
        double cosThetaSq = cosTheta * cosTheta;

        // pTheta: conjugate momentum for theta
        state.pTheta = constants.inertiaPerpendicular * state.thetaDot;

        // pPhi: conjugate momentum for phi
        state.pPhi = constants.inertiaAxis * cosTheta * state.psiDot
            + (constants.inertiaPerpendicular * sinThetaSq + constants.inertiaAxis * cosThetaSq) * state.phiDot;

        // pPsi: conjugate momentum for psi
        state.pPsi = constants.inertiaAxis * (cosTheta * state.phiDot + state.psiDot);
    }

    /**
     * Updates the angular velocities from the current momenta.
     */
    public void updateAngularVelocity() {
        // Compute sin and cos of theta once for readability
        var sinTheta = Math.sin(state.theta);
        var cosTheta = Math.cos(state.theta);
        var sinThetaSq = sinTheta * sinTheta;

        // Angular velocity thetaDot
        state.thetaDot = state.pTheta / constants.inertiaPerpendicular;

        // Angular velocity phiDot
        state.phiDot = (state.pPhi - state.pPsi * cosTheta) / (sinThetaSq * constants.inertiaPerpendicular);

        // Angular velocity psiDot
        var numerator = constants.inertiaPerpendicular * state.pPsi
                 + constants.inertiaAxis * state.pPsi * cosTheta * cosTheta
                 - constants.inertiaAxis * state.pPhi * cosTheta;
        var denominator = constants.inertiaPerpendicular * constants.inertiaAxis * sinThetaSq;
        state.psiDot = numerator / denominator;
    }

    /**
     * Right hand side of the Hamiltonian equations.
     *
     * @param x current state vector
     * @return time derivative of {@code x}
     */
    public RealVector F(RealVector x){
        // Unpack state vector
        var theta = x.getEntry(0);    
        var pTheta = x.getEntry(3);

        // Precompute trigonometric terms
        var sinTheta = Math.sin(theta);
        var cosTheta = Math.cos(theta);
        var sinThetaSq = sinTheta * sinTheta;

        // Use current state for pPhi and pPsi (constants of motion in this reduced system)
        var pPhi = state.pPhi;
        var pPsi = state.pPsi;

        // Auxiliary term for compactness
        var phiTerm = pPhi - pPsi * cosTheta;

        // Hamiltonian equations
        var dTheta = pTheta / constants.inertiaPerpendicular;
        var dPhi = phiTerm / (constants.inertiaPerpendicular * sinThetaSq);
        var dPsi = pPsi / constants.inertiaAxis - cosTheta * phiTerm / (constants.inertiaPerpendicular * sinThetaSq);

        var dPTheta =
            phiTerm * (
            -pPsi / (constants.inertiaPerpendicular * sinTheta)
            + phiTerm * cosTheta / (constants.inertiaPerpendicular * sinTheta * sinTheta * sinTheta)
            )
            + constants.mass * constants.g * constants.l * sinTheta;

        return new ArrayRealVector(new double[] { dTheta, dPhi, dPsi, dPTheta });
    }

    /**
     * Advances the state vector using a single fourth order Runge–Kutta step.
     *
     * @param x current state vector
     * @param h integration step size
     * @return increment to be added to {@code x}
     */
    public RealVector step(RealVector x, double h){
        // Compute Runge-Kutta increments
        var k1 = F(x);
        var k2 = F(x.add(k1.mapMultiply(h / 2.0)));
        var k3 = F(x.add(k2.mapMultiply(h / 2.0)));
        var k4 = F(x.add(k3.mapMultiply(h)));

        // Weighted sum of increments
        var increment = k1
            .add(k2.mapMultiply(2.0))
            .add(k3.mapMultiply(2.0))
            .add(k4)
            .mapMultiply(h / 6.0);

        return increment;
    }

    /**
     * Integrates the equations of motion between two time points.
     *
     * @param a start time
     * @param b end time
     */
    public void evolute(double a,double b){
        // Initialize state vector and integration parameters
        var x = new ArrayRealVector(new double[] {
            state.theta, state.phi, state.psi, state.pTheta
        });
        var t = a;
        var totalTime = b - a;
        var h = totalTime / 10.0;
        var hMax = totalTime / 2.0;
        var eps = 1E-12;
        var iter = 0;

        // Adaptive Runge-Kutta integration loop
        while (t < b && iter < itermax) {
            iter++;
            h = Math.min(h, b - t); // Ensure we don't step past the end
            var H = h;

            // Two half-steps
            var u = x.add(step(x, H / 2));
            u = u.add(step(u, H / 2));

            // One full step
            var v = x.add(step(x, H));

            // Estimate error
            var diff = u.subtract(v);
            var error = Math.sqrt(diff.dotProduct(diff)) / 15.0;

            // Accept step if error is within tolerance
            if (error < accuracy) {
                t += H;
                x = u;
            }

            // Adjust step size for next iteration
            if (error >= eps) {
                h = Math.min(
                    2 * H,
                    Math.min(0.9 * H * Math.pow(accuracy / error, 0.2), hMax)
                );
            } else {
                h = Math.min(2 * H, hMax);
            }
            if (h < eps) {
                h = 2 * eps;
            }
        }

        // Update state with integrated values
        state.theta = x.getEntry(0);
        state.phi = x.getEntry(1);
        state.psi = x.getEntry(2);
        state.pTheta = x.getEntry(3);
    }

    /**
     * Restores a set of default parameters and state.
     */
    public void resetInitialConditions() {

        // Set physical constants to typical spinning top values
        constants.inertiaPerpendicular = 0.75;   // Moment of inertia about axis 1 (kg·m²)
        constants.inertiaAxis = 5.10;   // Moment of inertia about axis 3 (kg·m²)
        constants.mass = 1.65;    // Mass (kg)
        constants.g = 9.81;    // Gravitational acceleration (m/s²)
        constants.l = 20;    // Distance from pivot to center of mass (m)

        // Set initial angles (radians)
        state.theta = 0.5 * Math.PI; // 90 degrees (horizontal)
        state.phi = 0.0;
        state.psi = 0.0;

        // Set initial angular velocities (rad/s)
        state.thetaDot = 30.0;
        state.phiDot = 22.0;
        state.psiDot = 56.0;

        // Update momenta and energy to match initial conditions
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
