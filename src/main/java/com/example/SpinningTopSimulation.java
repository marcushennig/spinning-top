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
    private final double dt = 0.001;

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
        double sinTheta = Math.sin(state.theta);
        double cosTheta = Math.cos(state.theta);

        double kineticTheta = state.pTheta * state.pTheta / (2 * constants.J1);
        double kineticPsi = state.pPsi * state.pPsi / (2 * constants.J3);
        double potential = constants.M * constants.g * constants.l * cosTheta;
        double phiTerm = state.pPhi - state.pPsi * cosTheta;
        double kineticPhi = (phiTerm * phiTerm) / (2 * constants.J1 * sinTheta * sinTheta);

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
        state.pTheta = constants.J1 * state.thetaDot;

        // pPhi: conjugate momentum for phi
        state.pPhi = constants.J3 * cosTheta * state.psiDot
            + (constants.J1 * sinThetaSq + constants.J3 * cosThetaSq) * state.phiDot;

        // pPsi: conjugate momentum for psi
        state.pPsi = constants.J3 * (cosTheta * state.phiDot + state.psiDot);
    }

    /**
     * Updates the angular velocities from the current momenta.
     */
    public void updateAngularVelocity() {
        // Compute sin and cos of theta once for readability
        double sinTheta = Math.sin(state.theta);
        double cosTheta = Math.cos(state.theta);
        double sinThetaSq = sinTheta * sinTheta;

        // Angular velocity thetaDot
        state.thetaDot = state.pTheta / constants.J1;

        // Angular velocity phiDot
        state.phiDot = (state.pPhi - state.pPsi * cosTheta) / (sinThetaSq * constants.J1);

        // Angular velocity psiDot
        double numerator = constants.J1 * state.pPsi
                 + constants.J3 * state.pPsi * cosTheta * cosTheta
                 - constants.J3 * state.pPhi * cosTheta;
        double denominator = constants.J1 * constants.J3 * sinThetaSq;
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
        double theta = x.getEntry(0);    
        double pTheta = x.getEntry(3);

        // Precompute trigonometric terms
        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);
        double sinThetaSq = sinTheta * sinTheta;

        // Use current state for pPhi and pPsi (constants of motion in this reduced system)
        double pPhi = state.pPhi;
        double pPsi = state.pPsi;

        // Auxiliary term for compactness
        double phiTerm = pPhi - pPsi * cosTheta;

        // Hamiltonian equations
        double dTheta = pTheta / constants.J1;
        double dPhi = phiTerm / (constants.J1 * sinThetaSq);
        double dPsi = pPsi / constants.J3 - cosTheta * phiTerm / (constants.J1 * sinThetaSq);

        double dPTheta =
            phiTerm * (
            -pPsi / (constants.J1 * sinTheta)
            + phiTerm * cosTheta / (constants.J1 * sinTheta * sinTheta * sinTheta)
            )
            + constants.M * constants.g * constants.l * sinTheta;

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
        RealVector k1 = F(x);
        RealVector k2 = F(x.add(k1.mapMultiply(h / 2.0)));
        RealVector k3 = F(x.add(k2.mapMultiply(h / 2.0)));
        RealVector k4 = F(x.add(k3.mapMultiply(h)));

        // Weighted sum of increments
        RealVector increment = k1
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
    public void NDSolve(double a,double b){
        // Initialize state vector and integration parameters
        RealVector x = new ArrayRealVector(new double[] {
            state.theta, state.phi, state.psi, state.pTheta
        });
        double t = a;
        double totalTime = b - a;
        double h = totalTime / 10.0;
        double hMax = totalTime / 2.0;
        double eps = 1E-12;
        int iter = 0;

        // Adaptive Runge-Kutta integration loop
        while (t < b && iter < itermax) {
            iter++;
            h = Math.min(h, b - t); // Ensure we don't step past the end
            double H = h;

            // Two half-steps
            RealVector u = x.add(step(x, H / 2));
            u = u.add(step(u, H / 2));

            // One full step
            RealVector v = x.add(step(x, H));

            // Estimate error
            RealVector diff = u.subtract(v);
            double error = Math.sqrt(diff.dotProduct(diff)) / 15.0;

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
        constants.J1 = 0.75;   // Moment of inertia about axis 1 (kg·m²)
        constants.J3 = 5.10;   // Moment of inertia about axis 3 (kg·m²)
        constants.M = 1.65;    // Mass (kg)
        constants.g = 9.81;    // Gravitational acceleration (m/s²)
        constants.l = 0.80;    // Distance from pivot to center of mass (m)

        // Set initial angles (radians)
        state.theta = 0.5 * Math.PI; // 90 degrees (horizontal)
        state.phi = 0.0;
        state.psi = 0.0;

        // Set initial angular velocities (rad/s)
        state.thetaDot = 0.0;
        state.phiDot = 0.0;
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
