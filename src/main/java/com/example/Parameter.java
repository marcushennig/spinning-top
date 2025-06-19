package com.example;

/**
 * Configuration parameters for the spinning top simulation.
 * It bundles together the current {@link SpinningTopState}
 * and the {@link PhysicalConstants} used in the model.
 */
public class Parameter {
    /** Current dynamic state. */
    public SpinningTopState state = new SpinningTopState();
    /** Physical constants of the system. */
    public PhysicalConstants constants = new PhysicalConstants();

    /** Angular velocity about \u03b8. */
    public double thetaDot;
    /** Angular velocity about \u03c6. */
    public double phiDot;
    /** Angular velocity about \u03c8. */
    public double psiDot;

    /** Current energy of the system. */
    public double energy;
    /** Maximum number of iterations for the integrator. */
    public int itermax;
    /** Desired accuracy for the integrator. */
    public double accuracy;
    /** Time step used by the viewer. */
    public double dt;
}
