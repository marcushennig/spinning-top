package com.example;

/**
 * Holds the dynamic state of the spinning top.
 * This includes the Euler angles and their
 * associated canonical momenta.
 */
public class SpinningTopState {
    /** Nutation angle. */
    public double theta;
    /** Precession angle. */
    public double phi;
    /** Spin angle. */
    public double psi;

    /** Momentum conjugate to {@code theta}. */
    public double pTheta;
    /** Momentum conjugate to {@code phi}. */
    public double pPhi;
    /** Momentum conjugate to {@code psi}. */
    public double pPsi;

    /** Angular velocity about \u03b8. */
    public double thetaDot;
    /** Angular velocity about \u03c6. */
    public double phiDot;
    /** Angular velocity about \u03c8. */
    public double psiDot;

    /** Current Hamiltonian energy of the system. */
    public double energy;
}
