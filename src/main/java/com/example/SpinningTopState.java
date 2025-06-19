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
    public double p_theta;
    /** Momentum conjugate to {@code phi}. */
    public double p_phi;
    /** Momentum conjugate to {@code psi}. */
    public double p_psi;
}
