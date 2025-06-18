# Spinning Top Simulation

This project visualises the motion of a symmetrical spinning top using Java and JOGL. The physics is based on the Hamiltonian formulation for a rigid body under gravity. The orientation of the top is described by the Euler angles `theta`, `phi` and `psi`. A fourth‑order Runge–Kutta method with adaptive step size integrates the equations of motion each frame and the result is rendered as a simple 3‑D cone.

## Physics background

For a symmetrical top the principal moments of inertia satisfy $I_1 = I_2 \neq I_3$. In a uniform gravitational field the Hamiltonian can be written as

```
H = p_theta^2/(2*I1) + p_phi^2/(2*I1*sin(theta)^2)
    + p_psi^2/(2*I3) - p_phi*p_psi*cos(theta)/(I1*sin(theta)^2)
    + M*g*l*cos(theta)
```
where $p_\theta$, $p_\phi$ and $p_\psi$ are the canonical momenta associated with the Euler angles. The simulation keeps track of these values and updates the orientation accordingly. Energy is conserved so monitoring it is a good way to check the integration accuracy.

Default parameter values such as the mass `M`, distance to the centre of mass `l`, and inertias `J1` and `J3` can be found in `SpinningTopSimulation.resetInitialConditions()`.

## Building

The project uses Maven. Compile and package the application with

```bash
mvn package
```

This produces `target/spinning-top-1.0-SNAPSHOT.jar`. Maven will download JOGL and its dependencies when building for the first time.

## Running

Execute the application with

```bash
java -jar target/spinning-top-1.0-SNAPSHOT.jar
```

A window opens showing a spinning cone. The arrow keys are not used; the simulation simply runs with the configured parameters. To experiment with different motions edit the initial conditions in `SpinningTopSimulation` and rebuild.

Make sure the JOGL native libraries are available on your system when running the program.

## Further reading

The PDF `A Symmetrical Spinning Top.pdf` in this repository contains a short derivation of the equations of motion and explains the numerical method in more detail.

