# Spinning Top Simulation

This project visualises the motion of a symmetrical spinning top using Java and LWJGL. The physics is based on the Hamiltonian formulation for a rigid body under gravity. The orientation of the top is described by the Euler angles `theta`, `phi` and `psi`. A fourth‑order Runge–Kutta method with adaptive step size integrates the equations of motion each frame and the result is rendered as a simple 3‑D cone.

## Building

The project uses Maven. Compile and package the application with

```bash
mvn package
```

This produces `target/spinning-top-1.0-SNAPSHOT.jar`. Maven will download LWJGL and its dependencies, including the macOS ARM64 native libraries, when building for the first time.

## Running

Execute the application with

```bash
java -XstartOnFirstThread -jar target/spinning-top-1.0-SNAPSHOT.jar
```

A window opens showing a spinning cone. The arrow keys are not used; the simulation simply runs with the configured parameters. To experiment with different motions edit the initial conditions in `SpinningTopSimulation` and rebuild.

LWJGL bundles the required native libraries so no additional setup is needed when running the program.

## Physics background

For a symmetrical top the principal moments of inertia satisfy $I_1 = I_2 \neq I_3$. In a uniform gravitational field the Hamiltonian can be written as

```
H = p_theta^2/(2*I1) + p_phi^2/(2*I1*sin(theta)^2)
    + p_psi^2/(2*I3) - p_phi*p_psi*cos(theta)/(I1*sin(theta)^2)
    + M*g*l*cos(theta)
```
where $p_\theta$, $p_\phi$ and $p_\psi$ are the canonical momenta associated with the Euler angles. The simulation keeps track of these values and updates the orientation accordingly. Energy is conserved so monitoring it is a good way to check the integration accuracy.

Default parameter values such as the mass `M`, distance to the centre of mass `l`, and inertias `J1` and `J3` can be found in `SpinningTopSimulation.resetInitialConditions()`.



## Further reading

Starting from the kinetic energy
$$
T = \tfrac{1}{2}\,(I_1 \omega_1^2 + I_2 \omega_2^2 + I_3 \omega_3^2),
$$
applying a Legendre transform yields the Hamiltonian
$$
H = \frac{p_\theta^2}{2I_1} + \frac{p_\phi^2}{2I_1 \sin^2\theta} + \frac{p_\psi^2}{2I_3} - \frac{p_\phi p_\psi \cos\theta}{I_1 \sin^2\theta} + M g l \cos\theta.
$$
The resulting first order system requires the six initial conditions $(\theta,\phi,\psi,p_\theta,p_\phi,p_\psi)$.

To integrate these equations the document proposes the classical fourth-order Runge--Kutta method,

$$
\begin{aligned}
  k_1 &= h\,f(t_n, y_n),\\
  k_2 &= h\,f\bigl(t_n + \tfrac{h}{2},\; y_n + \tfrac{k_1}{2}\bigr),\\
  k_3 &= h\,f\bigl(t_n + \tfrac{h}{2},\; y_n + \tfrac{k_2}{2}\bigr),\\
  k_4 &= h\,f(t_n + h,\; y_n + k_3),\\
  y_{n+1} &= y_n + \tfrac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4).
\end{aligned}
$$
Adaptive step size control is recommended so that smaller steps are taken when the solution changes rapidly.

Example parameter sets used in the PDF include

```
M = 1.05, g = 9.81, I1 = 0.15, I3 = 5.10,
θ0 = 10°,  φ0 = 0°,  ψ0 = 0°,
pθ0 = 0,   pφ0 = 5.30,  pψ0 = 4.80

M = 1.15, g = 9.81, I1 = 0.23, I3 = 5.10,
θ0 = 0.8°, φ0 = 24°, ψ0 = 30°,
pθ0 = 0,   pφ0 = 0.0,  pψ0 = 13.3

M = 1.65, g = 9.81, I1 = 0.75, I3 = 5.10,
θ0 = 0.5°, φ0 = 0°,  ψ0 = 0°,
pθ0 = 35,  pφ0 = 0.2,  pψ0 = 4.4
```
