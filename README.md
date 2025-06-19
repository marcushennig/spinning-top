# Spinning Top Simulation

This project visualises the motion of a symmetrical spinning top using Java and **jMonkeyEngine**. The physics is based on the Hamiltonian formulation for a rigid body under gravity. The orientation of the top is described by the Euler angles $\theta$, $\phi$ and $\psi$. A fourth‑order Runge–Kutta method with adaptive step size integrates the equations of motion each frame and the result is rendered as a simple 3‑D cone.

## Building

The project uses Maven. Compile and package the application with

```bash
mvn package
```

This produces `target/spinning-top-1.0-SNAPSHOT.jar`. Maven will download jMonkeyEngine and its dependencies the first time the project is built.

## Running

Execute the application with

```bash
java -jar target/spinning-top-1.0-SNAPSHOT.jar
```

A window opens showing a spinning cone. The arrow keys are not used; the simulation simply runs with the configured parameters. To experiment with different motions edit the initial conditions in `SpinningTopSimulation` and rebuild.

jMonkeyEngine handles the native libraries so no additional setup is needed when running the program.

## Physics background

This simulation models a **symmetrical spinning top** under gravity using the Hamiltonian formalism. The orientation of the top is described by the **Euler angles**:  
- $\theta$ (nutation): angle between the symmetry axis and the vertical  
- $\phi$ (precession): rotation around the vertical axis  
- $\psi$ (spin): rotation around the symmetry axis

The **principal moments of inertia** for a symmetric top satisfy $I_1 = I_2 \neq I_3$, where:
- $I_1, I_2$: moments of inertia perpendicular to the symmetry axis
- $I_3$: moment of inertia about the symmetry axis

### Kinetic Energy

The angular velocity components in the body frame are:
- $\omega_1 = \dot{\theta} \cos\psi + \dot{\phi} \sin\theta \sin\psi$
- $\omega_2 = -\dot{\theta} \sin\psi + \dot{\phi} \sin\theta \cos\psi$
- $\omega_3 = \dot{\psi} + \dot{\phi} \cos\theta$

The **kinetic energy** is:
$$
T = \frac{1}{2} \left[ I_1 (\omega_1^2 + \omega_2^2) + I_3 \omega_3^2 \right]
$$

### Potential Energy

The **potential energy** due to gravity is:
$$
V = M g l \cos\theta
$$
where:
- $M$: mass of the top
- $g$: gravitational acceleration
- $l$: distance from the pivot to the center of mass

### Lagrangian

The **Lagrangian** is:
$$
L = T - V
$$

### Generalized Momenta

The **generalized momenta** (canonical momenta) are:
$$
\begin{aligned}
p_\theta &= \frac{\partial L}{\partial \dot{\theta}} = I_1 \dot{\theta} \\
p_\phi   &= \frac{\partial L}{\partial \dot{\phi}} = I_1 \sin^2\theta\, \dot{\phi} + I_3 \cos\theta\, (\dot{\psi} + \dot{\phi} \cos\theta) \\
p_\psi   &= \frac{\partial L}{\partial \dot{\psi}} = I_3 (\dot{\psi} + \dot{\phi} \cos\theta)
\end{aligned}
$$

### Legendre Transform to the Hamiltonian

The **Hamiltonian** is obtained via the Legendre transform:
$$
H = \sum_i p_i \dot{q}_i - L
$$
where $q_i$ are the generalized coordinates $(\theta, \phi, \psi)$.

Expressing $\dot{\theta}, \dot{\phi}, \dot{\psi}$ in terms of the momenta:
$$
\begin{aligned}
\dot{\theta} &= \frac{p_\theta}{I_1} \\
\dot{\phi} &= \frac{p_\phi - p_\psi \cos\theta}{I_1 \sin^2\theta} \\
\dot{\psi} &= \frac{p_\psi}{I_3} - \frac{\cos\theta}{I_1 \sin^2\theta}(p_\phi - p_\psi \cos\theta)
\end{aligned}
$$

Plugging these into the Legendre transform yields:
$$
\boxed{
H = \frac{p_\theta^2}{2 I_1}
  + \frac{(p_\phi - p_\psi \cos\theta)^2}{2 I_1 \sin^2\theta}
  + \frac{p_\psi^2}{2 I_3}
  + M g l \cos\theta
}
$$

Or, expanded:
$$
H = \frac{p_\theta^2}{2 I_1}
  + \frac{p_\phi^2}{2 I_1 \sin^2\theta}
  + \frac{p_\psi^2}{2 I_3}
  - \frac{p_\phi p_\psi \cos\theta}{I_1 \sin^2\theta}
  + M g l \cos\theta
$$

### Summary of Variables

| Symbol         | Meaning                                      |
|----------------|----------------------------------------------|
| $\theta$       | Nutation angle                               |
| $\phi$         | Precession angle                             |
| $\psi$         | Spin angle                                   |
| $p_\theta$     | Canonical momentum conjugate to $\theta$     |
| $p_\phi$       | Canonical momentum conjugate to $\phi$       |
| $p_\psi$       | Canonical momentum conjugate to $\psi$       |
| $I_1, I_3$     | Principal moments of inertia                 |
| $M$            | Mass of the top                              |
| $g$            | Gravitational acceleration                   |
| $l$            | Distance from pivot to center of mass        |

### Numerical Integration

The simulation integrates the equations of motion using the **classical fourth-order Runge–Kutta method** with adaptive step size:
$$
\begin{aligned}
  k_1 &= h\,f(t_n, y_n) \\
  k_2 &= h\,f\left(t_n + \frac{h}{2},\; y_n + \frac{k_1}{2}\right) \\
  k_3 &= h\,f\left(t_n + \frac{h}{2},\; y_n + \frac{k_2}{2}\right) \\
  k_4 &= h\,f(t_n + h,\; y_n + k_3) \\
  y_{n+1} &= y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\end{aligned}
$$

**Energy conservation** is monitored to check the accuracy of the integration.


This approach provides a physically accurate and numerically stable simulation of a spinning top, visualized in real time.

