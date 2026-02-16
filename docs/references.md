# FLIP/APIC Reference Papers

This solver structure and terminology follow these references that are commonly used in graphics/CFD-oriented fluid solvers:

1. **Jiang, C., Schroeder, C., Selle, A., Teran, J., and Stomakhin, A. (2015).**
   *The Affine Particle-In-Cell Method.* ACM Transactions on Graphics (SIGGRAPH).

2. **Zhu, Y., and Bridson, R. (2005).**
   *Animating Sand as a Fluid.* ACM SIGGRAPH/Eurographics Symposium on Computer Animation.
   - The production FLIP style used in DCC tools is heavily based on this family of particle-grid methods.

3. **Bridson, R. (2015).**
   *Fluid Simulation for Computer Graphics, Second Edition.* CRC Press.
   - Practical pressure projection and boundary treatment details.

## Notes for this codebase

- Particle-to-grid transfer: APIC affine momentum transfer.
- Grid pressure solve: Poisson solve of incompressibility equation.
- Grid-to-particle transfer: FLIP/PIC blend (tunable ratio).
- Particle advection: midpoint (RK2-style) advection based on sampled grid velocity.
