# INTERSPACE — Physics and Mathematical Basis (Stage 12)

## Scope and evidence standard

This is a **numerical engineering research prototype**, not an experimentally validated device. Equations derived from established mechanics/electrodynamics support *conditional feasibility* of impact absorption and conversion **if a material, actuator, and electrical circuit meet their constitutive and power constraints**. They do not establish a free-standing field barrier or a means of storing energy in literal empty interatomic gaps. Each module below distinguishes physical law, constitutive assumption, discretization, and validation status.

## A. Geometric dome and mesh

The undeformed hemisphere is parameterized by

\[\mathbf X(\theta,\phi)=R(\sin\phi\cos\theta,\sin\phi\sin\theta,\cos\phi),\quad 0\le\phi\le\pi/2.\]

Its equator lies at \(z=0\). Mesh edges are a discrete network approximation, not a continuum shell constitutive law. The ring is fixed in the current baseline. Its reaction force can carry momentum, but a fixed boundary does **zero mechanical work** because its velocity is zero. Consequently, ring electrical harvesting requires a separate *relative-motion* transducer, not direct extraction from a perfectly rigid fixed node.

## B. Interatomic versus mesoscale interactions

For an ideal pair bond, the Morse energy is

\[V(r)=D_e(1-e^{-a(r-r_0)})^2-D_e,\qquad \mathbf F_i=-\nabla_i\sum V(r_{ij}).\]

This is a classical pair-potential approximation, not a full electronic-structure description. The later dome is a *mesoscale spring network*, not an atomic-resolution molecular dynamics simulation. A linear spring edge has

\[U_e=\tfrac12 k_e(\ell_e-\ell_{e0})^2,\quad \mathbf f_{ij}=\left[k_e(\ell_e-\ell_{e0})+c_e\dot\ell_e\right]\hat{\mathbf n}_{ij}.\]

The dashpot dissipation rate is \(P_{\rm heat}=\sum_e c_e\dot\ell_e^2\ge0\). Assigning arbitrary stiffness or damping to an electric field is **not** a derived microscopic interaction.

## C. 3D dynamics and contact

For each mobile node \(i\), Newton's law is

\[m_i\ddot{\mathbf x}_i=\mathbf f_i^{\rm internal}+\mathbf f_i^{\rm contact}+\mathbf f_i^{\rm actuator}.\]

A frictionless penalty contact approximation between a spherical projectile and a triangulated surface uses penetration \(\delta=\max(0,r_p-d(\mathbf x_p,\mathcal S))\) and

\[U_c=\tfrac12 k_c\delta^2,\qquad \mathbf F_c=-\nabla_{\mathbf x_p}U_c.\]

The surface receives equal-and-opposite force distributed to triangle vertices with barycentric weights. This is a numerical approximation; convergence depends on contact stiffness, timestep, and triangle quality. A penalty model does not establish no-penetration at finite stiffness. Total momentum conservation applies to the *closed system including ring supports and actuator reaction*.

## D. Wave propagation and transport

In a continuum elastic approximation, \(\rho\ddot{\mathbf u}=\nabla\cdot\boldsymbol\sigma+\mathbf b\), with \(\boldsymbol\sigma=\mathsf C:\boldsymbol\varepsilon\) and \(\boldsymbol\varepsilon=\tfrac12(\nabla\mathbf u+\nabla\mathbf u^T)\). For isotropic linear elasticity, longitudinal and shear speeds are

\[c_L=\sqrt{(\lambda+2\mu)/\rho},\qquad c_T=\sqrt{\mu/\rho}.\]

A discrete spring network need not reproduce those speeds until its effective elastic constants are calibrated. Mechanical energy flux in the continuum is \(\mathbf q_E=-\boldsymbol\sigma\cdot\dot{\mathbf u}\), and support energy-transfer measurements must use force times relative velocity across a transducer, not fixed-support force alone.

## E. Deployment and actuation

For a prescribed deployment coordinate \(q\), a generalized actuator obeys

\[M(q)\ddot q+C(q,\dot q)+\partial U/\partial q=\tau_{\rm act},\quad P_{\rm act}=\tau_{\rm act}\dot q.\]

Electrical energy demand obeys \(E_{\rm supply}\ge\int \max(P_{\rm act},0)/\eta_{\rm act}\,dt\), and actuator constraints require \(|\tau_{\rm act}|\le\tau_{\max}\), \(|P_{\rm act}|\le P_{\max}\). Prescribed trajectories that exceed either limit are *infeasible under those assumptions*. Earlier Stage 9/11 examples exceeded illustrative limits; they must not be advertised as viable deployment hardware.

## F. Piezoelectric conversion and electrical energy

A linear lumped piezoelectric coupling is

\[m\ddot x+c\dot x+kx+\theta V=F(t),\qquad C\dot V+V/R=\theta\dot x.\]

Multiplying by \(\dot x\) and \(V\) gives

\[\frac{d}{dt}\left(\tfrac12m\dot x^2+\tfrac12kx^2+\tfrac12CV^2\right)=F\dot x-c\dot x^2-V^2/R.\]

The reciprocal coupling terms cancel. The resistor dissipates \(\int V^2/R\,dt\): **resistive load energy is not battery-stored energy**. Real storage needs a rectifier/converter, finite state of charge, voltage and current constraints, and a coupled electrical load. Idealized post-processing of harvested power is only an upper-bound ledger.

## G. Energy and momentum audit

For a complete system, including sources, support, field and heat,

\[E_{\rm mech}(t)+E_{\rm elec}(t)+E_{\rm heat}(t)+E_{\rm other}(t)-E_{\rm total}(0)=W_{\rm external}(0,t).\]

Numerical residual is \(r_E=E_{\rm accounted}(t)-E_{\rm accounted}(0)-W_{\rm external}\). A small \(r_E\) demonstrates bookkeeping/integration consistency **only**, not empirical validity. Momentum balance includes the time integral of external support and actuator forces.

## H. Field-coupling hypothesis and its limitation

For a neutral linearly polarizable particle in a prescribed quasistatic field, \(U=-\tfrac12\alpha E^2\), and \(\mathbf F=\tfrac12\alpha\nabla E^2\). Positive \(\alpha\) attracts it to high-field regions. A static conservative field cannot permanently harvest its energy. No Maxwell-consistent field-source and microscopic coupling has been demonstrated for the proposed self-supporting barrier. The claim that energy 'fills interatomic gaps' should be interpreted as excitation energy in bonds, phonons, electromagnetic modes or transducer states, **not** energy forming a new solid surface in vacuum.

## I. What would be required to claim a buildable device

1. Select a material and geometry; determine density, modulus, damping, fracture and fatigue properties experimentally.
2. Calibrate a deployable linkage/membrane and show actuator forces, power, time and work within specifications.
3. Calibrate contact and wave dynamics against measured impacts, including rebound and permanent damage.
4. Measure transducer force–voltage characteristics, storage conversion and back-action under realistic loading.
5. Establish temporal and spatial convergence with constitutive parameters held fixed.
6. Compare flat, fixed-dome and deployable designs at matched mass, footprint, impact conditions and net energy consumption.

**Conclusion:** Classical mechanics and electromechanics permit *material-supported* energy-absorbing, energy-harvesting domes. The present simulations do not prove this particular deployable configuration is buildable, advantageous, or capable of forming a free-standing field barrier.
