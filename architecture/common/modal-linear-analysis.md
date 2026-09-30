# Sparse Modal Linear Analysis

Status: implemented for planar and spatial models

## 1. Operating point

Modal analysis linearizes the active implicit equations at one instant. The
operating point must satisfy the model's position and velocity constraints, but
it need not be a static equilibrium. The ordinary simultaneous initializer
solves the acceleration, reaction, and applied-force variables at the specified
configuration and velocity before linearization. Consequently, a model can be
linearized at its declared initial state or at a configuration transferred from
a saved result.

With `initialization = "static_equilibrium"`, static equilibrium is solved
first and all physical velocity and acceleration values are set to zero before
the dynamic consistency solve. This produces the conventional modes about a
stationary equilibrium.

## 2. Complete implicit linearization

Let the active canonical variables be $y$, including acceleration, velocity,
position, reaction, applied-geometry, applied-rate, and applied-load variables.
The assembled equations are

$$
F(t,y,\dot y)=0.
$$

At the operating point, their perturbation is

$$
J\,\delta y+E\,\delta\dot y=0,
\qquad
J=\frac{\partial F}{\partial y},
\qquad
E=\frac{\partial F}{\partial\dot y}.
$$

The ordinary analytical sparse Jacobian supplies both matrices. Evaluation
with derivative coefficient zero gives $J$; subtracting that result from an
evaluation with coefficient one gives $E$. No force element implements a
second, modal-specific equation set.

For a spatial body, the operating-point Euler parameters still evaluate its
finite orientation. Their four perturbation variables and the three
orientation-bridge equations plus normalization equation are removed from the
modal pencil. The body-fixed pseudo angles provide the three independent local
rotation columns instead. Mechanical orientation partials are therefore the
same local partials used by spatial dynamics. After solving a mode, the viewer
receives the equivalent tangent Euler-parameter perturbation

$$
\delta p=\frac{1}{2}Q(p)\delta\vartheta^b.
$$

This reconstructed value is for display; Euler parameters are not independent
modal coordinates.

The spatial regression examples include a one-coordinate torsional pendulum
with a closed-form frequency and a three-link pendulum checked against an
independently assembled $3\times3$ mass, damping, and stiffness model. Both
retain the complete sparse mechanism equations during the modal solve. A
spring-driven spatial four-bar also verifies that complete redundant
constraint families remain removed during linearization.

For $\delta y=\hat y e^{st}$,

$$
(J+sE)\hat y=0.
$$

Algebraic definitions remain in this descriptor pencil. Their modal values are
therefore solved together with body motion and reactions. Spring extensions,
force magnitudes, marker geometry, and other algebraic quantities appear in a
mode shape but introduce no finite modes because their columns in $E$ are zero.

## 3. Sparse shift-invert solution

For shift $\sigma$, define

$$
U=(J+\sigma E)^{-1}E.
$$

If $\theta$ is a nonzero eigenvalue of $U$, the corresponding physical
eigenvalue is

$$
s=\sigma-\frac{1}{\theta}.
$$

Only selected position and velocity variables have nonzero columns in $E$.
Let $c$ contain those column indices and let

$$
U_c=(J+\sigma E)^{-1}E_{:,c}.
$$

The nonzero eigenvalues of the full shift-invert operator are exactly the
eigenvalues of the small matrix $U_c[c,:]$. The implementation therefore:

1. assembles sparse $J$ and $E$ once;
2. factors the sparse shifted matrix once with UMFPACK;
3. solves simultaneously for the $2n_s$ nonzero columns of $E$, where $n_s$
   is the number of selected mechanical states; and
4. solves a dense eigenproblem of order $2n_s$, not of canonical-system order.

The full canonical eigenvector is recovered by a matrix multiplication and one
division by $\theta$. This retains all algebraic and reaction components while
avoiding a dense physical-coordinate nullspace. It is efficient for the
present model sizes and calculates the complete finite spectrum. A later very
large model can use a Krylov method on the same sparse shift-invert operator to
request only part of that spectrum.

Each reported mode is phase adjusted and normalized so its largest scaled
physical displacement has unit magnitude. The solver records the relative
error in the original equation

$$
J\hat y+sE\hat y=0.
$$

## 4. Interpretation

For $s=\alpha+i\omega_d$, the stored values are

$$
f_n=\frac{|s|}{2\pi},
\qquad
f_d=\frac{|\omega_d|}{2\pi},
\qquad
\zeta=-\frac{\alpha}{|s|}.
$$

At static equilibrium these have the usual small-vibration interpretation. At
a moving or accelerating point they are instantaneous frozen-time properties.
The present implementation includes velocity-dependent force derivatives and
all operating-point reactions through the complete implicit Jacobian. A
nonsmooth force transition, such as the precise instant of contact engagement,
does not have a unique linearization and should not be used as an operating
point.

## 5. Result storage and viewing

A modal `.simp` file stores one canonical operating point and one complex
canonical vector per requested mode. Real and imaginary parts are separate HDF5
datasets. Frequencies, damping ratios, equation errors, and the shift are also
stored.

The viewer creates display samples only when the file is opened. It applies a
phase rotation to the selected complex mode and adds the scaled perturbation to
the operating point, so animation does not enlarge the stored result. The
viewer provides a mode selector; the optional second command-line argument
chooses its initial mode number.
