# Rotating Flexible Blade

The rotating-blade example checks whether a mechanism assembled from the
floating-reference beam elements develops spin stiffening. Each individual
beam has constant small-deformation mass and stiffness matrices. It does not
contain an added local geometric-stiffness matrix. The complete assembly can
still stiffen because rotation develops axial tension in the connected beams,
and the finite-position mechanism equations couple that tension to transverse
motion.

The example is
[`rotating_flexible_blade.jl`](../../examples/spatial/rotating_flexible_blade.jl).
It builds a four-meter blade from rectangular beam segments, joins adjacent
end markers with fixed joints, and attaches the root to ground with a
revolute. Gravity first bends the stationary blade downward. A smooth motion
generator then raises the rotor speed from zero to its final value over four
seconds. The beam damping time scale is 0.03 s. This damping removes
high-frequency transients without hiding the change in static deflection or
modal frequency.

The four-segment sweep prefers the natural body-local elastic rates and uses
four free-end floating-frame velocities to complete the state set. The
eight-segment convergence case leaves state selection automatic. It monitors
the reciprocal condition estimate already produced by the sparse iteration-
matrix factorization and performs a new constraint QR before the partition
becomes singular. Each change retains the complete BDF history.

Run the complete sweep from the repository root:

```bash
julia --project=. examples/spatial/rotating_flexible_blade.jl
```

The run writes
`results/examples/spatial/rotating-flexible-blade.simp` and
`results/examples/spatial/rotating-flexible-blade-eight-segment.simp`. The
first contains the four-segment blade and the second retains the smoother
geometric representation of the eight-segment comparison. Open either result
with `bin/simpView`. A single diagnostic case can also be run with

```bash
julia --project=. examples/spatial/rotating_flexible_blade.jl \
    --single 8 7.5 results/scratch/rotating-flexible-blade.simp
```

## Results

The speed sweep uses four beam segments. After each dynamic spin-up, modal
analysis linearizes the complete constrained system at the final moving
operating point.

| Rotor speed (rad/s) | Tip height after static (m) | Rotating tip height (m) | First bending frequency (Hz) |
|---:|---:|---:|---:|
| 0.0 | -0.465552 | -0.465552 | 0.910892 |
| 2.5 | -0.465552 | -0.380371 | 1.004349 |
| 5.0 | -0.465552 | -0.244883 | 1.244977 |
| 7.5 | -0.465552 | -0.152992 | 1.563981 |

The blade rises as centrifugal tension increases. The first bending frequency
rises by about 72 percent between rest and 7.5 rad/s. These two independent
changes show the same stiffening effect.

The stored result uses the smooth four-segment history. At 7.5 rad/s, an
eight-segment comparison gives a tip height of -0.151087 m and a
first bending frequency of 1.574673 Hz. The frequencies differ by less than
one percent, and the final tip deflections differ by about 1.3 percent.
This agreement is sufficient for the present demonstration without making the
routine speed sweep unnecessarily expensive.

The complete four-segment run required 699 accepted and 4 rejected steps. The
eight-segment run required 710 accepted and 24 rejected steps. It made 19
early state-partition changes without a corrector failure or history restart.
An experimental selector based directly on the modal response also found
nonsingular partitions, but it excited small high-frequency elastic responses
and required several times as many steps. Modal information is therefore best
kept as a diagnostic here; the inexpensive constraint QR and physically
preferred local beam coordinates performed better during integration.

The modal calculation is a frozen-operating-point analysis. In a rotating
system, forward and backward branches and strongly damped roots may also be
present. The reported frequency in this experiment is the lowest oscillatory
out-of-plane bending mode, not simply an assertion that every first-listed
root has the same physical identity at every speed.

## Scope

This experiment does not turn the individual beam into a general nonlinear
finite-element formulation. It does not add deformation-dependent inertia or
local stress stiffness inside a beam. The stiffening is produced by the
finite-motion assembly of several floating-reference members and the axial
loads transmitted through their fixed joints. More segments improve the
geometric representation at additional computational cost.
