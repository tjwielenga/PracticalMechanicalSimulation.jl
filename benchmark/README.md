# Planar Performance Benchmarks

These scripts exercise the supported TOML model path without adding large
generated models to the source or paper test suites. Timing results are machine
and Julia-version dependent; they are baselines rather than acceptance tests.

## Spatial allocation benchmark

The spatial allocation benchmark measures two representative Sim3D models
after compilation and model loading. It requests only the initial and final
output samples so that the reported allocation principally measures equation
evaluation and integration rather than saved-result storage.

```bash
julia --project=. benchmark/spatial_allocation_benchmark.jl
```

The default reports the median-time run from three repetitions. A different
positive repetition count may be supplied on the command line. The table
includes cumulative allocated bytes, allocation count, garbage-collection
time, solver work, and a final-state norm. The last quantities help detect an
optimization that changes the numerical path instead of merely reducing
temporary storage. Allocated MiB is cumulative allocation during the run, not
the program's peak resident memory.

An initial comparison was made September 25, 2026 with Julia 1.12.6 on the
arm64 development Mac. The second measurement converted the fixed-size
spatial rotation, quaternion, marker, and rotation-Jacobian arithmetic to
static arrays.

| model | version | run s | allocated MiB | allocations |
|:--|:--|--:|--:|--:|
| Spatial four-bar | Ordinary arrays | 0.041540 | 160.4 | 3,194,463 |
| Spatial four-bar | Static spatial algebra | 0.009005 | 27.7 | 269,428 |
| Bristle rolling tire | Ordinary arrays | 0.536167 | 1,789.4 | 35,468,251 |
| Bristle rolling tire | Static spatial algebra | 0.063094 | 199.6 | 1,691,417 |

The four-bar became 4.6 times faster and allocated 83 percent fewer bytes.
The tire became 8.5 times faster and allocated 89 percent fewer bytes. Its
allocation count fell by about 95 percent. Accepted and rejected steps,
residual and Jacobian evaluations, Newton iterations, and final-state norms
were identical before and after the change.

## Pendulum chains

Run the default 10-, 25-, and 50-link benchmark:

```bash
julia --project=. benchmark/pendulum_chain_benchmark.jl
```

The script generates each TOML model in memory, warms the complete two-link
path before measuring, and reports model loading separately from the complete
run. The run timing includes a fresh load. It also reports canonical system
size, state count, sparse Jacobian nonzeros, solver work, and maximum revolute
joint position error.

Additional sizes can be listed on the command line:

```bash
julia --project=. benchmark/pendulum_chain_benchmark.jl 5 20 100
```

To save and view one generated chain:

```bash
julia --project=. benchmark/run_pendulum_chain.jl 10 \
    results/benchmarks/pendulum-chain-10.simp
bin/simpView
```

The saved result embeds the generated TOML, so it can also be extracted with
the ordinary model-extraction command.

## Initial generic-LU baseline

Measured September 4, 2026 with Julia 1.12.6 on the arm64 development Mac.
Each model ran for 0.25 seconds with 11 requested output samples. These are
single post-warm-up measurements; they are sufficient to expose scaling but
are not statistical timing claims.

| links | variables | states | Jacobian nnz | load s/MiB | run s/MiB | variables/wall s | accepted/rejected | residual evals | Jacobians/factorizations | Newton iterations | max order | max joint error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 110 | 10 | 498 | 0.001037/1.280 | 0.6013/67.24 | 182.9 | 95/2 | 235 | 97/97 | 137 | 5 | 5.427e-11 |
| 25 | 275 | 25 | 1278 | 0.001887/4.025 | 11.27/243.5 | 24.40 | 123/2 | 281 | 125/125 | 155 | 5 | 6.675e-12 |
| 50 | 550 | 50 | 2578 | 0.003930/11.09 | 112.2/812.1 | 4.902 | 147/2 | 313 | 149/149 | 163 | 5 | 2.033e-12 |

Canonical size and Jacobian nonzeros grow approximately linearly, and loading
through configuration correction and automatic QR selection remains small.
The solver work counts also grow modestly: from 10 to 50 links, accepted steps
increase by about 1.55 times and Newton iterations by about 1.19 times. Run
time, however, increases by about 187 times. The loss is therefore in the cost
of each implicit step rather than the number of steps required.

A sampling profile of the 25-link case places most execution samples in
`generic_lufact!`, reached from `lu!` on the sparse Newton matrix. Scalar sparse
indexing inside that generic in-place factorization is also prominent. The
first performance experiment therefore replaced `lu!` with Julia's ordinary
sparse `lu` dispatch, which uses UMFPACK for these `Float64` CSC matrices.

## Dedicated sparse-LU result

The same benchmark after that one-line solver change produced:

| links | variables | states | Jacobian nnz | load s/MiB | run s/MiB | variables/wall s | GC % | accepted/rejected | residual evals | Jacobians/factorizations | Newton iterations | max order | max joint error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 110 | 10 | 498 | 0.001038/1.264 | 0.01381/60.82 | 7965 | 0.00 | 95/2 | 235 | 97/97 | 137 | 5 | 5.427e-11 |
| 25 | 275 | 25 | 1278 | 0.001859/4.025 | 0.04187/186.5 | 6568 | 0.00 | 123/2 | 281 | 125/125 | 155 | 5 | 6.675e-12 |
| 50 | 550 | 50 | 2578 | 0.003926/11.09 | 0.1070/422.7 | 5140 | 10.13 | 147/2 | 313 | 149/149 | 163 | 5 | 2.033e-12 |

The solver work and errors are identical. Run-time speedups are approximately
44 times for 10 links, 269 times for 25 links, and 1,049 times for 50 links.
The optimized 50-link case is about 7.7 times slower than the 10-link case for
a fivefold increase in canonical size, instead of 187 times slower.

Allocation remains substantial even though it is no longer masking the sparse
factorization behavior. No collection occurred during the measured 10- and
25-link intervals; garbage collection occupied about 10 percent of the
50-link interval. Reusing residual, scaling, selection-map, and factorization
workspaces is therefore a reasonable later optimization, but it is now a
secondary cost rather than the primary algorithmic bottleneck.

## Reused symbolic factorization

UMFPACK separates analysis of the sparse pattern from factorization of its
changing numerical values. Retaining that symbolic result across DDASSL steps
gives the following measurements. The last entry in the LU column is the
number of symbolic analyses.

| links | variables | states | Jacobian nnz | load s/MiB | run s/MiB | variables/wall s | GC % | accepted/rejected/corrector failures | residual evals | Jacobians/numeric LU/symbolic LU | Newton iterations | max order | max joint error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 110 | 10 | 498 | 0.001076/1.280 | 0.01350/55.68 | 8148 | 0.00 | 95/2/0 | 235 | 97/97/1 | 137 | 5 | 5.427e-11 |
| 25 | 275 | 25 | 1278 | 0.001988/4.025 | 0.03890/169.2 | 7069 | 0.00 | 123/2/0 | 281 | 125/125/1 | 155 | 5 | 6.675e-12 |
| 50 | 550 | 50 | 2578 | 0.003857/11.09 | 0.09945/382.1 | 5530 | 7.64 | 147/2/0 | 313 | 149/149/1 | 163 | 5 | 2.033e-12 |

Each smooth run performs one symbolic analysis and reuses it for every later
numerical factorization. Relative to fresh UMFPACK analysis at every step,
50-link time falls by about 7 percent and allocation by about 10 percent. The
integrator discards the symbolic result after two consecutive corrector
failures. A changed state selection stays in the same integrator call, retains
the complete BDF history, and starts a new symbolic analysis for the replacement
state-equation rows.

## Reused numerical factorization

The level-scaled iteration matrix is designed to change slowly as the BDF step
size changes. Retaining its numerical factors turns later correctors into
modified-Newton iterations. The factors are refreshed when the scaled BDF
coefficient leaves the DDASSL ratio window, when a reused matrix cannot
converge, after five attempted steps, or after an event. A failed stale-matrix
corrector is restarted from its predictor with fresh numeric values before the
time step is rejected.

With the DDASSL correction multiplier applied to reused factors, the same
benchmark produced:

| links | variables | states | Jacobian nnz | load s/MiB | run s/MiB | variables/wall s | GC % | accepted/rejected/corrector failures | residual evals | Jacobians/numeric LU/symbolic LU | Newton iterations | max order | max joint error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 110 | 10 | 498 | 0.001059/1.280 | 0.009932/48.33 | 11075 | 0.00 | 95/2/0 | 327 | 20/20/1 | 229 | 5 | 5.416e-11 |
| 25 | 275 | 25 | 1278 | 0.001884/4.025 | 0.02871/146.2 | 9579 | 0.00 | 123/2/0 | 392 | 25/25/1 | 266 | 5 | 6.669e-12 |
| 50 | 550 | 50 | 2578 | 0.003953/11.09 | 0.08028/340.6 | 6851 | 12.98 | 147/2/0 | 454 | 30/30/1 | 304 | 5 | 2.033e-12 |

The 50-link run now needs 30 Jacobian evaluations and numerical factorizations
for 149 attempted steps. Relative to symbolic-only reuse, its measured time
falls by about 19 percent and allocation by about 11 percent. Newton
back-solves increase because an older matrix is approximate, but they are much
cheaper than sparse numerical factorization. The original step sequence and
constraint error are retained.

## Ten-second damped-chain experiment

The short benchmarks above are useful for measuring scaling, but 0.25 seconds
does not allow a horizontal pendulum chain to develop much motion. A second
experiment therefore ran the 10-link chain from rest for 10 seconds with 1,001
saved samples. Each link was 0.2 m long, had a mass of 1 kg, and used the
slender-link planar inertia $mL^2/12$. The links initially lay along the
positive horizontal axis under gravity. The integration settings retained the
$10^{-7}$ relative tolerance, $10^{-9}$ absolute tolerance, and 0.005 s maximum
step used by the generated benchmark model.

The undamped chain produced physically plausible but visually violent motion.
Potential energy was transferred into rapid local spinning, with a sampled
maximum body angular speed of about 92.2 rad/s and a maximum relative joint
speed of about 154.6 rad/s. Adding 0.001 N m s/rad of viscous damping to every
revolute joint reduced those sampled maxima to about 85.1 and 121.1 rad/s,
respectively, but the motion remained energetic. Increasing every joint
damping coefficient to 0.01 N m s/rad produced an informative but visually
well-behaved response. Its sampled maxima were about 54.9 rad/s for a body and
93.0 rad/s across a joint. This last value is a useful visualization default;
it is a modeling choice for the benchmark rather than a proposed universal
damping value.

### State-coordinate comparison

The 0.01 N m s/rad model was then run with three state choices:

1. Automatic scaled, pivoted-QR selection. For this configuration it selected
   six body angular velocities and four body vertical velocities.
2. Ten preferred body angular velocities, one for each link.
3. Ten preferred revolute-joint relative angular velocities. This case enabled
   optional rotation coordinates on every revolute joint, adding explicit
   relative angle, angular-velocity, and angular-acceleration variables and
   their defining equations.

All three choices were accepted without fallback and all produced similar
qualitative motion. The body-coordinate systems contained 120 canonical
variables, including the ten damper load variables. Adding explicit joint
rotation coordinates increased the canonical system to 150 variables. Each
formulation still had ten selected states.

The following values are means of three full runs in one Julia process after
warming all three paths. Their order was alternated to reduce ordering bias.
Run time includes model loading, which is approximately one millisecond. The
measurements were made September 4, 2026 with Julia 1.12.6 on the arm64
development Mac.

| state choice | run s | allocated MiB | accepted/rejected | residual evaluations | Jacobians/numeric LU | Newton iterations | corrector failures |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| automatic QR | 0.8548 | 3526.1 | 6289/39 | 24661 | 1434/1434 | 18332 | 0 |
| body angular velocities | 0.6995 | 2951.7 | 6496/19 | 20047 | 1303/1303 | 13531 | 0 |
| relative joint angular velocities | 0.9010 | 3665.8 | 7356/21 | 21934 | 1476/1476 | 14556 | 0 |

For this open serial chain, preferred body angular velocities were about 18
percent faster than the automatic selection, allocated about 16 percent less,
and required about 26 percent fewer Newton iterations. Explicit relative joint
coordinates were about 5 percent slower than automatic selection and allocated
about 4 percent more. That modest cost is reasonable because the relative
formulation enlarges the unreduced system by 25 percent and exposes coordinates
that are often more meaningful to a mechanism analyst.

The experiment supports keeping automatic QR as the general default. Preferred
body coordinates can improve performance when the analyst knows a natural
choice, while optional relative coordinates provide a useful physical
interface without requiring the entire program to use reduced coordinates.
Long, nonlinear chain trajectories are sensitive to small numerical
differences, so qualitative agreement and solver statistics are more useful
here than pointwise comparison late in the run.

### State selection and initialization cost

The state-selection and initialization stages were timed separately using
warmed median measurements. The body-coordinate constraint partial matrix was
$20\mathbin{\times}30$; the relative-coordinate version was
$30\mathbin{\times}40$. A "QR pair" includes the column-pivoted factorization
used to divide dependent and independent velocity candidates and the separate
row-pivoted factorization used to detect redundant constraint equations.

| state choice | QR pair per selection pass | complete selection per pass | full model load | implicit initialization | total before integration |
| --- | ---: | ---: | ---: | ---: | ---: |
| automatic QR | 7.7 microseconds | 33.2 microseconds | 0.820 ms | 0.217 ms | 1.04 ms |
| body angular velocities | 7.8 microseconds | 33.1 microseconds | 0.806 ms | 0.208 ms | 1.01 ms |
| relative joint angular velocities | 17.7 microseconds | 49.4 microseconds | 0.896 ms | 0.239 ms | 1.13 ms |

The loader uses two passes. It first selects states, then rebuilds the model
with all candidate state-equation pairs allocated and only the selected pairs
active. It evaluates the selection diagnostics again on the final system.
Consequently, total QR work during loading is about twice the per-pass value in
the table. Even so, state selection is negligible beside a 0.7--0.9 second
simulation.

The horizontal starting geometry already satisfied its position equations, so
configuration projection performed no correction; its residual check took
about 8--9 microseconds. All three formulations required one simultaneous
implicit Newton correction to establish consistent accelerations, reactions,
and other algebraic variables.

State selection deliberately uses dense linear algebra. Its matrix contains
only velocity-level constraint rows and candidate velocities, so it is much
smaller than the full unreduced mechanical system. It uses dense pivoted QR,
and preferred-state validation uses a dense SVD of the complementary dependent
block. Position correction likewise uses a dense rectangular Jacobian,
pivoted QR for independent-row detection, and a weighted dense minimum-motion
solve when a correction is needed.

In contrast, simultaneous implicit initialization assembles sparse Jacobian
blocks and solves the resulting square sparse system through Julia's sparse
backslash operation, which dispatches to UMFPACK for these Float64 CSC
matrices. It currently performs a fresh factorization for each initialization
Newton correction. That is immaterial for this already-consistent case, but it
may deserve attention for large models requiring many initialization
iterations. DDASSL itself retains the sparse symbolic factorization and reuses
numerical factors according to its modified-Newton refresh policy.

The practical conclusions are:

- Dense rank-revealing methods are appropriate for the small state-selection
  and configuration-projection matrices encountered so far.
- Sparse factorization is essential for the full implicit system and was the
  decisive improvement in the chain-size benchmarks.
- Automatic selection and initialization are not current performance
  bottlenecks.
- State choice can nevertheless change integration work and allocation, so
  representative preferred-coordinate cases should remain in performance
  studies.
- Allocation remains substantial--roughly 3 GiB for this 10-second, 10-link
  calculation--and is a better next optimization target than replacing the
  small dense QR factorizations.

## Closed-loop parallelogram chains

The next generated benchmark is a connected sequence of parallelogram
four-bar cells. Every cell adds one horizontal coupler and shares its two
ground-pivoted rockers with its neighbors. The complete chain therefore has
one degree of freedom regardless of its length. This gives a substantially
different sparse system from the open pendulum chain: it contains many ideal
closed loops, but its motion remains simple enough to check visually and
kinematically.

Run the default 10-, 25-, and 50-cell benchmark with:

```bash
julia --project=. benchmark/parallelogram_chain_benchmark.jl
```

Alternative cell counts may be supplied on the command line. A viewable
10-second result can be generated with:

```bash
julia --project=. benchmark/run_parallelogram_chain.jl \
    10 results/benchmarks/parallelogram-chain-10-rocker-omega-10s.simp 10 1001
bin/simpView
```

The initial geometry places identical 0.2 m rockers at an angle of 0.7 rad and
connects their tips with 0.2 m couplers. Every body has a mass of 1 kg and the
appropriate slender-link inertia. Gravity drives the motion, while a
0.01 N m s/rad rotational damper at every ground pivot prevents persistent
undamped oscillation. The angular velocity of the first rocker,
`rocker00.omega`, is the preferred state.

An initial measurement on September 4, 2026 used the same 0.25-second interval
and solver tolerances as the open-chain benchmarks:

| cells | variables | states | Jacobian nnz | load s/MiB | run s/MiB | variables/wall s | GC % | accepted/rejected/corrector failures | residual evals | Jacobians/numeric LU/symbolic LU | Newton iterations | max order | max joint error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 262 | 1 | 1254 | 0.002143/4.371 | 0.02806/128.6 | 9336 | 0.00 | 67/1/0 | 322 | 23/23/1 | 253 | 5 | 5.920e-10 |
| 25 | 637 | 1 | 3084 | 0.008333/16.08 | 0.08450/327.4 | 7539 | 12.09 | 66/1/0 | 326 | 22/22/1 | 258 | 5 | 6.638e-10 |
| 50 | 1262 | 1 | 6134 | 0.03006/48.38 | 0.1777/638.8 | 7101 | 11.98 | 66/1/0 | 326 | 22/22/1 | 258 | 5 | 6.887e-10 |

Variables and sparse Jacobian nonzeros again grow linearly. Because all cells
are kinematically coupled into the same motion, accepted steps, Jacobian
evaluations, and Newton iterations remain essentially unchanged with size.
The 50-cell model is about 4.8 times the canonical size of the 10-cell model
and takes about 6.3 times as long for the complete load and run.

This topology also gives a useful qualification to the dense-QR conclusion
above. Automatic state selection operates on a $62\mathbin{\times}63$ dense
constraint partial matrix for 10 cells and a
$302\mathbin{\times}303$ matrix for 50 cells. Model loading grows from about
2.1 to 30.1 ms, faster than canonical system size, although it remains small
in absolute terms. The closed-loop benchmark should therefore be retained as
a way to determine when a sparse rank-revealing method becomes worthwhile.

The first version allowed automatic QR to select `coupler01.V_y` as the single
physical state and retained that choice throughout the simulation. That
translational component becomes a poor coordinate in portions of the motion.
Selecting `rocker00.omega` instead gives a nonsingular physical coordinate
through an entire rocker rotation. The preferred choice was accepted without
fallback for every benchmark size. In the 10-cell, 10-second run it also
reduced rejected steps from 61 to 1, Newton iterations from 7,795 to 5,336,
allocation from 3,738 to 2,807 MiB, and warmed run time from 0.895 to 0.672 s.
Maximum sampled joint-position error fell from
$7.13\mathbin{\times}10^{-9}$ m to $3.40\mathbin{\times}10^{-9}$ m. This is a
useful example of why automatic state selection may need to be revisited when
a long trajectory moves far from the configuration where its QR analysis was
performed.

## Staggered bouncing-ball contact bank

The third generated benchmark exercises abrupt force-law transitions rather
than ideal joint constraints. It creates independent spherical bodies above a
common horizontal ground plane. Their initial heights correspond to free-fall
impact times distributed from 0.25 to 0.75 s, preventing all contacts from
changing stiffness simultaneously. Each 0.1 m-radius, 1 kg ball uses the
existing 10,000 N/m plane-contact law with a damping factor of 0.15 and falls
under gravity from rest.

Run the default 10-, 25-, and 50-ball cases with:

```bash
julia --project=. benchmark/bouncing_ball_bank_benchmark.jl
```

Generate and view the 10-ball result with:

```bash
julia --project=. benchmark/run_bouncing_ball_bank.jl \
    10 results/benchmarks/bouncing-ball-bank-10.simp 1.5 301
bin/simpView
```

The contact table adds root-function evaluations, located events, history
restarts, wall time per event, and maximum penetration observed at the 301
requested output samples. Rejected steps are reported separately from
corrector failures because an event-rich calculation may reject steps for
error control or transition localization without failing Newton iteration.

### Hard-restart baseline

The first implementation discarded all BDF history, returned to order one,
and reduced the next step by a factor of four at every located transition.

| balls | variables | Jacobian nnz | run s/MiB | variables/wall s | accepted/rejected/failures | root evaluations | events/restarts | ms/event | max sampled penetration, m |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 140 | 530 | 0.1973/933.0 | 709.5 | 2077/489/0 | 3529 | 45/45 | 4.385 | 0.05642 |
| 25 | 350 | 1325 | 1.024/5179 | 341.7 | 3717/1163/0 | 7478 | 114/114 | 8.984 | 0.05642 |
| 50 | 700 | 2650 | 3.572/20460 | 196.0 | 6033/2099/0 | 13409 | 225/225 | 15.87 | 0.05642 |

Canonical variables and Jacobian nonzeros are exactly linear in ball count,
and the Newton corrector succeeds throughout. Run time and allocation are not
linear, however. Increasing the model from 10 to 50 balls multiplies canonical
size by five, run time by about 18, and allocation by about 22. Every located
transition restarts the BDF history of the complete system, and every root
evaluation scans all active contacts. Consequently, time per event also rises
with model size even though the mechanical Jacobian is block diagonal.

The identical sampled maximum penetration at each size is expected because
the balls are mechanically independent and the same range of impact speeds is
present in every bank. The relatively large 56 mm value reflects the chosen
10,000 N/m benchmark stiffness and the highest free-fall speed; it is not a
recommended physical contact stiffness. Increasing stiffness would reduce
penetration while making the force transition more demanding numerically.

This baseline confirms that the event path is robust--225 transitions are
located and restarted without a corrector failure in the 50-ball case--but it
also makes allocation and global event-restart cost much more visible than the
smooth linkage benchmarks. Reusable residual, root, and event workspaces are
therefore stronger optimization candidates before attempting a coupled-contact
model such as balls striking a common bushing-supported platform.

### Soft restart at order two

The compliant contact force is continuous at entry, exit, and the damping
cutoff even though a force derivative changes abruptly. The revised internal
root policy therefore retains recent solution history, caps the BDF order at
two, refreshes the numerical iteration matrix, and reduces the next step by a
factor of two. The mechanical runner selects this response automatically from
the compliant contact semantics, so no TOML option is added. At the lower
level, the former Boolean restart keyword is replaced by one `root_restart`
policy accepting `:hard`, `:soft`, or `:none`; generic DDASSL retains `:hard`
as its conservative default.

| balls | variables | Jacobian nnz | run s/MiB | variables/wall s | accepted/rejected/failures | root evaluations | events/restarts | ms/event | max sampled penetration, m |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 140 | 530 | 0.1425/682.3 | 982.7 | 1592/227/0 | 3045 | 45/45 | 3.166 | 0.05642 |
| 25 | 350 | 1325 | 0.6728/3370 | 520.2 | 2448/495/0 | 6144 | 112/112 | 6.007 | 0.05642 |
| 50 | 700 | 2650 | 2.103/11850 | 332.8 | 3472/875/0 | 10845 | 224/224 | 9.389 | 0.05642 |

Relative to the hard restart, the soft policy reduces run time by about 28,
34, and 41 percent for 10, 25, and 50 balls. Allocation falls by about 27, 35,
and 42 percent, while rejected steps fall by 54--58 percent. Maximum sampled
penetration is unchanged and no corrector fails. The event totals differ
slightly in the two larger models because accumulated trajectory differences
move late rebound transitions relative to the fixed 1.5 s endpoint.

For an accuracy check, both 10-ball trajectories were compared with a run at
$10^{-10}$ relative tolerance, $10^{-12}$ absolute tolerance, and a 0.0001 s
maximum step. The hard and soft results had nearly identical worst-case
differences: approximately 0.109 mm in position and 0.0102 m/s in velocity.
The soft result's RMS differences were somewhat larger--0.0205 versus
0.0126 mm in position and 0.000467 versus 0.000329 m/s in velocity. Thus the
soft restart gives a substantial performance improvement without changing the
worst observed error, but it is not accuracy-neutral in every norm.

This result supports a second-order restart for continuous compliant forces.
An event that applies an impulse or otherwise changes a state discontinuously
would still require a true order-one restart and explicit post-event state
initialization, which the current contact model does not perform.

## Torsional rotor trains

The fourth generated benchmark is a smooth stiff system rather than a
mechanism with changing geometry or contact events. It consists of identical
fixed-center rotors, each supported by a revolute joint to ground. Neighboring
rotors are coupled by torsional spring-dampers, and the first rotor is coupled
to ground by one more identical element. The last rotor starts 0.25 rad from
equilibrium and all other positions and velocities start at zero.

Run the default 10-, 25-, and 50-rotor cases with:

```bash
julia --project=. benchmark/rotor_train_benchmark.jl
```

Generate and view the 10-rotor result with:

```bash
julia --project=. benchmark/run_rotor_train.jl \
    10 results/benchmarks/rotor-train-10.simp 1.0 201
bin/simpView
```

Each rotor has a mass of 1 kg, a pitch radius of 0.1 m, and a polar inertia of
0.005 kg m$^2$. Every coupling has a stiffness of 20 N m/rad and a damping
time scale of 0.001 s, so its damping coefficient is 0.02 N m s/rad. The body
angular velocities are specified as the preferred states. The benchmark runs
for 1 s with a 0.002 s maximum step and compares all sampled rotor angles and
angular velocities with an independent modal solution of

$$
M\ddot{\theta} + \tau K\dot{\theta} + K\theta = 0.
$$

The following warmed results were measured on September 4, 2026:

| rotors | variables | states | Jacobian nnz | run s/MiB | variables/wall s | accepted/rejected/failures | Jacobians/numeric/symbolic | max angle error, rad | max omega error, rad/s | frequency range, rad/s | final energy ratio |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 120 | 10 | 367 | 0.06069/291.5 | 1977 | 908/11/0 | 185/185/1 | 5.989e-7 | 6.147e-5 | 9.453--125.1 | 0.03332 |
| 25 | 300 | 25 | 922 | 0.1529/673.0 | 1962 | 867/1/0 | 175/175/1 | 5.942e-7 | 6.782e-5 | 3.895--126.3 | 0.03356 |
| 50 | 600 | 50 | 1847 | 0.2771/1363 | 2165 | 827/0/0 | 166/166/1 | 7.823e-7 | 8.797e-5 | 1.967--126.4 | 0.03352 |

Canonical variables, selected states, Jacobian nonzeros, run time, and
allocation all scale nearly linearly. Increasing the train from 10 to 50
rotors multiplies the variable count by five, run time by 4.6, and allocation
by 4.7; the variables-per-wall-second measure remains close to 2,000. This is
the clearest linear-scaling result among the benchmarks so far.

The highest natural frequency approaches about 126.5 rad/s as rotors are
added, while the lowest frequency falls because a longer train has a softer
collective mode. Consequently, the modal frequency range widens from about
13:1 to 64:1. The integrator handles that increasing stiffness ratio without
a corrector failure, and one symbolic factorization serves each complete run.
Numerical factors are reused across several accepted steps rather than being
rebuilt for every Newton correction.

The maximum errors remain essentially independent of train size and are
consistent with the requested tolerances. Because the damping matrix is
proportional to the stiffness matrix, the modal reference is exact apart from
floating-point evaluation and gives this benchmark a substantially stronger
accuracy check than a constraint-error measurement alone. The nearly equal
final energy ratios also show that the three sizes reproduce the intended
physical damping consistently.
