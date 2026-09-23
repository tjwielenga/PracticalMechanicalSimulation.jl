# Tire bristle lift-off experiment

This experiment exercises the optional lumped-bristle tire on a single wheel.
It does not change the Large Van model. The input is
[`bristle-tire-liftoff.toml`](../../models/spatial/bristle-tire-liftoff.toml).
Its patch-length and stiffness curves are illustrative, not measured tire
data.

The wheel first settles on a road tilted ten degrees sideways. A torsional
spring holds the wheel spin, while the tire's lateral bristle shear balances
the downhill component of gravity. After static relaxation, a smooth upward
force pulse acts for 0.12 s. The wheel leaves the road, falls back, and
bounces again. This checks both loss and re-establishment of contact without
changing tire models between static and dynamic analysis.

At the dynamic handoff, normal force is about 96.6 N and lateral force about
17.0 N, consistent with a 10 kg assembly on a ten-degree slope. Normal and
tangential tire forces become exactly zero during flight. By 0.32 s, the
stored lateral shear is effectively zero. The first sampled re-contact is
near 0.45 s; its lateral force develops from new slip, not from pre-lift-off
shear. A second touchdown occurs near 0.67 s. The run produces 61 output
frames over one second and completes without a tire event restart.

The main numerical finding was at the zero-load boundary. With normal load
and both tangential trial forces equal to zero, automatic differentiation of
the friction-ellipse norm had produced non-finite Jacobian entries. The
unloaded branch now uses its zero derivative, while retaining the loaded
ellipse and the exact zero-force law. A separate smooth shear-rate blend was
tested and removed because it would slowly dissipate shear in a parked,
loaded tire. The result still has a constitutive corner at lift-off, but the
implicit integrator crosses it in this experiment.

Run and inspect the example with:

```bash
julia --project=. bin/simp3d models/spatial/bristle-tire-liftoff.toml \
  --output results/examples/spatial/bristle-tire-liftoff.simp --overwrite
bin/simpView
```

The spatial tire-bristle regression checks also cover the zero-load sparse
Jacobian, the friction-ellipse limit at every saved sample, and shear loss
when a loaded contact patch shrinks. This is a single-wheel numerical test,
not a tire validation. We still need measured curves and wider speed/load
tests before using the model in a vehicle.

## Fore–aft companion with a braked wheel

The [fore–aft model](../../models/spatial/bristle-tire-fore-aft-liftoff.toml)
tilts the road about its $y$ axis instead. Its revolute wheel has an applied
torsional spring-damper acting as a brake, rather than a prescribed rotation.
At the static handoff, the longitudinal tire force is $-17.03$ N, lateral
force is essentially zero, and the brake applies $-8.41$ N·m. The wheel
rotates about $0.0042$ rad against the spring. The brake torque agrees with
the tire-force moment about the axle to within the static solve tolerance.

The same upward pulse lifts this wheel clear of the road. The first sampled
re-contact is at 0.45 s; stored longitudinal shear just before it is
effectively zero. The largest sampled longitudinal force is about 114 N,
well below the load-dependent friction limit throughout the run. Dynamic
integration completes with 404 accepted and 85 rejected steps, without a
forced event restart.

```bash
julia --project=. bin/simp3d models/spatial/bristle-tire-fore-aft-liftoff.toml \
  --output results/examples/spatial/bristle-tire-fore-aft-liftoff.simp --overwrite
bin/simpView
```
