# Planar friction examples

Four models exercise the scalar carried-shear friction law:

- [`sliding-block-surface-friction.toml`](../../models/planar/sliding-block-surface-friction.toml)
  adds one-sided tangential friction to a compliant sphere-plane contact.

- [`revolute-bearing-friction.toml`](../../models/planar/revolute-bearing-friction.toml)
  lets a rotor coast to rest, stick, and break away under a later torque.
- [`translational-guide-friction.toml`](../../models/planar/translational-guide-friction.toml)
  does the same for a slider in a compound translational joint.
- [`inplane-friction-block.toml`](../../models/planar/inplane-friction-block.toml)
  applies tangential friction to a standalone bilateral inplane primitive.

Run any model with `simp2d`, for example:

```bash
./bin/simp2d models/planar/sliding-block-surface-friction.toml \
    --output results/examples/planar/sliding-block-surface-friction.simp \
    --overwrite
bin/simpview-web
```

The shear state preserves the sticking force at zero velocity. Each example
starts with slip, remains nearly stationary after stopping, and later exceeds
the static capacity under a ramped load. SimpView displays the equal-and-
opposite force or torque arrows and exposes shear, slip, normal load, and
friction magnitude for plotting.
