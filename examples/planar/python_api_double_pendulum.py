#!/usr/bin/env python3
"""Hierarchical Sim2D model built with the Python API."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
try:
    from practical_mechanical_simulation import Sim2D
except ModuleNotFoundError:
    sys.path.insert(0, str(ROOT / "python" / "src"))
    from practical_mechanical_simulation import Sim2D


def pendulum_link(model, name, *, length, mass, center, color, width=0.08):
    body = Sim2D.rigid_body(
        model,
        name,
        mass=mass,
        inertia=mass * (length**2 + width**2) / 12,
        position=center,
    )
    inner = Sim2D.marker(body, "inner", position=[-length / 2, 0.0])
    outer = Sim2D.marker(body, "outer", position=[length / 2, 0.0])
    Sim2D.graphics(body, show_default=False, color=color)
    Sim2D.graphic(
        body,
        "member",
        shape="cylinder",
        markers=[inner, outer],
        radius=width / 2,
    )
    return {"body": body, "inner": inner, "outer": outer}


def double_pendulum(model, name, ground_pin):
    first = pendulum_link(
        model,
        f"{name}.first_link",
        length=1.0,
        mass=1.0,
        center=[0.5, 0.0],
        color="steelblue",
    )
    second = pendulum_link(
        model,
        f"{name}.second_link",
        length=0.8,
        mass=0.7,
        center=[1.4, 0.0],
        color="darkorange",
    )
    base_joint = Sim2D.revolute(
        model,
        f"{name}.base_joint",
        markers=[first["inner"], ground_pin],
        rotation_coordinates=True,
    )
    elbow_joint = Sim2D.revolute(
        model,
        f"{name}.elbow_joint",
        markers=[second["inner"], first["outer"]],
        rotation_coordinates=True,
    )
    Sim2D.gravity(
        model,
        f"{name}.gravity",
        acceleration=[0.0, -9.81],
        bodies=[first["body"], second["body"]],
    )
    return {"base_joint": base_joint, "elbow_joint": elbow_joint}


def build_model():
    model = Sim2D.Model(
        "hierarchical_double_pendulum",
        title="Hierarchical Sim2D Python API double pendulum",
        source_directory=ROOT,
    )
    Sim2D.analysis(model, mode="automatic")
    Sim2D.simulation(
        model,
        end_time=4.0,
        frames_per_second=60,
        relative_tolerance=1.0e-7,
        absolute_tolerance=1.0e-9,
        maximum_step=0.01,
    )
    Sim2D.graphics(model, background="white", body_palette="colorblind")

    ground = Sim2D.ground(model, "ground")
    ground_pin = Sim2D.marker(ground, "pin")
    pendulum = double_pendulum(model, "pendulum", ground_pin)
    Sim2D.state_selection(
        model,
        method="preferred",
        preferred_velocities=[
            Sim2D.variable(pendulum["base_joint"], "omega"),
            Sim2D.variable(pendulum["elbow_joint"], "omega"),
        ],
        allow_fallback=False,
    )
    return model


if __name__ == "__main__":
    output = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else ROOT / "results" / "examples" / "planar" / "python-api-double-pendulum.simp"
    )
    run = Sim2D.run(build_model(), output=output, overwrite=True)
    print(run.stdout, end="")
    print(f"Wrote {output}")
    print("Open it with bin/simpView")
