#!/usr/bin/env python3
"""Hierarchical Sim3D model built with the Python API."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
try:
    from practical_mechanical_simulation import Sim3D
except ModuleNotFoundError:
    sys.path.insert(0, str(ROOT / "python" / "src"))
    from practical_mechanical_simulation import Sim3D


def pendulum_link(model, name, *, length, mass, center, color, width=0.08):
    transverse_inertia = mass * (length**2 + width**2) / 12
    axial_inertia = mass * width**2 / 6
    body = Sim3D.rigid_body(
        model,
        name,
        mass=mass,
        inertia=[axial_inertia, transverse_inertia, transverse_inertia],
        position=center,
    )
    Sim3D.graphics(
        body,
        shape="box",
        size=[length, width, width],
        color=color,
    )
    inner = Sim3D.marker(body, "inner", position=[-length / 2, 0.0, 0.0])
    outer = Sim3D.marker(body, "outer", position=[length / 2, 0.0, 0.0])
    return {"body": body, "inner": inner, "outer": outer}


def double_pendulum(model, name, ground_pin):
    first = pendulum_link(
        model,
        f"{name}.first_link",
        length=1.0,
        mass=1.0,
        center=[0.5, 0.0, 0.0],
        color="steelblue",
    )
    second = pendulum_link(
        model,
        f"{name}.second_link",
        length=0.8,
        mass=0.7,
        center=[1.4, 0.0, 0.0],
        color="darkorange",
    )
    base_joint = Sim3D.revolute(
        model,
        f"{name}.base_joint",
        markers=[first["inner"], ground_pin],
        rotation_coordinates=True,
    )
    elbow_joint = Sim3D.revolute(
        model,
        f"{name}.elbow_joint",
        markers=[second["inner"], first["outer"]],
        rotation_coordinates=True,
    )
    Sim3D.gravity(
        model,
        f"{name}.gravity",
        acceleration=[0.0, -9.81, 0.0],
        bodies=[first["body"], second["body"]],
    )
    return {"base_joint": base_joint, "elbow_joint": elbow_joint}


def build_model():
    model = Sim3D.Model(
        "hierarchical_double_pendulum",
        title="Hierarchical Sim3D Python API double pendulum",
        source_directory=ROOT,
    )
    Sim3D.analysis(model, mode="automatic")
    Sim3D.simulation(
        model,
        end_time=4.0,
        frames_per_second=60,
        relative_tolerance=1.0e-7,
        absolute_tolerance=1.0e-9,
        maximum_step=0.01,
    )
    Sim3D.graphics(model, background="white", body_palette="colorblind")

    ground = Sim3D.ground(model, "ground")
    ground_pin = Sim3D.marker(ground, "pin")
    Sim3D.graphics(
        ground_pin,
        shape="xy_frame",
        axis_length=0.35,
        plane_size=0.16,
        plane_color="gray65",
        label="ground.pin",
    )
    pendulum = double_pendulum(model, "pendulum", ground_pin)
    Sim3D.state_selection(
        model,
        method="preferred",
        preferred_velocities=[
            Sim3D.variable(pendulum["base_joint"], "omega"),
            Sim3D.variable(pendulum["elbow_joint"], "omega"),
        ],
        allow_fallback=False,
    )
    return model


if __name__ == "__main__":
    output = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else ROOT / "results" / "examples" / "spatial" / "python-api-double-pendulum.simp"
    )
    run = Sim3D.run(build_model(), output=output, overwrite=True)
    print(run.stdout, end="")
    print(f"Wrote {output}")
    print("Open it with bin/simpView")
