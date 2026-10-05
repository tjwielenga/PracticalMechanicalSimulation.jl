"""Python model-building interface for the spatial Sim3D solver."""

from __future__ import annotations

from typing import Any

from ._builder import (
    ElementRef,
    SolverError,
    SolverRun,
    VariableRef,
    analysis,
    beam_marker,
    document,
    element,
    graphic,
    graphics,
    initial_conditions,
    marker,
    parameters,
    run_model,
    set_properties,
    simulation,
    state_selection,
    variable,
    write_model,
)
from ._builder import Model as _Model


class Model(_Model):
    def __init__(self, name: object, **properties: Any) -> None:
        super().__init__(name, dimension="spatial", **properties)


_ELEMENTS = {
    "ground": "ground",
    "rigid_body": "rigid_body",
    "flexible_beam": "flexible_beam",
    "gravity": "gravity",
    "applied_force": "applied_force",
    "directed_torque": "directed_torque",
    "applied_torque": "applied_torque",
    "spanning_force": "spanning_force",
    "spherical": "spherical",
    "perp": "perp",
    "cv_phase": "cv_phase",
    "inplane": "inplane",
    "inline": "inline",
    "hinge": "hinge",
    "orient": "orient",
    "revolute": "revolute",
    "fixed": "fixed",
    "constant_velocity": "constant_velocity",
    "cylindrical": "cylindrical",
    "translational": "translational",
    "rotational_motion": "rotational_motion",
    "translational_motion": "translational_motion",
    "spanning_motion": "spanning_motion",
    "span": "span",
    "directed_distance": "directed_distance",
    "bushing": "bushing",
    "plane_contact": "plane_contact",
    "curve": "curve",
    "curve_contact": "curve_contact",
    "flat_follower_contact": "flat_follower_contact",
    "surface_friction": "surface_friction",
    "revolute_friction": "revolute_friction",
    "translational_friction": "translational_friction",
    "inplane_friction": "inplane_friction",
    "rolling_tire": "rolling_tire",
    "coupler": "coupler",
    "gear_pair": "gear_pair",
    "rack_and_pinion": "rack_and_pinion",
    "pulley": "pulley",
    "belt": "belt",
    "belt_span": "belt_span",
    "equation_component": "equation_component",
}


def _builder(kind: str):
    def build(model: Model, name: object, **properties: Any) -> ElementRef:
        return element(model, kind, name, **properties)

    return build


for _name, _kind in _ELEMENTS.items():
    globals()[_name] = _builder(_kind)
    globals()[_name].__name__ = _name


def run(model: Model, **options: Any) -> SolverRun:
    return run_model(model, **options)


simulate = run

__all__ = [
    "Model",
    "ElementRef",
    "VariableRef",
    "SolverRun",
    "SolverError",
    "document",
    "element",
    "marker",
    "graphic",
    "graphics",
    "set_properties",
    "analysis",
    "simulation",
    "state_selection",
    "initial_conditions",
    "parameters",
    "variable",
    "write_model",
    "beam_marker",
    "run",
    "simulate",
    *_ELEMENTS,
]
