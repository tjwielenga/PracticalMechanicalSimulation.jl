"""Python model-building interface for the planar Sim2D solver."""

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
        super().__init__(name, dimension="planar", **properties)


_ELEMENTS = {
    "ground": "ground",
    "rigid_body": "rigid_body",
    "flexible_beam": "flexible_beam",
    "floating_marker": "floating_marker",
    "revolute": "revolute",
    "inplane": "inplane",
    "perp": "perp",
    "translational": "translational",
    "fixed": "fixed",
    "span": "span",
    "distance_coordinate": "distance_coordinate",
    "coupler": "coupler",
    "gear_pair": "gear_pair",
    "rack_and_pinion": "rack_and_pinion",
    "pulley": "pulley",
    "belt": "belt",
    "belt_span": "belt_span",
    "gravity": "gravity",
    "applied_force": "applied_force",
    "applied_torque": "applied_torque",
    "torsional_spring_damper": "torsional_spring_damper",
    "spanning_force": "spanning_force",
    "bushing": "bushing",
    "curve": "curve",
    "plane_contact": "plane_contact",
    "curve_contact": "curve_contact",
    "flat_follower_contact": "flat_follower_contact",
    "surface_friction": "surface_friction",
    "rotational_motion": "rotational_motion",
    "translational_motion": "translational_motion",
    "revolute_friction": "revolute_friction",
    "translational_friction": "translational_friction",
    "inplane_friction": "inplane_friction",
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
