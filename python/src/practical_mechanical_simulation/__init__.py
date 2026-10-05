"""Python model builders for Practical Mechanical Simulation.

The capitalized module aliases intentionally match the Julia ``Sim2D`` and
``Sim3D`` namespaces::

    from practical_mechanical_simulation import Sim2D, Sim3D
"""

from . import simp2d as Sim2D
from . import simp3d as Sim3D
from ._builder import ElementRef, SolverError, SolverRun, VariableRef

__all__ = [
    "Sim2D",
    "Sim3D",
    "ElementRef",
    "VariableRef",
    "SolverRun",
    "SolverError",
]
