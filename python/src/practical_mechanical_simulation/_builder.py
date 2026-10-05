"""Shared hierarchical builder and solver bridge for Sim2D and Sim3D."""

from __future__ import annotations

import copy
import math
import os
import shutil
import subprocess
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from numbers import Integral, Real
from pathlib import Path
from typing import Any

from . import _toml


_RESERVED_TABLES = {
    "model",
    "parameters",
    "analysis",
    "simulation",
    "state_selection",
    "initial_conditions",
    "graphics",
}


def _model_name(name: object) -> str:
    text = str(name)
    if not text:
        raise ValueError("model and element names cannot be empty")
    if any(not segment for segment in text.split(".")):
        raise ValueError("model and element names cannot contain empty path segments")
    return text


@dataclass(frozen=True)
class ElementRef:
    """Reference to a named model element."""

    model: "Model"
    name: str
    kind: str

    def __str__(self) -> str:
        return self.name


@dataclass(frozen=True)
class VariableRef:
    """Reference to a named variable belonging to a model element."""

    model: "Model"
    name: str

    def __str__(self) -> str:
        return self.name


@dataclass(frozen=True)
class SolverRun:
    """Completed invocation of an existing Sim2D or Sim3D solver."""

    dimension: str
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    output_path: Path | None


class SolverError(RuntimeError):
    """Raised when an invoked Sim2D or Sim3D process fails."""

    def __init__(self, run: SolverRun):
        details = run.stderr.strip() or run.stdout.strip() or "no solver output"
        super().__init__(
            f"{run.dimension} solver exited with status {run.returncode}: {details}"
        )
        self.run = run


def _normalized_number(value: Real) -> int | float:
    if isinstance(value, Integral):
        return int(value)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("generated model values must be finite")
    return number


def _normalize(model: "Model | None", value: Any) -> Any:
    if isinstance(value, (ElementRef, VariableRef)):
        if model is not None and value.model is not model:
            api_name = "Sim2D" if model.dimension == "planar" else "Sim3D"
            raise ValueError(
                f"a {api_name} reference belongs to a different model"
            )
        return value.name
    if isinstance(value, Enum):
        return _normalize(model, value.value)
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, Real):
        return _normalized_number(value)
    if isinstance(value, Mapping):
        return {str(key): _normalize(model, item) for key, item in value.items()}
    if hasattr(value, "tolist") and callable(value.tolist):
        return _normalize(model, value.tolist())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        items = [_normalize(model, item) for item in value]
        if items and all(isinstance(item, Real) and not isinstance(item, bool) for item in items):
            if any(isinstance(item, float) for item in items):
                return [float(item) for item in items]
        return items
    raise TypeError(f"unsupported model value of type {type(value).__name__}")


class Model:
    """Mutable hierarchical model description shared by the Python APIs."""

    def __init__(
        self,
        name: object,
        *,
        dimension: str,
        title: str | None = None,
        source_directory: str | os.PathLike[str] | None = None,
        **properties: Any,
    ) -> None:
        if dimension not in {"planar", "spatial"}:
            raise ValueError("model dimension must be 'planar' or 'spatial'")
        name_text = _model_name(name)
        model_table: dict[str, Any] = {
            "name": name_text,
            "title": name_text if title is None else str(title),
            "dimension": dimension,
        }
        for key, value in properties.items():
            if value is not None:
                model_table[str(key)] = _normalize(None, value)
        self.dimension = dimension
        self.source_directory = Path(source_directory or Path.cwd()).resolve()
        self._document: dict[str, Any] = {"model": model_table}

    def __repr__(self) -> str:
        name = self._document["model"]["name"]
        return f"Sim{'2D' if self.dimension == 'planar' else '3D'}.Model({name!r})"


def _table_at(
    model: Model, path: Sequence[str], *, create: bool = True
) -> dict[str, Any] | None:
    table = model._document
    for segment in path:
        child = table.get(segment)
        if child is None and create:
            child = {}
            table[segment] = child
        if child is None:
            return None
        if not isinstance(child, dict):
            raise ValueError(f"'{segment}' is not a model table")
        table = child
    return table


def _update(model: Model, path: Sequence[str], properties: Mapping[str, Any]) -> Model:
    table = _table_at(model, path)
    assert table is not None
    for key, value in properties.items():
        if value is not None:
            table[str(key)] = _normalize(model, value)
    return model


def document(model: Model) -> dict[str, Any]:
    """Return a detached copy of the generated model document."""

    return copy.deepcopy(model._document)


def element(model: Model, kind: object, name: object, **properties: Any) -> ElementRef:
    """Add any element supported by the dimension-specific TOML loader."""

    qualified_name = _model_name(name)
    segments = qualified_name.split(".")
    if segments[0] in _RESERVED_TABLES:
        raise ValueError(f"'{qualified_name}' conflicts with a reserved model table")
    parent = _table_at(model, segments[:-1])
    assert parent is not None
    if segments[-1] in parent:
        raise ValueError(f"model element '{qualified_name}' is already defined")
    table: dict[str, Any] = {"type": str(kind)}
    for key, value in properties.items():
        if value is not None:
            table[str(key)] = _normalize(model, value)
    parent[segments[-1]] = table
    return ElementRef(model, qualified_name, str(kind))


def marker(owner: ElementRef, name: object, **properties: Any) -> ElementRef:
    """Add a body/ground marker or reference a generated beam marker."""

    if owner.kind not in {"ground", "rigid_body", "flexible_beam"}:
        raise ValueError(
            "markers must be owned by ground, a rigid body, or a flexible beam"
        )
    local_name = _model_name(name)
    if "." in local_name:
        raise ValueError("a marker's local name cannot contain a period")
    if owner.kind == "flexible_beam":
        if properties:
            raise ValueError(
                "generated flexible-beam markers do not accept marker properties"
            )
        return beam_marker(owner, local_name)
    return element(owner.model, "marker", f"{owner.name}.{local_name}", **properties)


def variable(owner: ElementRef, name: object) -> VariableRef:
    """Return a qualified variable reference belonging to an element."""

    local_name = _model_name(name)
    if "." in local_name:
        raise ValueError("a variable's local name cannot contain a period")
    return VariableRef(owner.model, f"{owner.name}.{local_name}")


def set_properties(reference: ElementRef, **properties: Any) -> ElementRef:
    table = _table_at(reference.model, reference.name.split("."), create=False)
    if table is None:
        raise ValueError(f"model element '{reference.name}' is not defined")
    _update(reference.model, reference.name.split("."), properties)
    return reference


def analysis(model: Model, **properties: Any) -> Model:
    return _update(model, ["analysis"], properties)


def state_selection(model: Model, **properties: Any) -> Model:
    return _update(model, ["state_selection"], properties)


def initial_conditions(model: Model, **properties: Any) -> Model:
    return _update(model, ["initial_conditions"], properties)


def parameters(model: Model, **properties: Any) -> Model:
    return _update(model, ["parameters"], properties)


def graphics(target: Model | ElementRef, **properties: Any) -> Model | ElementRef:
    if isinstance(target, Model):
        return _update(target, ["graphics"], properties)
    _update(target.model, [*target.name.split("."), "graphics"], properties)
    return target


def graphic(reference: ElementRef, name: object, **properties: Any) -> ElementRef:
    local_name = _model_name(name)
    if "." in local_name:
        raise ValueError("a graphic's local name cannot contain a period")
    _update(
        reference.model,
        [*reference.name.split("."), "graphics", local_name],
        properties,
    )
    return reference


def simulation(
    model: Model, *, frames_per_second: Real | None = None, **properties: Any
) -> Model:
    _update(model, ["simulation"], properties)
    if frames_per_second is not None:
        fps = float(frames_per_second)
        if fps <= 0:
            raise ValueError("frames_per_second must be positive")
        table = model._document["simulation"]
        start_time = float(table.get("start_time", 0.0))
        end_time = float(table.get("end_time", 1.0))
        if end_time < start_time:
            raise ValueError("simulation end_time must not precede start_time")
        table["output_samples"] = round((end_time - start_time) * fps) + 1
    return model


def write_model(
    path: str | os.PathLike[str], model: Model, *, overwrite: bool = False
) -> Path:
    """Write the generated portable model document as TOML."""

    destination = Path(path)
    if destination.exists() and not overwrite:
        raise FileExistsError(f"model file already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(_toml.dumps(model._document), encoding="utf-8")
    return destination


def beam_marker(beam: ElementRef, node: object) -> ElementRef:
    """Return one of a flexible beam's generated end_i, cm, or end_j markers."""

    if beam.kind != "flexible_beam":
        raise ValueError("beam_marker requires a flexible_beam handle")
    marker_name = str(node)
    if marker_name not in {"end_i", "cm", "end_j"}:
        raise ValueError("flexible beam marker must be end_i, cm, or end_j")
    return ElementRef(beam.model, f"{beam.name}.{marker_name}", "marker")


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _solver_path(dimension: str, explicit: str | os.PathLike[str] | None) -> Path:
    executable_name = "simp2d" if dimension == "planar" else "simp3d"
    if explicit is not None:
        path = Path(explicit).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"solver executable does not exist: {path}")
        return path

    environment_name = f"SIMP{2 if dimension == 'planar' else 3}D_EXECUTABLE"
    if os.environ.get(environment_name):
        return _solver_path(dimension, os.environ[environment_name])

    root = Path(
        os.environ.get("PRACTICAL_MECHANICAL_SIMULATION_ROOT", _repository_root())
    )
    repository_solver = root / "bin" / executable_name
    if repository_solver.is_file():
        return repository_solver.resolve()

    found = shutil.which(executable_name)
    if found:
        return Path(found).resolve()
    raise FileNotFoundError(
        f"could not find {executable_name}; set {environment_name} or "
        "PRACTICAL_MECHANICAL_SIMULATION_ROOT"
    )


def run_model(
    model: Model,
    *,
    output: str | os.PathLike[str] | None = None,
    duration: Real | None = None,
    samples: int | None = None,
    overwrite: bool = False,
    solver: str | os.PathLike[str] | None = None,
    model_path: str | os.PathLike[str] | None = None,
    cwd: str | os.PathLike[str] | None = None,
    capture_output: bool = True,
    check: bool = True,
) -> SolverRun:
    """Write TOML and invoke the existing dimension-specific solver."""

    if samples is not None and duration is None:
        raise ValueError("duration is required when samples is supplied")
    executable = _solver_path(model.dimension, solver)
    working_directory = Path(cwd or Path.cwd()).resolve()
    output_path = None if output is None else Path(output)
    if output_path is not None and not output_path.is_absolute():
        output_path = (working_directory / output_path).resolve()

    temporary_path: Path | None = None
    if model_path is None:
        model.source_directory.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".toml",
            prefix=f".{model._document['model']['name']}-",
            dir=model.source_directory,
            encoding="utf-8",
            delete=False,
        )
        try:
            handle.write(_toml.dumps(model._document))
        finally:
            handle.close()
        input_path = temporary_path = Path(handle.name)
    else:
        input_path = write_model(model_path, model, overwrite=overwrite)

    command = [str(executable), str(input_path.resolve())]
    if duration is not None:
        command.append(str(float(duration)))
    if samples is not None:
        command.append(str(int(samples)))
    if output_path is not None:
        command.extend(["--output", str(output_path)])
    if overwrite:
        command.append("--overwrite")

    try:
        completed = subprocess.run(
            command,
            cwd=working_directory,
            text=True,
            capture_output=capture_output,
            check=False,
        )
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    run = SolverRun(
        dimension=model.dimension,
        command=tuple(command),
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
        output_path=output_path,
    )
    if check and run.returncode != 0:
        raise SolverError(run)
    return run
