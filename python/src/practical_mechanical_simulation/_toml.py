"""Small dependency-free TOML writer for generated Sim2D/Sim3D models."""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from typing import Any


_BARE_KEY = re.compile(r"^[A-Za-z0-9_-]+$")


def _key(value: str) -> str:
    return value if _BARE_KEY.fullmatch(value) else json.dumps(value)


def _path(parts: Sequence[str]) -> str:
    return ".".join(_key(part) for part in parts)


def _inline_table(value: Mapping[str, Any]) -> str:
    fields = ", ".join(
        f"{_key(str(key))} = {_value(item)}" for key, item in value.items()
    )
    return "{ " + fields + " }"


def _value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("generated model values must be finite")
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, Mapping):
        return _inline_table(value)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return "[" + ", ".join(_value(item) for item in value) + "]"
    raise TypeError(f"cannot write {type(value).__name__} as TOML")


def dumps(document: Mapping[str, Any]) -> str:
    """Serialize the model document using ordinary TOML tables."""

    lines: list[str] = []

    def emit_table(path: list[str], table: Mapping[str, Any]) -> None:
        if path:
            if lines and lines[-1] != "":
                lines.append("")
            lines.append(f"[{_path(path)}]")

        children: list[tuple[str, Mapping[str, Any]]] = []
        for key, value in table.items():
            if isinstance(value, Mapping):
                children.append((str(key), value))
            else:
                lines.append(f"{_key(str(key))} = {_value(value)}")

        for key, child in children:
            emit_table([*path, key], child)

    emit_table([], document)
    return "\n".join(lines).rstrip() + "\n"
