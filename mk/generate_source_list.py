#!/usr/bin/env python3

# ======================= python ======================= #
#                         ____ _        ___        __    #
#  _ __   __ _ _ __ __ _ / ___| |      / \ \      / /__  #
# | '_ \ / _` | '__/ _` | |   | |     / _ \ \ /\ / / __| #
# | |_) | (_| | | | (_| | |___| |___ / ___ \ V  V /\__ \ #
# | .__/ \__,_|_|  \__,_|\____|_____/_/   \_\_/\_/ |___/ #
# |_|                                                    #
# ====================================================== #

"""Generate the combined source list for the base and user directories."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


MODULE_RE = re.compile(
    r"^\s*module\s+([a-z]\w*)\s*$",
    re.IGNORECASE,
)
SUBMODULE_RE = re.compile(
    r"^\s*submodule\s*\([^)]*\)\s*([a-z]\w*)",
    re.IGNORECASE,
)


def declared_units(path: Path) -> set[str]:
    units = set()
    with path.open(encoding="utf-8", errors="replace") as source:
        for line in source:
            code = line.split("!", 1)[0]
            for pattern in (MODULE_RE, SUBMODULE_RE):
                match = pattern.match(code)
                if match:
                    units.add(match.group(1).casefold())
    return units


def make_words(paths: list[Path]) -> str:
    escaped = []
    for path in paths:
        value = str(path)
        value = value.replace("$", "$$").replace("#", r"\#")
        value = re.sub(r"([ \t])", r"\\\1", value)
        escaped.append(value)
    return " ".join(escaped)


def generate(user_dir: Path, base_dir: Path, base_sources: list[Path]) -> str:
    user_dir = user_dir.resolve()
    base_dir = base_dir.resolve()
    base_units: dict[str, Path] = {}
    resolved_base_sources = []
    for source in base_sources:
        source = source.expanduser().resolve()
        if not source.is_file():
            raise ValueError(f"Base source does not exist: {source}")
        resolved_base_sources.append(source)
        for unit in declared_units(source):
            previous = base_units.get(unit)
            if previous is not None and previous.resolve() != source.resolve():
                raise ValueError(
                    f"Fortran unit {unit!r} is declared by multiple base sources: "
                    f"{previous} and {source}"
                )
            base_units[unit] = source
    matches: dict[Path, Path] = {}
    user_sources = set()
    duplicate_user_sources = set()
    build_dir = (user_dir / "build").resolve()
    for source in sorted(user_dir.rglob("*")):
        if not source.is_file() or source.suffix.casefold() != ".f90":
            continue
        resolved = source.resolve()
        if (
            resolved == build_dir
            or build_dir in resolved.parents
            or resolved == base_dir
            or base_dir in resolved.parents
        ):
            continue
        user_sources.add(resolved)

        user_units = declared_units(source)
        replaced_sources = {
            base_units[unit]
            for unit in user_units
            if unit in base_units
            and base_units[unit].resolve() != resolved
        }
        if len(replaced_sources) > 1:
            names = ", ".join(str(path) for path in sorted(replaced_sources))
            raise ValueError(
                f"{source} declares units from multiple base sources: {names}"
            )
        if replaced_sources:
            base_source = replaced_sources.pop()
            previous = matches.get(base_source)
            if previous is not None and previous.resolve() != resolved:
                if previous.read_bytes() == resolved.read_bytes():
                    duplicate_user_sources.add(resolved)
                    continue
                raise ValueError(
                    f"Multiple user sources replace {base_source}: "
                    f"{previous} and {source}"
                )
            matches[base_source] = source

    replaced_user_sources = set(matches.values())
    base_final_sources = [
        matches.get(source, source)
        for source in resolved_base_sources
        if source not in matches
    ]
    additional_user_sources = sorted(
        user_sources - replaced_user_sources - duplicate_user_sources
    )
    final_sources = []
    for source in resolved_base_sources:
        final_sources.append(matches.get(source, source))
    final_sources.extend(additional_user_sources)

    return (
        f"FINAL_SRCS := {make_words(final_sources)}\n"
        f"BASE_FINAL_SRCS := {make_words(base_final_sources)}\n"
        f"USER_FINAL_SRCS := "
        f"{make_words(sorted(replaced_user_sources | set(additional_user_sources)))}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--user-dir", type=Path, required=True)
    parser.add_argument("--base-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--base-source", type=Path, action="append", default=[])
    args = parser.parse_args()

    try:
        content = generate(args.user_dir, args.base_dir, args.base_source)
    except ValueError as error:
        parser.error(str(error))
    if not args.output.exists() or args.output.read_text(encoding="utf-8") != content:
        args.output.write_text(content, encoding="utf-8")


if __name__ == "__main__":
    main()
