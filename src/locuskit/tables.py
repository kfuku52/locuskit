"""Strict TSV boundaries shared by the LOCUSKIT commands."""

from __future__ import annotations

import csv
import gzip
import io
import math
import os
import re
import tempfile
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

_CONTROL = re.compile(r"[\x00-\x1f\x7f]")
_BAD_IDENTIFIER = re.compile(r"[\s\x00-\x1f\x7f]")
_INTEGER = re.compile(r"-?[0-9]+")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


@contextmanager
def open_output(path: Path) -> Iterator[TextIO]:
    if path.suffix != ".gz":
        with path.open("w", encoding="utf-8", newline="") as stream:
            yield stream
        return
    with (
        path.open("wb") as raw,
        gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped,
    ):
        with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
            yield stream


def read_table(
    path: Path, required: tuple[str, ...], *, exact: bool = False
) -> Iterator[dict[str, str]]:
    with open_text(path) as stream:
        reader = csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar=None)
        header = next(reader, [])
        if (
            not header
            or len(header) != len(set(header))
            or not all(header)
            or not set(required).issubset(header)
            or (exact and tuple(header) != required)
        ):
            raise ValueError(f"{path}: invalid header, expected {', '.join(required)}")
        for line, row in enumerate(reader, 2):
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"{path}:{line}: expected {len(header)} columns, got {len(row)}")
            if _CONTROL.search(" ".join(row)):
                raise ValueError(f"{path}:{line}: control character in TSV field")
            yield dict(zip(header, row, strict=True))


@contextmanager
def staged_text(path: Path, inputs: Iterable[Path]) -> Iterator[TextIO]:
    path = path.absolute()
    for source in inputs:
        if path.resolve() == source.resolve() or (
            path.exists() and source.exists() and path.samefile(source)
        ):
            raise ValueError(f"Output would overwrite an input: {path}")
    if path.is_symlink():
        raise ValueError(f"Output is a symlink: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=path.suffix
    )
    os.close(descriptor)
    staged = Path(name)
    try:
        with open_output(staged) as stream:
            yield stream
        os.replace(staged, path)
    finally:
        staged.unlink(missing_ok=True)


def write_table(
    path: Path, header: tuple[str, ...], rows: Iterable[Iterable[object]], inputs: Iterable[Path]
) -> None:
    path = path.absolute()
    with staged_text(path, inputs) as stream:
        writer = csv.writer(
            stream, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None
        )
        writer.writerow(header)
        for row in rows:
            fields = ["" if value is None else str(value) for value in row]
            if len(fields) != len(header) or any(_CONTROL.search(value) for value in fields):
                raise ValueError(f"{path}: output row has invalid width or control characters")
            writer.writerow(fields)


def identifier(value: str, field: str) -> str:
    if not value or _BAD_IDENTIFIER.search(value):
        raise ValueError(f"{field} must be a nonempty identifier without whitespace: {value!r}")
    return value


def integer(value: str, field: str, *, minimum: int | None = None) -> int:
    if not _INTEGER.fullmatch(value):
        raise ValueError(f"{field} must be an integer: {value!r}")
    number = int(value)
    if minimum is not None and number < minimum:
        raise ValueError(f"{field} must be at least {minimum}: {value!r}")
    return number


def evalue(value: str, field: str) -> float:
    try:
        result = float(value)
    except ValueError:
        raise ValueError(f"{field} must be a finite nonnegative E-value: {value!r}") from None
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{field} must be a finite nonnegative E-value: {value!r}")
    return result
