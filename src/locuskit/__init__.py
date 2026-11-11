"""LOCUSKIT: genomic neighborhoods, collinearity and fractionation analysis."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("locuskit")
except PackageNotFoundError:
    __version__ = "0+unknown"
