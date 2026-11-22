"""Read-only adapters for existing ClockBase audit reports."""

from .measurement_pilot import adapt_measurement_pilot
from .protected_repair import adapt_protected_repair

__all__ = ["adapt_measurement_pilot", "adapt_protected_repair"]
