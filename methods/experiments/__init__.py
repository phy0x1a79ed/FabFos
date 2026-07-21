"""Experiment specs. One module per experiment, each defining SPEC: ExperimentSpec.

A spec names inputs and targets. It imports `canon` and computes; it never
transcribes a basis value (see ../canon.py, and the tripwire that enforces it).
"""
