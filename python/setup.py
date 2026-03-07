"""Backward-compatible setup.py for editable installs.

Prefer ``pip install -e .`` which uses pyproject.toml via PEP 517.
"""

from setuptools import setup, find_packages

setup(
    name="mplus-registration",
    version="2.0.0",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=["numpy>=1.21"],
)
