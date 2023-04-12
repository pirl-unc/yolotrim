set -e

maturin build
pip install -e .
