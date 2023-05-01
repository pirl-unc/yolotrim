set -e

maturin build --release
pip install -e .
