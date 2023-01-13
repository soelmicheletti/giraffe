from setuptools import find_packages, setup

setup(
    name = "giraffe",
    packages = find_packages(exclude=[]),
    version = "0.0.1",
    license = "MIT",
    description = "GIRAFFE: biologically informed inference of gene regulatory networks",
    author = "Soel Micheletti",
    author_email = "smicheletti@hsph.harvard.edu",
    url = "https://github.com/soelmicheletti/giraffe",
    keywords=[
        "gene regulatory networks",
        "matrix factorization",
    ],
    install_requires=["torch>=1.1", "numpy>=1.19.2", "matplotlib>=3.3.4", "pandas>=1.4.4"],
)