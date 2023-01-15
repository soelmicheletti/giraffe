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
    install_requires=["matplotlib>=3.3.4",
                      "netZooPy>=0.8.1",
                      "numpy>=1.19.2",
                      "pandas>=1.4.4",
                      "plotly>=5.10",
                      "sklearn>=0.0",
                      "torch>=1.1"],
)