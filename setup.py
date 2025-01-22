from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'Simulate moire lattice systems in both real and momentum space and calculate various related observables.'
with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()
with open("requirements.txt", "r") as f:
    INSTALL_REQUIRES = f.read().splitlines()

# Setting up
setup(
    name="moirepy",
    version=VERSION,
    author="Aritra Mukhopadhyay, Jabed Umar",
    author_email="amukherjeeniser@gmail.com, jabedumar12@gmail.com",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    url="https://github.com/jabed-umar/MoirePy",
    keywords=['python', 'moire', 'lattice', 'physics', 'materials', 'condensed matter'],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Education",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ]
)
