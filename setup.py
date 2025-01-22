import subprocess
from setuptools import setup, find_packages


VERSION = subprocess.run(["git", "describe", "--tags"], stdout=subprocess.PIPE).stdout.decode("utf-8").strip()

if VERSION == "":
    pwd = subprocess.run(["pwd"], stdout=subprocess.PIPE).stdout.decode("utf-8").strip()
    VERSION = pwd.split("/")[-1][8:]

if "-" in VERSION:
    # when not on tag, git describe outputs: "1.3.3-22-gdf81228"
    # pip has gotten strict with version numbers
    # so change it to: "1.3.3+22.git.gdf81228"
    # See: https://peps.python.org/pep-0440/#local-version-segments
    v, i, s = VERSION.split("-")
    VERSION = v + "+" + i + ".git." + s


DESCRIPTION = 'Simulate moire lattice systems in both real and momentum space and calculate various related observables.'
with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()
with open("requirements.txt", "r") as f:
    INSTALL_REQUIRES = f.read().splitlines()


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
