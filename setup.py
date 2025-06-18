import subprocess
from setuptools import setup, find_packages


VERSION = subprocess.run(["git", "describe", "--tags"], stdout=subprocess.PIPE).stdout.decode("utf-8").strip()

if VERSION == "":
    pwd = subprocess.run(["pwd"], stdout=subprocess.PIPE).stdout.decode("utf-8").strip()
    VERSION = pwd.split("/")[-1][8:]  # fallback to the last part of the current directory name, e.g., "moirepy-1.3.3" -> "1.3.3"

if "-" in VERSION:
    # when not on tag, git describe outputs: "1.3.3-22-gdf81228"
    # pip has gotten strict with version numbers
    # so change it to: "1.3.3+22.git.gdf81228"
    # See: https://peps.python.org/pep-0440/#local-version-segments
    v, i, s = VERSION.split("-")  # v = version, i = number of commits since tag, s = sha
    VERSION = v + "+" + i + ".git." + s

with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()
with open("requirements.txt", "r") as f:
    INSTALL_REQUIRES = f.read().splitlines()

setup(
    name="moirepy",
    version=VERSION,  # obtained automatically through github release tags
    license="MIT",
    author="Aritra Mukhopadhyay, Jabed Umar",
    author_email="amukherjeeniser@gmail.com, jabedumar12@gmail.com",
    description='Simulate moire lattice systems in both real and momentum space and calculate various related observables.',
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,  # read from README.md
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,  # read from requirements.txt
    url="https://github.com/jabed-umar/MoirePy",
    keywords=[
        'python', 'moire', 'moiré', 'moire lattice', 'twistronics',
        'bilayer graphene', 'tight binding', 'lattice simulation',
        'physics', 'material science', 'condensed matter', 'k-space', 'real-space'
    ],
    python_requires=">=3.7",
    project_urls={
        "Documentation": "https://jabed-umar.github.io/MoirePy/",
        "Source Code": "https://github.com/jabed-umar/MoirePy",
        "Bug Tracker": "https://github.com/jabed-umar/MoirePy/issues",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
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
