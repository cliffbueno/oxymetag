from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="oxymetag",
    version="1.1.2",
    author="Clifton P. Bueno de Mesquita",
    author_email="cliff.buenodemesquita@colorado.edu",
    description="Oxygen metabolism profiling from metagenomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cliffbueno/oxymetag",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    entry_points={
    "console_scripts": [
        "oxymetag=oxymetag.cli:main",
    ],
},
    package_data={
        "oxymetag": [
            "data/*",
            "scripts/*",
        ],
    },
    include_package_data=True,
)
