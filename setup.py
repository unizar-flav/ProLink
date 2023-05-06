
import os
from setuptools import setup


# get text of README.md
current_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(current_path, "README.md")) as f:
    readme_text = f.read()

setup(
    name="ProLink",
    version="0.1.0",
    description="Execute multiple proteomic analysis tools automatically",
    long_description=readme_text,
    long_description_content_type="text/markdown",
    url="https://github.com/unizar-flav/ProLink",
    author="Víctor Sanz, Sergio Boneta",
    license="GPLv3",
    python_requires='>=3.6',
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    install_requires=["biopython"],
    packages=["ProLink", "ProLink.modules"],
    package_data={"ProLink": ["parameters.yaml", "mega_configs/*"]},
    include_package_data=True,
    entry_points={
        "console_scripts": ["prolink=ProLink.__main__:main"]
    }
)
