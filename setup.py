import subprocess
from codecs import open
from os import path

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install


class PostDevelopCommand(develop):
    """Post-installation for development mode."""

    def run(self):
        subprocess.check_call(['npm', 'i'])
        subprocess.check_call(['npm', 'run-script', 'build'])
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):
        subprocess.check_call(['npm', 'i'])
        subprocess.check_call(['npm', 'run-script', 'build'])
        install.run(self)


with open(path.join(path.abspath(path.dirname(__file__)), "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

requires = [
    "anndata",
    "CacheControl",
    "flask",
    "flask-compress",
    "fsspec",
    "gunicorn",
    "loompy",
    "numpy",
    "pandas>=1.0",
    "pyarrow",
    "pymongo",
    "scipy",
    "zarr"
]

setup(
    name="cirrocumulus",
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    # cmdclass={
    #         'develop': PostDevelopCommand,
    #         'install': PostInstallCommand,
    # },
    # version="2.0.2",
    description="Single-cell visualization application",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/klarman-cell-observatory/cirrocumulus",
    author="Joshua Gould",
    author_email='jgould@broadinstitute.org',
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization"
    ],
    extras_require=dict(
        test=['pytest', 'scanpy']
    ),
    keywords="single cell/nucleus genomics visualization",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    python_requires=">=3.6",
    package_data={
    },
    entry_points={"console_scripts": ["cirro=cirrocumulus.__main__:main"]},
)
