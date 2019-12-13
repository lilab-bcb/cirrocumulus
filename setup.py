from codecs import open
from os import path

from setuptools import setup, find_packages


with open(path.join(path.abspath(path.dirname(__file__)), "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

requires = [
        "anndata",
        "flask",
        "flask-compress",
        "fsspec",
        "natsort",
        "ujson"
]

setup(
    name="cirrocumulus",
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    description="scRNA-Seq visualization tool",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/klarman-cell-observatory/cirrocumulus",
    author="Joshua Gould",
    author_email='jgould@broadinstitute.org',
    classifiers=[  # https://pypi.python.org/pypi?%3Aaction=list_classifiers
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Framework :: Jupyter",
            "Natural Language :: English",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: POSIX :: Linux",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    extras_require=dict(
        test=[
                'pytest'
        ]
    ),
    keywords="single cell/nucleus genomics visualization",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    python_requires="~=3.5",
    package_data={
    },
    entry_points={"console_scripts": ["cirro=cirro.__main__:main"]},
)
