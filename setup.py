from setuptools import setup, find_packages

#!/usr/bin/env python

import importlib
import os
import subprocess
from textwrap import dedent
from types import SimpleNamespace

from setuptools import Command, find_packages, setup

GIT_REPO = "aliquot-maf-tools"
PACKAGE = "aliquotmaf"

PYPI_REPO = "bioinf-{}".format(PACKAGE)
GIT_REPO_URL = "https://github.com/NCI-GDC/{}".format(GIT_REPO)

INSTALL_REQUIRES = []

TESTS_REQUIRES = [
    'mock',
    'pytest',
    'pytest-cov',
]

DEV_REQUIRES = [
    'detect-secrets==0.13.1',
    'isort',
    'flake8',
    'pre-commit',
]


GIT_COMMANDS = SimpleNamespace(
    branch=["git", "rev-parse", "--abbrev-ref", "HEAD"],
    commit=["git", "rev-list", "--count", "HEAD"],
    hash=["git", "rev-parse", "HEAD"],
    shorthash=["git", "rev-parse", "--short", "HEAD"],
)



try:
    # Set versions if version file exists
    mod = importlib.import_module("{}".format(PACKAGE))
    __pypi_version__ = mod.__version__
except Exception:
    # Set defaults otherwise
    __pypi_version__ = '0.0.0'


class PrintVersion(Command):
    description = "Print out specified version, default long version."
    user_options = [
        ("pypi", None, "Print package version."),
        ("short", None, "Print semantic version."),
        ("hash", None, "Print commit hash."),
    ]

    def initialize_options(self):
        self.pypi = False
        self.short = False
        self.hash = False

    def finalize_options(self):
        pass

    def run(self):
        if self.pypi:
            print(__pypi_version__)
        elif self.short:
            print(__short_version__)
        elif self.hash:
            try:
                commit_hash = call_subprocess(GIT_COMMANDS.hash)
            except Exception:
                print('')
            else:
                print(commit_hash)
        else:
            print(__long_version__)




def call_subprocess(cmd: list):
    """Return stdout of given command."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    stdout, _ = p.communicate()
    return stdout.decode().strip()


class Requirements(Command):
    description = "foobar"
    user_options = [
        ("install", None, "Bundles only install requirements"),
        ("test", None, "Bundles only install requirements"),
        ("dev", None, "Bundles all requirements"),
    ]

    def initialize_options(self):
        self.install = False
        self.test = False
        self.dev = False

    def finalize_options(self):
        assert self.install + self.test + self.dev == 1, "Provide only one arg"

    def run(self):
        path = os.path.join(".", "requirements.in")
        if self.dev:
            reqs = INSTALL_REQUIRES + TESTS_REQUIRES + DEV_REQUIRES
        elif self.test:
            reqs = INSTALL_REQUIRES + TESTS_REQUIRES
        elif self.install:
            reqs = INSTALL_REQUIRES
        self.write_requirements(path, reqs)
        return

    def write_requirements(self, path, reqs):
        with open(path, "w") as fh:
            fh.write("\n".join(reqs) + "\n")


setup(
    name=PYPI_REPO,
    description="Mutation Annotation Format (MAF) library",
    version=__pypi_version__,
    url=GIT_REPO_URL,
    python_requires=">=3.6",
    setup_requires=['setuptools_scm'],
    use_scm_version={
        "write_to": os.path.join(f"{PACKAGE}", "_version.py"),
        "fallback_version": __pypi_version__,
    },
    packages=find_packages(),
    package_data={'maflib': ['schemas/*.json']},
    install_requires=INSTALL_REQUIRES,
    tests_require=TESTS_REQUIRES,
    cmdclass={
        "capture_requirements": Requirements,
        "print_version": PrintVersion,
    },
    include_package_data=True,
    entry_points= {
        'console_scripts': [
            'aliquot-maf-tools = aliquotmaf.__main__:main'
        ]
    }
)

# __END__
