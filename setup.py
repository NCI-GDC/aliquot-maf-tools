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

INSTALL_REQUIRES = ["bioinf-maflib>=1.5.0", "pysam", "pandas", "pyarrow==4.0.*"]

TESTS_REQUIRE = [
    "mock",
    "pytest",
    "pytest-cov",
]

DEV_REQUIRES = [
    "detect-secrets==0.13.1",
    "isort",
    "flake8",
    "pre-commit",
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
    __pypi_version__ = "0.0.0"


class PrintVersion(Command):
    description = "Print out specified version, default long version."
    user_options = [
        ("pypi", None, "Print package version."),
        ("docker", None, "Print Docker-friendly package version."),
        ("hash", None, "Print commit hash."),
    ]

    def initialize_options(self):
        self.pypi = False
        self.docker = False
        self.hash = False

    def finalize_options(self):
        pass

    def run(self):
        if self.pypi:
            print(__pypi_version__)
        elif self.docker:
            print(__pypi_version__.replace("+", "."))
        elif self.hash:
            try:
                commit_hash = call_subprocess(GIT_COMMANDS.hash)
            except Exception:
                print("")
            else:
                print(commit_hash)
        else:
            print(__pypi_version__)


def call_subprocess(cmd: list):
    """Return stdout of given command."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    stdout, _ = p.communicate()
    return stdout.decode().strip()


class Requirements(Command):
    description = "Write out dev-requirements.in file."
    user_options = [
        ("dev", None, "Bundles all requirements"),
    ]

    def initialize_options(self):
        self.dev = False

    def finalize_options(self):
        pass

    def run(self):
        REQUIREMENT = ["-c requirements.txt"]
        if self.dev:
            reqs = REQUIREMENT + DEV_REQUIRES + TESTS_REQUIRE
            path = "dev-requirements.in"
        else:
            raise ValueError("Choose one of install, test, or dev")
        self.write_requirements(path, reqs)
        return

    def write_requirements(self, path, reqs):
        with open(path, "w") as fh:
            fh.write("\n".join(reqs) + "\n")


setup(
    name=PYPI_REPO,
    description="Tools for creating and filtering aliquot-level MAFs",
    version=__pypi_version__,
    url=GIT_REPO_URL,
    python_requires=">=3.6",
    setup_requires=["setuptools_scm"],
    use_scm_version={
        "write_to": os.path.join(PACKAGE, "_version.py"),
        "fallback_version": __pypi_version__,
    },
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    tests_require=TESTS_REQUIRE,
    cmdclass={"capture_requirements": Requirements, "print_version": PrintVersion},
    include_package_data=True,
    scripts=[os.path.join(os.path.dirname(__file__), 'bin', PACKAGE)],
    entry_points={"console_scripts": ["aliquot-maf-tools = aliquotmaf.__main__:main"]},
)

# __END__
