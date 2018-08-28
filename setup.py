from setuptools import setup, find_packages

setup(
    name = "aliquot-maf-tools",
    author = "Kyle Hernandez",
    version = 0.1,
    description = "Tools for creating and filtering aliquot-level MAFs",
    license = "Apache 2.0",
    packages = find_packages(),
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
    entry_points= {
        'console_scripts': [
            'aliquot-maf-tools = aliquotmaf.__main__:main'
        ]
    }
)
