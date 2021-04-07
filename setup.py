#!/usr/bin/env python

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Requirements
with open("requirements.txt", "r") as fh:
    install_requires = fh.read()

setuptools.setup(
    name="simple_scvi",
    version="0.0.1",
    author="Juho Timonen",
    author_email="juho.timonen@iki.fi",
    description="Simple scVI.",
    zip_safe=False,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jtimonen/simple-scvi",
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    setup_requires=["pip>=19.0.3"],
    license="BSD-3-Clause",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
