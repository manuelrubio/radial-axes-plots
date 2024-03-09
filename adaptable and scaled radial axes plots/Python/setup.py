import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="calculate_mapping",
    version="0.0.1",
    author="Lucas Ramos and Manuel Rubio",
    author_email="manuel.rubio@urjc.es",
    description="Computes the mappings of several radial axes methods",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)