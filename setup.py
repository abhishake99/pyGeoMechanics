import setuptools
 
with open("README.md", "r") as fh:
    long_description = fh.read()
 
setuptools.setup(
    name="pyGeoMechanics",
    version="0.0.1",
    author="Abhishek Ramawat",
    author_email="ramawatabhishek22@gmail.com",
    description="Package for Geomechanics",
    license="Apache License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=[

        'scipy',
        'pandas',
        'numpy',
        'matplotlib',
        'lasio',
        'wellpathpy'

    ],
    python_requires='>=3.8'
)