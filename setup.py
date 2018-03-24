from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

version_file = open(os.path.join(here, 'hydroid_gui', 'VERSION'))
version = version_file.read().strip()

# Get the long description from the relevant file
with codecs.open(os.path.join(here, 'DESCRIPTION.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='hydroid_gui',

    # Versions should comply with PEP440. For single-sourced versioning, see
    # http://packaging.python.org/en/latest/tutorial.html#version
    version=version,

    description='GUI for HYDROID - Python package for analyzing hydroxyl-radical footprinting experiments of protein-DNA complexes',
    long_description=long_description,

    # The project URL.
    url='https://github.com/intbio/HYDROID_GUI',

    # Author details
    author='Grigoriy A. Armeev, Alexey K. Shaytan, David Landsman, Anna R. Panchenko',
    author_email='armeev@intbio.org',

    # Choose your license
    license='Public Domain',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        # Project maturity. 
        'Development Status :: 3 - Alpha',

        # Intended audience
        'Intended Audience :: Science/Research',

        # Topic
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # License should match "license" above
        'License :: GPL2',

        # Python versions supported
        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='science hydroxyl-radical footprinting image analysis',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['data', 'hydroid_gui_conda']),

    # Run-time package dependencies. These will be installed by pip when your
    # project is installed.
    install_requires=[
        'hydroid',
        'pyqt',
        'billiard',
    ],
    # Data files included in your packages. If using Python 2.6 or less, 
    # then these have to be included in MANIFEST.in as well.
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'HYDROID_GUI=hydroid_gui:main,
        ],
    },
    python_requires='==2.7.*',
    # MANIFEST.in included entries should be included as package data and
    # installed into site-packages 
    include_package_data=False,

    # Default to False unless you specifically intend the package to be
    # installed as an .egg
    zip_safe=False,
)
