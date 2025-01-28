from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="autoQuant",
    version="1.0",
    author="Jorn Steen",
    author_email="jorn.steen@kemi.uu.se",
    description="GUI for quantum yield calculation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JornSteen/autoQuant",
    install_requires=[
        'numpy=>2.2.2',
        'scipy=>1.15.1',
        'lmfit=>1.3.2',
        'pandas=>2.2.3',
        'PyQt5',
        'matplotlib=>3.10.0'
    ],
    
    #extras_require={
    #},
    
    packages=find_packages(
		where='.',
		include=[
		'PowerProcessing',
		'QY',
		'tools',
		'UIs']),
    package_dir={'':'.'},
  
    ## try out later: to make a command-line script (need def main() in main.py)
    #entry_points={
    #    'console_scripts': [
    #        'autoquant=main:main',
    #    ],
    #},
    
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    python_requires='=>3.12',
    keywords=["science", "quantum yield", "QY"],
)
