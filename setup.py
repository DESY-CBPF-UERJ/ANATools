from setuptools import setup

from src import __version__

setup(
    name='anatools',
    version=__version__,
    author='Gabriel Moreira, Gilson Correia, Matheus Costa, Raphael Souza',
    author_email='gabrielmscampos@gmail.com, gilson.correiasilva@gmail.com, mataiascost@gmail.com, trivium.raphael@gmail.com',
    packages=[
    	'anatools',
    	'anatools.analysis',
    	'anatools.data'
    ],
    package_dir={
    	'anatools': 'src',
    	'anatools.analysis': 'src/analysis',
    	'anatools.data': 'src/data'
    },
    license='LICENSE',
    description='HEP Analysis tools.',
    keywords=['HEP', 'Analysis'],
    install_requires=[
        "awkward0>=0.15.5",
        "cachetools>=4.2.1",
        "cycler>=0.10.0",
        "kiwisolver>=1.3.1",
        "matplotlib>=3.4.1",
        "mplhep>=0.3.2",
        "mplhep-data>=0.0.2",
        "numpy>=1.20.2",
        "packaging>=20.9",
        "pandas>=1.2.3",
        "Pillow>=8.2.0",
        "pyparsing>=2.4.7",
        "python-dateutil>=2.8.1",
        "pytz>=2021.1",
        "six>=1.15.0",
        "uhi>=0.2.1",
        "uproot3>=3.14.4",
        "uproot3-methods>=0.10.1",
        "tqdm>=4.59.0",
        "mpmath>=1.1.0",
        "futures>=3.1.1",
        "scikit-learn>=0.24.2",
        #"json>=0.0.0",
    ],
)

with open("./install.sh", "w") as f:
    f.write(f"pip3 install ./dist/anatools-{__version__}.tar.gz")
