from setuptools import setup

setup(
    name='anatools',
    version='0.1.0',
    author='Gabriel Moreira, Gilson Correia, Matheus Costa, Raphael Souza',
    author_email='gabrielmscampos@gmail.com, gilson.correiasilva@gmail.com, mataiascost@gmail.com, trivium.raphael@gmail.com',
    packages=[
    	'anatools',
        'anatools.plot'
    ],
    package_dir={
    	'anatools': 'src',
        'anatools.plot': 'src/plot'
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
    ],
)