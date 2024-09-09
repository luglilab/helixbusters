from setuptools import setup, find_packages

setup(
    name='helixbusters',
    version='0.1.0',
    description='A package for identifying and repairing double-strand breaks in DNA',
    author='Simone Puccio',
    author_email='simone.puccio@humanitasresearch.it',
    url='https://github.com/luglilab/helixbusters',
    packages=find_packages(),
    install_requires=[
        "pandas==2.2.2",
        "numpy==2.0.0",
        "xlrd>=2.0.1",
        "biopython>=1.79",
        "cutadapt==4.9",
        "requests==2.22.0",
        "pysam"
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
)