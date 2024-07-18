from setuptools import setup, find_packages

setup(
    name='helixbusters',
    version='0.1.0',
    description='A package for identifying and repairing double-strand breaks in DNA',
    author='Simone Puccio',
    author_email='simone.puccio@humanitasresearch.it',
    url='https://github.com/yourusername/helixbusters',
    packages=find_packages(),
    install_requires=[
        # List of dependencies
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)