from setuptools import find_packages, setup

setup(
    name='htools',
    packages=find_packages(include=['htools','htools.*']),
    version='0.1.0',
    description='Research Tools for Hollis Akins',
    author='Hollis Akins',
    license='MIT',
    install_requires=['astropy','numpy','matplotlib','scipy','tqdm','pandas','rich','dotmap'],
)
