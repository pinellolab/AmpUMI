from setuptools import setup, find_packages

setup(
    name='AmpUMI',
    version='1.0.0',
    description='Design and analysis of unique molecular indentifiers for deep amplicon sequencing', 
    url = 'https://github.com/pinellolab/AmpUMI.git',
    author='Kendell Clement',
    author_email = 'k.clement.dev@gmail.com',
    packages=find_packages(),
    install_requires=['sympy','mpmath','numpy']
)
