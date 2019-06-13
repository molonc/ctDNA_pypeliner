from setuptools import setup, find_packages
import versioneer

setup(
	name='ctdna_pypeliner',
	packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
	description='ctDNA pipeline using multiple sequence analysis tools',
	author='Patricia Ye',
	author_email='patriciaye99@gmail.com',
	url='https://github.com/shahcompbio/ctDNA_pypeliner',
	entry_points={'console_scripts': ['ctdna_pypeliner = ctDNA.run:main']},
    package_data={'':['scripts/*.py', 'scripts/*.R', 'scripts/*.npz', 'scripts/*.pl', "config/*.yaml", "data/*"]}
	)