from setuptools import setup

setup(
	name = "ngg2",
	packages = ['ngg2'],
	version = "1.3.0",
	description = 'Tool to identify NGGNGG Cas9 gRNA sites in any indexed FASTA file.',
	author = "Elisha Roberson",
	author_email = 'dr.eli.roberson@gmail.com',
	url = 'https://github.com/RobersonLab/ngg2',
	license = 'MIT',
	classifiers=[
	"Development Status :: 5 - Production/Stable",
	"Environment :: Console",
	"Intended Audience :: Science/Research",
	"License :: OSI Approved :: MIT License",
	"Topic :: Scientific/Engineering :: Bio-Informatics",
	"Programming Language :: Python :: 2.7",
	"Programming Language :: Python :: 3",
	"Programming Language :: Python :: 3.2",
	"Programming Language :: Python :: 3.3",
	"Programming Language :: Python :: 3.4",
	"Programming Language :: Python :: 3.5"
	],
	keywords="CRISPR Cas9 3'GG site finder",
	install_requires = ['pyfaidx', 'regex>=2016.01.10', 'six'],
	entry_points = {'console_scripts':["ngg2 = ngg2.__main__:run"]},
	test_suite = 'nose.collector',
	tests_require = ['nose']
)
