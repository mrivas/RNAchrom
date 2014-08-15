from setuptools import setup

setup(name='RNAchrom',
	version='0.1',
	description='Tools for RNA-chromatin interactions',
	url='http://mrivas.hasdocs.com/RNAchrom/',
	author='Marcelo Rivas-Astroza',
	license='GPL',
	packages=['RNAchrom'],
	entry_points={'console_scripts':['annotateBAM=annotateBAM.command_line:main'],},
	install_requires=['numpy','sys','argparse','HTSeq','scipy',],
	dependency_links=['https://github.com/noteed/python-pickle/archive/master.zip'],
	zip_safe=False)
