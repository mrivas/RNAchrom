from setuptools import setup

setup(name='funniest',
	version='0.1',
	description='Tools for RNA-chromatin interactions',
	url='https://github.com/mrivas/RNAchrom',
	author='Marcelo Rivas-Astroza',
	license='GPL',
	packages=['RNAchrom'],
	scripts=['bin/annotateBAM.py'],
	install_requires=['argparse','sys','pysam','HTSeq','numpy','re'],
	zip_safe=False)
