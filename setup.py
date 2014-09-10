from setuptools import setup

setup(name='RNAchrom',
	version='0.1',
	description='Tools for RNA-chromatin interactions',
	url='http://mrivas.hasdocs.com/RNAchrom/',
	author='Marcelo Rivas-Astroza',
	license='GPL',
	packages=['RNAchrom'],
	scripts=['bin/annotateBAM','bin/countHits','bin/corrWind','bin/sj','bin/specificity','bin/removeOverlapps','bin/detectLinks'],
	install_requires=['numpy','sys','argparse','HTSeq','scipy',],
	zip_safe=False)
