from setuptools import setup

setup(name='RNAchrom',
	version='0.1',
	description='Tools for RNA-chromatin interactions',
	url='http://mrivas.hasdocs.com/RNAchrom/',
	author='Marcelo Rivas-Astroza',
	license='GPL',
	packages=['RNAchrom'],
	scripts=['bin/annotateBAM.py','bin/countHits.py','bin/corrWind.py','bin/sj.py','bin/specificity.py','bin/removeOverlapps.py'],
	install_requires=['numpy','sys','argparse','HTSeq','scipy',],
	zip_safe=False)
