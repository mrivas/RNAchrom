from setuptools import setup

setup(name='RNAchrom',
	version='0.1',
	description='Tools for RNA-chromatin interactions',
	url='https://szbio.ucsd.edu/public/RNAchrom/',
	author='Marcelo Rivas-Astroza',
	license='GPL',
	packages=['RNAchrom'],
	scripts=['bin/annotateBAM','bin/countHits','bin/corrWind','bin/sj','bin/specificity','bin/removeOverlapps','bin/detectLinks','bin/filterBAM','bin/countReads','bin/bamFragDistribution','bin/fixBAMmates','bin/countFragments','bin/splitBAM','bin/countReadTypes','bin/geneTargets','bin/countRestrSites','bin/plotCounts','bin/connections','bin/annotateConn','bin/vennDiagram','bin/geneNames','bin/annotateInteractions','bin/countInteractions','bin/addGO'],
	install_requires=['numpy','pysam','argparse','HTSeq','scipy','pybedtools','pandas',],
	zip_safe=False)
