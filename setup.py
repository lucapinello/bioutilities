from distutils.core import setup
setup(name='bioutilities',
      version='0.2',
      py_modules=['bioutilities'],
      description='A collection of tools for biological computation',
      install_requires=['numpy','scipy','xlrd','xlwt','blist','pysam'],
      )