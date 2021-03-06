from setuptools import setup, find_packages

setup(name='ogrdbweb',
      version='0.1.0',
      description='Immunoglobin germline web interface',
      author='Chris Warth <cwarth@fredhutch.org>',
      packages=find_packages(exclude=['contrib', 'docs', 'tests']),
      install_requires=[
          'flask','frozen-flask','Flask-Breadcrumbs','numpy','biopython'
      ],

      package_data={'ogrdbweb': ['templates/*', 'static/*']}
)
