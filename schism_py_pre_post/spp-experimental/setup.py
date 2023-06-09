from setuptools import setup

with open('requirements.txt', 'r') as f:
  requirements = f.read().splitlines()

  setup(
    name='spp-essentials',
    version='0.0.1',
    description='An sub package of the experimental package of "Python tools for pre/post-processing SCHISM models"',
    license='MIT',
    packages=[
      'spp-essentials',
    ],
    package_data={},
    install_requires=requirements,
  )

