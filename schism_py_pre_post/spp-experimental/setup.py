from setuptools import setup

with open('requirements.txt', 'r') as f:
  requirements = f.read().splitlines()

  setup(
    name='spp-experimental',
    version='0.0.1',
    description='An experimental package of "Python tools for pre/post-processing SCHISM models"',
    packages=['my_package'],
    license='MIT',
    packages=[
      'spp-experimental',
      'spp-experimental.essentials',
    ],
    package_data={},
    install_requires=requirements,
  )

)
