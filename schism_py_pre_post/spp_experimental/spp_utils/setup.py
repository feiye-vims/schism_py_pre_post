from setuptools import setup

with open('requirements.txt', 'r') as f:
  requirements = f.read().splitlines()

  setup(
    name='spp_utils',
    version='0.0.1',
    description='misc utilities',
    license='MIT',
    packages=['spp_utils'],
    package_data={},
    install_requires=requirements,
  )

