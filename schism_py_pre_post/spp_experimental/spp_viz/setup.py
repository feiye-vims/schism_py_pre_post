from setuptools import setup, find_namespace_packages

with open('requirements.txt', 'r') as f:
  requirements = f.read().splitlines()

  setup(
    name='spp_viz',
    version='0.0.1',
    description='A package of modules used for visualization',
    license='MIT',
    packages=find_namespace_packages(include=['spp.spp_viz.*']),
    package_data={},
    install_requires=requirements,
  )

