from setuptools import setup

with open('requirements.txt', 'r') as f:
  requirements = f.read().splitlines()

  setup(
    name='spp_essentials',
    version='0.0.1',
    description='The essentials of spp that require minimum dependencies',
    license='MIT',
    packages=['spp_essentials'],
    package_data={},
    install_requires=requirements,
  )

