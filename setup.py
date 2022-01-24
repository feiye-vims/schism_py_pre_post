import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
  long_description = fh.read()

  setuptools.setup(
  name='schism_py_pre_post',
  version='0.0.1',
  author='Fei Ye',
  author_email='feiye@vims.edu',
  description='Python tools for pre/post-processing SCHISM models',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='',
  project_urls = {
    "Issues": ""
  },
  license='MIT',
  packages=['schism_py_pre_post', 'Download', 'Plot', 'Timeseries', 'Shared_modules'],
  install_requires=['climata'],
)
