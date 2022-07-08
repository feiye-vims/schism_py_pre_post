import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
  long_description = fh.read()

  setuptools.setup(
  name='schism_py_pre_post',
  version='0.0.4.3',
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
  packages=[
    'schism_py_pre_post',
    'schism_py_pre_post.Grid',
    'schism_py_pre_post.Plot',
    'schism_py_pre_post.Timeseries',
    'schism_py_pre_post.Shared_modules',
    'schism_py_pre_post.Download',
    'schism_py_pre_post.Geometry',
  ],
  package_data={'schism_py_pre_post': ['Datafiles/*']},
  install_requires=[
    'numpy',
    'pandas',
    'xarray',
    'climata==0.5.0',
    'gsw==3.4.0',
    'pyshp==2.1.3'
  ],
)
