from setuptools import setup, find_packages


setup(
  name='epcsaftpy',
  license='MIT',
  version='0.0.09',
  description='e-PC-SAFT',
  author='Esteban Cea-Klapp',
  author_email='estebancea@udec.cl',
  url='https://github.com/estebancea/epcsaftpy',
  download_url='https://github.com/estebancea/epcsaftpy.git',
  long_description="epcsaft is an open-source python package of e-PC-SAFT",
  packages=find_packages(),
  install_requires=['numba', 'numpy', 'scipy', 'pandas'],
  platforms=["Windows", "Linux", "Mac OS", "Unix"],
  keywords=['PC-SAFT'],
  zip_safe=False
)
