from setuptools import setup, find_packages


setup(
    name='fusion_connect',
    license='MIT',
    version='0.1.0',
    description='Lightweight Python client for the Autodesk Platform Services (APS) REST API used by Autodesk Fusion.',
    long_description=(
        'fusion_connect is a small Python package that authenticates against '
        'Autodesk Platform Services (APS, formerly Forge) and exposes a thin '
        'wrapper around the Data Management API used by Autodesk Fusion.'
    ),
    packages=find_packages(),
    install_requires=['requests>=2.25'],
    python_requires='>=3.8',
    platforms=['Windows', 'Linux', 'Mac OS', 'Unix'],
    keywords=['autodesk', 'fusion', 'aps', 'forge', 'cad'],
    zip_safe=False,
)
