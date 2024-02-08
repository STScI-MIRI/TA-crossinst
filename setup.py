from setuptools import setup

setup(
     # Needed to silence warnings (and to be a worthwhile package)
    name='jwst_ta',
    url='https://github.com/STScI-MIRI/TA-crossinst',
    # Needed to actually package something
    packages=['jwst_ta'],
    # Needed for dependencies
    install_requires=[
        'numpy',
        'astropy',
        'matplotlib',
    ],
    # *strongly* suggested for sharing
    version='0.1',
    # The license can be anything you like
    license='MIT',
    description='Cross-instrument TA algorithm simulator',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.md').read(),   
)
