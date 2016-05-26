from distutils.core import setup

setup(
    name='4SUDRB-seq_hmm',
    version='1.0',
    url='https://github.com/eranroz/4suDRB-seq',
    license='MIT',
    author='eranroz',
    author_email='eranroz@cs.huji.ac.il',
    description='simple HMM for 4SUDRB-seq. see readme',
    dependency_links=[
        'https://github.com/eranroz/hmm',
    ],
    install_requires=['pandas']
)
