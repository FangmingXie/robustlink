from setuptools import setup

setup(
    name='robustlink',
    version='0.1.0',
    description='A python package for linking genes and enhancers',
    url='https://github.com/FangmingXie/robustlink',
    author='Fangming Xie, Eran Mukamel',
    author_email='emukamel@ucsd.edu',
    license='MIT',
    packages=['robustlink'],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.8',
    ],
    install_requires=[
        'numpy',
    ]
    )