import os
import sys
from setuptools import setup


def main():
        # if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
                # sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
                # sys.exit(1)

        setup(name="IRTools",
              version="1.2.2.1",
              description="a computational toolset for detection and analysis of intron retention from RNA-Seq libraries",
              author='Zhouhao Zeng, Andy Cao',
              author_email='zzhlbj23@gwmail.gwu.edu, acao3@email.gwu.edu',
              url='https://github.com/WeiqunPengLab/IRTools/',
              package_dir={'IRTools' : 'IRTools'},
              packages=['IRTools'],   
              package_data={'IRTools': ['data/*.gtf']},
              include_package_data=True,
              scripts=['bin/IRTools'],
              classifiers=[
                      'Development Status :: 4 - Beta',
                      'Environment :: Console',
                      'Intended Audience :: Developers',
                      'Intended Audience :: Science/Research',              
                      'License :: OSI Approved :: GNU General Public License (GPL)',
                      'Operating System :: MacOS :: MacOS X',
                      'Operating System :: POSIX',
                      'Topic :: Scientific/Engineering :: Bio-Informatics',
                      'Programming Language :: Python',
                      ],
              install_requires=[
                      'numpy',
                      'scipy',
                      'pandas',
                      'networkx',
                      'bx-python',
                      'HTSeq==0.13.5',
                      'pysam==0.19.1'],
              )

if __name__ == '__main__':
        main()
