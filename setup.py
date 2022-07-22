import os
import setuptools

def read(fname, mode='r', encoding='utf-8'):
    """ Wrapper for open() """
    return open(os.path.join(os.path.dirname(__file__), fname),
                mode=mode, 
                encoding=encoding)

with read('README.md') as f:
    long_description = f.read()


setuptools.setup(
    name='segy-crooked2D-to-pseudo3D',
    version='0.1.1',
    author='Fynn Warnke',
    author_email='fwrnke@mailbox.org',
    description='Create pseudo-3D cube from crooked 2D SEG-Y file',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    url='https://github.com/fwrnke/segy-crooked2D-to-pseudo3D',
    project_urls={
        'Bug Tracker': 'https://github.com/fwrnke/segy-crooked2D-to-pseudo3D/issues',
    },
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering'
    ],
    packages=setuptools.find_packages(),
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.18.0',
        'segyio>=1.8.0',
    ],
    entry_points={
    'console_scripts': [
        'create_pseudo3D=segy_crooked2D_to_pseudo3D.crooked2D_to_pseudo3D:main',
        ],
    },
)
