# segy-crooked2D-to-pseudo3D

This package allows the user to create a pseudo-3D cube from a crooked 2D input SEG-Y file. The 2D input data (traces x samples) is expanded in the third dimension to create a 3D cube (inlines x crosslines x samples) where each inline of the cube is identical to the original 2D input data. Every crossline (xline) therefore corresponds to the original trace or CDP from the 2D input.

**NOTE:** Input SEG-Y coordinates should be in a projected coordinate reference system (e.g. UTM)!

## Installation

Generally, I recommend using a `conda` environment with `pip` installed. To install this package directly from GitHub, run the following line:

```python
pip install git+https://github.com/fwrnke/segy-crooked2D-to-pseudo3D.git
```

You can also download the repository as a ZIP archive. To install `segy-crooked2D-to-pseudo3D` locally using `pip` after unzipping the downloaded archive:

```bash
>>> cd ./segy-crooked2D-to-pseudo3D  # root directory of unzipped package
>>> pip install [-e] .         # -e: setuptools "develop mode"
```

## Usage

### Command line interface (CLI)

`segy-crooked2D-to-pseudo3D` was designed to be mainly used from the command line and provides an interface that can be used by running:

```bash
>>> create_pseudo3D {input_path} [optional parameters]
```

This script requires one positional argument:

- `input_path`: can be either (**A**) input SEG-Y file path or (**B**) input directory of several SEG-Y files to process

Optionally, the following parameter can be used:

- `--help`, `-h`: show help
- `--output_dir`, `-o`: Output directory for pseudo-3D SEG-Y cube(s). Identical to input directory if not specified.
- `--suffix`, `-s`: File suffix. Only used when ``input_path`` is a directory
- `--txt_suffix`: Suffix to append to output filename (default: "_pseudo-3D")
- `--n_ilines`: Number of inlines for pseudo-3D cube (should be **odd** number!)
- `--spacing`: Inline spacing (in units of CRS, e.g. meter)
- `--overwrite`: Overwrite output SEG-Y file(s) if existing
- `--verbose`, `-V`: Print additional information to stdout


