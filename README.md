# LOVAMAP

## Overview

LOVAMAP (LOcal Void Analysis using Medial Axis by Particle configuration) is an analytical tool designed to characterize aspects of packed particles (i.e., assemblies, scaffolds, domains) - in particular, the void space between particles. It identifies and subcategorizes the medial axis of the void space and uses these elements to segment the space into 3D pores. LOVAMAP returns many descriptive metrics of the particle assembly, including features of the particles, 3D pores, paths through the void space, and more. In addition, data can be used to plot 3D pores and other features of the structure.

## Objective

LOVAMAP analyzes 3D particle domains and extracts detailed information about the void space between particles. The software follows a systematic multi-step process:

1. **Input Processing**: LOVAMAP accepts particle domain data in two formats:
   - **.txt/.csv/.dat files**: For particle assemblies comprising spheres, list each particle center coordinate (x, y, z) and radius. Begin the file with a 4-line description of the file. Particle data should begin on line 5.
   - **.json files**: For pre-voxelized domains, follow our JSON format (see 'JSON format' below).

2. **Methods** For a detailed explanation of LOVAMAP's methodology, please read the Methods section of [Riley et al. *Nat. Comp. Sci.* 2023](https://www.nature.com/articles/s43588-023-00551-x).

5. **Descriptors**: For a list and description of all metrics computed by LOVAMAP, see [lovamap.com/learn](https://lovamap.com/learn).

6. **Excel Outputs**: LOVAMAP includes an option to export all descriptors as an Excel file (see Excel file format).

### Hardware requirements

No special hardware is required to run `LOVAMAP`. MATLAB's minimum hardware requirements, found [here](https://www.mathworks.com/support/requirements/matlab-system-requirements.html), should be sufficient for data sets of modest size.

### Operating system requirements

### Software dependencies

The following MATLAB add-ons are necessary:
- [Image Processing Toolbox](https://www.mathworks.com/help/images/)
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/)
- [SLM - Shape Language Modeling](https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling)
- [circlefit3d - fit circle to three points in 3d space](https://www.mathworks.com/matlabcentral/fileexchange/34792-circlefit3d-fit-circle-to-three-points-in-3d-space)
- [Inhull](https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull)
- [Generate maximally perceptually-distinct colors](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
- [GetFullPath](https://www.mathworks.com/matlabcentral/fileexchange/28249-getfullpath)

If you plan on using `LOVAMAP` on a Windows machine **and** do not already have a C++-14 compliant compiler, you will also need the following MATLAB add-on:
- [MATLAB Support for MinGW-w64 C/C++ compiler](https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler)

## Installation guide

We include the script "setup_lovamap.m" to automatically download all external dependencies required to run LOVAMAP, i.e.,
- Ensure the appropriate MATLAB toolboxes are installed
- Install third-party dependencies
- Compile MEX files required to run LOVAMAP
- Add all dependency directories to the path

If issues occur while running this script, you may setup LOVAMAP manually following the instructions below.

### MATLAB add-ons

All of the add-ons listed above can be installed via MATLAB's add-on manager, which is accessible from MATLAB itself. The following are all released by MathWorks,
- [Image Processing Toolbox](https://www.mathworks.com/help/images/)
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/) 
- [MATLAB Support for MinGW-w64 C/C++ compiler](https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler) (only applicable to Windows machines)

so there will only be the option to install them. MATLAB will probably restart during each of these installations.

For the following, which are not released by MathWorks, select the option 'Add to MATALB (download and add to path)':
- [SLM - Shape Language Modeling](https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling)
- [circlefit3d - fit circle to three points in 3d space](https://www.mathworks.com/matlabcentral/fileexchange/34792-circlefit3d-fit-circle-to-three-points-in-3d-space)
- [Inhull](https://www.mathworks.com/matlabcentral/fileexchange/10226-inhull)
- [Generate maximally perceptually-distinct colors](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
- [GetFullPath](https://www.mathworks.com/matlabcentral/fileexchange/28249-getfullpath)

![Add to MATLAB (download and add to path)](docs/img/add_to_matlab.png)

### Compiling mex files

`LOVAMAP` relies on MEX files (located in `mex-lovamap`) that must be compiled first, which requires a C++-14 compliant compiler. Once that has been setup,

1. Navigate into the `mex-lovamap` directory in this repository.
2. Run the `compile_mex` script.
3. If successful, you should see `*.mex` file(s) appear in the `mex-lovamap` directory.

And you're done! Now, all methods that depend on C/C++ code are ready to go.


## Instructions for use
### Quick start

The easiest way to run LOVAMAP is to use the `runLOVAMAP.m` script, which serves as a template for configuring and executing an analysis. Here's the basic workflow:

1. **Place your input file** in the `particle_domains/` directory or adjust the `file_path` variable to point to your domain file.

2. **Configure analysis parameters** in `runLOVAMAP.m`:
   ```matlab
	% FILE PARAMETERS
	filename          = 'particle_assembly.txt';
	file_path         = './particle_domains/';
	excel_path        = './outputs/';

	% INPUT PARAMETERS
	%%%%%%%%%%%%%%%%% Only applicable to spherical data (.txt,.csv,.dat file inputs) %%%%%%%%%%%%%%%%%
	voxel_size        = 2;
	voxel_range       = [1e7, 1e8];
	%%%%%%%%%%%%%%%%%
	combine_edge_subs = true;

	% OUTPUT PARAMETERS
	generate_raw_data = true;
	interior_only     = true;
   ```

3. **Run the script**:
   ```matlab
   runLOVAMAP
   ```
   The analysis will execute and output results to `outputs/` directory.

### Understanding the parameters

| Parameter | Purpose | Guidance |
|-----------|---------|----------|
| `file_name` | Name of particle assembly file (.txt, .csv, .dat, .json) | Must be surrounded by single quotation marks. |
| `file_path` | Name of path containing your particle assembly file | Must be surrounded by single quotation marks. Must end with a forward slash. May begin with '/. notation if subfolders lie within current folder. |
| `excel_path` | Name of path/folder where LOVAMAP will deposit output files (.xlsx) | Must be surrounded by single quotation marks. Must end with a forward slash. May begin with '/. notation if output folder is to lie within current folder. |
| `voxel_size` | Target mesh size, i.e., side-length of cubic voxel, in micrometers | This parameter is only necessary when inputting spherical-particle data (.txt, .csv, .dat). If a labeled domain file (.json) is inputted, this parameter value will be ignored. LOVAMAP will select the final mesh size based on balancing the target `voxel_size` with the target `voxel_range`. Smaller `voxel_size` values generate finer resolution but increase computation time.|
| `voxel_range` | Target range of total number of voxels in the system | Range must be written within brackets and separated by a comma, e.g., `[min, max]`. This parameter is only necessary when inputting spherical-particle data (.txt, .csv, .dat). If a labeled domain file (.json) is inputted, this parameter value will be ignored. LOVAMAP will select the final mesh size based on balancing the target `voxel_range` with the target `voxel_size`. Larger `min` and `max` values of `voxel_range` generate more voxels (and therefore finer resolution) but increase computation time. |
| `combine_edge_subs` | Set to `true` to combine surface pores that share the same open space into and out of the packing | When set to `true`, LOVAMAP will combine 3D pores that extend to the same entrance (or exit) opening. Setting this parameter to `false` will likely lead to over-segmentation of pores at the surface of the packing, but this tradeoff may be necessary if the packing contains few particles. |
| `generate_raw_data` | Set to `true` to output Excel data | Setting this parameter to `false` will not inhibit LOVAMAP from running or generating the `data` variable in the Workspace; however, no Excel file will be deposited in the `excel_path` folder. |
| `interior_only` | Set to `true` to only output interior pore data for the 3D-pore portion of the output file | Setting this parameter to `false` will output descriptor data for all 3D pores, i.e., interior pores and pores that extend to the surface of the packing. The value of this parameter will not affect the output of the non-3D-pore descriptors. If `generate_raw_data` equals `false`, the value of this parameter will be ignored. |

### Input file formats

#### TXT/CSV/DAT format (for spherical particles)
Plain text file containing 4 commented-out or empty rows followed by rows with 4 columns, one particle per row, e.g.:
```
% Monodisperse Packing Assembly
% Generated on Jan. 12 2026
% Software by Houdini
% x_center   y_center   z_center   radius
100.5      250.3      150.7      50
102.1      252.8      152.2      50
...
```

#### JSON format (for voxelized data)
Structured data containing domain information and labeled particle regions. The JSON file must contain the following properties:

| Property | Description |
|----------|-------------|
| `domain_size` | The x, y, z dimensions of the domain in micrometers |
| `voxel_size` | The length of each cubic voxel in micrometers |
| `bead_count` | The total number of particles in the scaffold |
| `beads` | A dictionary where the key corresponds to the particle's randomly assigned numeric id (e.g. "1", "2", etc.) and the value is a list of inclusive 1D voxel index ranges that represent that particle's occupied voxels. |

```json
{
  "domain_size": [500, 500, 500],
  "voxel_size": 2.0,
  "bead_count": 125,
  "beads": {
    "1": [[1, 100], [200, 250]],
    "2": [[101, 199], [251, 350]],
    ...
  }
}
```

### Output files

LOVAMAP generates an Excel file (`.xlsx`) containing:
- **Summary statistics**: Overall void metrics, void fraction, ridge counts
- **Pore subunit data**: Individual rows for each identified pore with geometry and topology metrics
- **Ridge information**: Details on 1D and 2D ridge networks, connectivity
- **Time log**: Computation time for each processing step (useful for performance optimization)

The Excel filename is derived from your input filename. For example, `beadInfo_domain.dat` produces `stats_beadInfo_domain.xlsx`.

### Advanced usage: Direct function calls

For programmatic access, you can call LOVAMAP directly:

```matlab
[data, time_log] = LOVAMAP(domain_file, voxel_size, voxel_range, crop_percent, ...
    dip_percent, hall_cutoff, shell_thickness, num_2D_slices, combine_edge_subs);
```

**Outputs:**
- `data`: Structure containing all computed metrics (pore volumes, ridge lengths, void fraction, etc.)
- `time_log`: Struct array logging execution time for each processing step

Access specific results via `data.Subunits`, `data.Ridges1D`, `data.Ridges2D`, etc.

### Tips for best results

1. **Voxel size**: Choose a voxel size that captures relevant pore features without excessive memory usage. Test with smaller domains first.
2. **Hall cutoff**: Use values on the order of your particle radius to meaningfully segment pores. Experiment with ±20% variations to assess sensitivity.
3. **Cropping**: For large domains with boundary effects, set `crop_percent` to 0.8-0.9 to exclude edge artifacts.
4. **Ligand accessibility**: Adjust `shell_thickness` based on biological context (e.g., cell reach distance in tissue engineering applications).
5. **Monitoring**: Review `time_log` output to identify bottlenecks, especially for large datasets.