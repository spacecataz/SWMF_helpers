# SWMF_helpers
A set of scripts and tools for creating, executing, and analyzing SWMF simulations.

These scripts are useful bits of code for performing common tasks with input
and output of the Space Weather Modeling Framework.

Some of these are highly developed and powerful, others are quickly-assembled
shell scripts.  Read the individual help information for full descriptions
(use either the `-h` flag when calling or use `more` to read the file header).

Many rely on Spacepy, sometimes on *development versions* of Spacepy.
Users should be aware of this limitation when running into issues.

## Module Summary

SWMF_helpers includes a set of modules that support scripts.

| Module Name | Description |
| --------------------------|------------------------|
| sw_tools.py | Tools for fetching and processing of solar wind files. |

## Script Summary

| Script Name | Description |
| --------------------------|------------------------|
| calc_rim_I.py | Calculate integrated FAC from a series of RIM output files. |
| CatLog.py  | Carefully concatenate partial log files, virtual mag/sat files, and others into a single output file.  |
| cleanup.pl | Remove output, temporary, and restart files from an SWMF run directory. |
| concat_outs.py | Concatenate many `*.out` files into `*.outs` movie files. |
| convert_mags.py | Convert SWMF virtual magnetometer files into CCMC-like output files. |
| countfiles.py | Count files in the current working directory, list in order from most to fewest. |
| countfiles.sh | Count files in the current working directory, list in alphabetical order (faster than Python version). |
| create_iono_frames.py | Create a series of PNG files from a RIM simulation. |
| create_mhd_movieframes.py | Create a series of PNG files from a geospace simulation. |
| ezquota.py | Check file quota status on NASA's Pleiades supercomputer. |
| gen_ace_realtime.py | Create and continuously update a solar wind input file from real-time data. |
| l1_propagate.py | Ballistically propagate solar wind parameters from L1 to the SWMF upstream boundary. |
| mag_compare.py | Generate quick-look data-model comparisons for SWMF and SuperMag magnetometer data. |
| maggrid_extract.py | Extract a time series of values from a series of magnetometer grid output files. |
| make_movie.py | Using FFMPEG, turn a series of PNG files into a movie file. |
| plot_sim_summary.py | Create a summary plot for a geospace simulation. |
| process_sme.py | Create SME indexes from a set of magnetometer grid outputs. |
| rangescan.py | Determine the range of values inside of an SWMF output set. |
| RepairLog.py | Repair logfile-like files that have overlapping time entries. |
| runeff.py | Check the run efficiency of a simulation using its log file. |
| swmf_orbit_fetcy.py | Fetches satellite trajectories and creates SWMF virtual satellite input files. |
| unpack_outs.py | Unpack an `*.outs` movie file into a set of single-frame `*.out` files. |


## Dependencies
Dependencies are not consistent across all scripts.  For example,
*countfiles.sh* is a simple shell script and works on any BASH interface.
*create_mhd_movieframes.py*, however, requires Python 3 and related libraries
and the [Spacepy Python library](https://github.com/spacepy/spacepy
"Spacepy Repository).

This table quickly summarizes what is needed:
| Library/Software Name | Description |
| --------------------------|------------------------|
|Python 3.X  | Most scripts use Python 3.x |
|Numpy  1.16.X | Requirement for Spacepy |
|Matplotlib 3.1.X | All visualization done with MPL. |
|Scipy 1.3.X | Requirement for Spacepy |
|NASA CDF | >=3.6.X | Requirement for Spacepy |
|Spacepy >=0.2.3| Handles SWMF output, expedites visualization. |
|ffmpeg | Used to convert images into movie files. |
|hapiclient >= 0.2.6 | Used for fetching data from CDAweb, other sources. |
