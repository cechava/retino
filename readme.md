# Retino package


The code contained in this repository implements the stimulation and analysis stages of a phase-encoded protocol to map out retinotopic organization across the visual cortex of rodents. This software was developed and used by Cesar Echavarria for use in his thesis project in the Cox Lab at Harvard University. 

The stimulation portion of the software presents visual stimuli and controls the acquisition of images of cortical tissue. 

\
 **Sample Cortical Tissue Image**
 \
 <img src="./sample_data/sample_frame.png" width = "218" height = "164">
 



 The analysis portion of the software takes in a stack of timestamped images of the cortex and outputs images reflecting how different regions of the screen map out across the exposed cortical tissue.

\
 **Sample Retinotopic Map**
 \
 <img src="./sample_output/maps/azimuth_stimulation_masked_phase_map.png" width = "218" height = "164">
 

Visual stimulation leads to a roughly 2% change in the pixel values of the cortical tissue image. This is too small to be appreciated by looking at the raw images. I normalize the change in pixel values and visualize the 'traveling waves' of neuronal activity to produces moves such as [**these ones**](./sample_output/movie).

## Materials

Visual stimuli were presented on a [72-inch LED TV](https://www.lg.com/us/tvs/lg-72LM9500-led-tv). Cortical tissue images were acquired with a [tandem-lens epifluoresence macroscope](https://www.sciencedirect.com/science/article/abs/pii/0165027091900382) coupled to a [CCD camera](https://www.alliedvision.com/en/products/cameras/detail/Manta/G-033.html)

## Intial Setup

1. create environments. Separate environment for acquisition (retino_acq) and analysis (retino_analysis). 

```
conda create env -f [filename].yml
```

* for some reason libtiff not properly installing from retino_acq.yml file.... will troubleshoot later


## Stimulation and Data Acquisition

The stimulation protocol consists of a bar moving across the screen in front of the animal in a periodic fashion. Sample stimulation movies can be found [**here**](./sample_stimulation⁩/movie/). A thread to simultaneously acquire images of the cortical tissue is initiated in the same sript.


```
python acquisition/Retinotopy_phaseEncoding_imageBar_constantImage.py -i [animal ID] -S [session] --save-images --output-path [path to folder]
```


## Analysis Workflow

**Note** These scripts are intended to be used for data collected within the Cox Lab so they expect a certain folder architecture. An update of the scripts for more general cases will be done later.

1. Analyze individual runs separately

```
python analysis/analyze_runs.py -i [animal ID] -S [session] -r [comma-separated list of runs] -m [Boolean for motion correction]
-n [Boolean for interpolation of data points to constant rate] -g [Boolean for removal of rolling mean] -w [integer indicating size of boxcar window for timecourse averaging of each pixel]
```

Typical example:
```
python analysis/analyze_runs.py -i JC026 -S 20181207 -r 'run1, run2, run3, run4, run5, run6' -m True -g True -w 11
```

This script analyzes individual runs separately. It also gets the image of the surface, performs a basic quality control check so the user can make sure that there wasn't an awful lot of movement between runs or shit data was acquired. The script also performs motion registration and correction, if indicated. Outputs unsmoothed maps for each run, which can be smoothed and thresholded with the script below.


2. Visualize single run results

```
python analysis/visualize_runs.py -i [animal ID] -S [session] -r [comma-separated list of runs] -m [Boolean for motion correction]
-n [Boolean for interpolation of data points to constant rate] -g [Boolean for removal of rolling mean] -w [integer indicating size of boxcar window for timecourse averaging of each pixel] -f [full-width at half-max size of kernel for smoothing] -t [magnitude ratio threshold]
```

Typical example:

```
python analysis/visualize_runs.py -i JC026 -S 20181207 -r 'run1, run2, run3, run4, run5, run6' -m True -g True -w 11 -f 7 -t .02
```

This script smooths phase map and threshold values based on magnitude ratio values per pixel for each run.

* At this point the user can look through individual run results and choose best runs for averaging in the subsequent steps. This tends to yield maps of better quality.


3. Average multiple runs and analyze

```
python analysis/average_and_analyze_runs.py -i [animal ID] -S [session] -r [comma-separated list of runs] -m [Boolean for motion correction]
-n [Boolean for interpolation of data points to constant rate] -g [Boolean for removal of rolling mean] -w [integer indicating size of boxcar window for timecourse averaging of each pixel]
```

Typical example:

```
python analysis/average_and_analyze_runs.py -i JC026 -S 20181207 -r 'run1, run2, run3, run4' -m True -g True -w 11
```

This script averages the timecourse of multiple runs and analyzes the magnitude and phase of the resulting timecourse.  Outputs unsmoothed maps for each condition, which can be smoothed and thresholded with the script below.

4. Visualize analysis results from multiple-run average

```
python analysis/visualize_average_run.py -i [animal ID] -S [session] -r [comma-separated list of runs] -m [Boolean for motion correction]
-n [Boolean for interpolation of data points to constant rate] -g [Boolean for removal of rolling mean] -w [integer indicating size of boxcar window for timecourse averaging of each pixel] -f [full-width at half-max size of kernel for smoothing] -t [magnitude ratio threshold]
```

Typical example:

```
python analysis/visualize_average_run.py -i JC026 -S 20181207 -r 'run1, run2, run3, run4' -m True -g True -w 11 -f 7 -t .02
```

This script smooths phase map and threshold values based on magnitude ratio values per pixel for each condition.

