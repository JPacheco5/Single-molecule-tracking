# Single-molecule-tracking
Single molecule tracking analysis in 1 and 2 colors
## Environment
The use of codes was standardized on:
- macOS Monterey v12.1
- ImageJ v2.0.0-rc-69/1.52p, TrackMate v6.0.3
- Python 3.9.16
- Spyder 5.1.5
- Anaconda 2.4.0

## Running time
The running time will depend on the number of trajectories per TrackMate file. For files ~10000 trajectories, each code must run in less than 5 min

## Tracking Files Names (1-color)
- Spots in tracks statistics.csv
- Track statistics.csv

## Tracking Files Names (2-color)
For 2-color tracking, the smaller numbers indicate the shorter wavelength:
- 1Spots in tracks statistics.csv
- 1Track statistics.csv
- 2Spots in tracks statistics.csv
- 2Track statistics.csv

## Area of Cell (µm²)
- Results.csv

## Binary Masks
Files saved with .tif and .txt extensions.
- 8-bit images in binary format with gray values of 0 and 255.

## Diffusion and Motion
### Use example dataset 2

### ALL-diffusion-alpha.py
1. Check the variable "files" to select the appropriate file directory.
2. Run the code.
3. The file "Diffusion-1-4" contains diffusion, alpha, R (correlation coefficient), and trajectory ID measurements.

### ALL-Radial-displacements.py
1. Check the variable "files" to select the appropriate file directory.
2. Run the code.
3. The file "radial resultsR" contains radial displacement measurements at 4 time-lags, which are used to generate cumulative distribution functions later fit with equations 3 and 4 to calculate diffusion and motion.

### Threshold-D.py
This code filters diffusion of PTH-TMR based on the threshold.
1. Check the variable "files" to select the appropriate file directory.
2. Ensure the input file is correctly named ('Diffusion-1-4.csv').
3. The file "Diffusion-1-4p” contains diffusion, alpha, R (correlation coefficient), and trajectory ID measurements with diffusion corrected values.

### ALL-diffusion-alpha-INT.py
Calculates the diffusion of colocalizing trajectories.
1. First, run the code "ALL-Free-interact.py".
2. Select the channel to analyze: XSpots in tracks statistics.csv and XTrack statistics.csv.
3. Generates a file with diffusion and motion of interactive trajectories.

## Colocalization Trajectories
### Use example dataset 1

### ALL-Free-interact.py
1. Check the variable "files" to select the appropriate file directory.
2. Generates the file "P_coloc" with colocalization data: Frames of interaction, Trajectories IDs for each channel, and coordinates of colocalization.
3. Generates files ResultsCH for each channel, quantifying the time pre-interaction and time of colocalization.

## Localization Codes
### Use example dataset 3

### localization_time1CH.py
This code quantifies the localization of molecules inside and outside domains.
- Requires a mask file saved as a .txt file with the coordinates of the binary mask.
- Requires the area of the cell saved in a file named Results.csv.

### Localization by intervals.py
Similar to the previous code, but for longer recordings, avoiding errors due to dynamic movement of the mask. Set for analysis every 30 seconds.

### Largest_empty_circles_voronoi.py
Calculates the average empty circle and the largest one.
1. Paste the binary mask.tif corresponding to the clathrin distribution on the Desktop. The image must have dimensions of 100x100px.
2. Generates a file named Radius, with all radii of circles and the image indicating the largest empty circle.

### molecules_per_frame.py
Counts the number of molecules, normalized by area and time.
1. Check the variable "files" to select the appropriate file directory.

## Clathrin Dynamics Codes
### Use example dataset 4

### CLC-Duration.py
Quantifies the duration of individual spots of clathrin.
1. Run Clathrin movies on TrackMate to generate trajectories in files Track statistics.csv and Spots in tracks statistics.csv.
2. Check the variable "files" to select the appropriate file directory.
3. Make histograms of trajectory duration.

### Initiation_rate.py
Calculates the number of Clathrin spots normalized by area and by minute.
1. Check the variable "files" to select the appropriate file directory.

