Simulation of pyroclastic avalnche performed over the 2014 Digital Elevation Model of Mount Etna (De Beni et al., 2015).
The original topography has been modifiedby excavating a detaching volume with an ellipsoidic shape
oriented towards the local maximum slope and with a prescribed total volume. The initial avalanche volume is defined by the
difference between the original and the modified topography.

The avalanche rheological parameters (Voellmy-Salm model) for this example have been selected in ranges consistent with previous studies of geophysical granular avalanches (e.g., Bartelt et al., 1999). In particular:

These values can be easily modified in the input file (IMEX_SloW2D.inp).

To run the simluation, first create a simbolic link of the executable in this folder:

>> ln -s ../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output to plot.
Usage example:

>> ./plot_p_2d.py example4_0075.p_2d

An additional Python script to plot the runout vs time is provided. The runout output file is written when the input flag OUTPUT_RUNOUT_FLAG is set to .TRUE. Usage example:

>> ./plot_runout.py example6_runout.txt

Finally, a Python script to create an animation of the simulation is also provided. This script plot the topography and animate the flow over it. An animated gif is saved by the script. The script requires two argument: run name (as given in IMEX_SfloW2D.inp) and duration of the animated gif in seconds

Usage example:

>> ./animate_asc.py example4 10


REFERENCES

De Beni, E., Behncke, B., Branca, S., Nicolosi, I., Carluccio, R., D’Ajello Caracciolo, F., and Chiappini, M.: The continuing story of Etna’s New Southeast Crater (2012–2014): Evolution and volume calculations based on field surveys and aerophotogrammetry, 303 (2015),
175–186, 2015.

Bartelt, P., Salm, L. B., and Gruberl, U.: Calculating dense-snow avalanche runout using a Voellmyfluid model with active/passive longitudinal straining, Journal of Glaciology, 45, 242–254, https://doi.org/10.3189/002214399793377301, 1999.
