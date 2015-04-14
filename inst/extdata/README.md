## Units

The units of all coordinates are in microns (um).

## bcbp2

The BC/BCBP data.  Each data file contains three lines (x,y,c).  x and
y is the spatial coordinate and c determines the class of neuron (1
for BC, 2 for BCBP).  For each data file ending `.dat`, there is also
a file ending `.w` that contains four numbers (xl xh yl yh).  These
four numbers define the rectangular bounding box from (xl,yl) to (xh,
yh).  All points should be located within the bounding box.

## fth2

The ferret dopaminergic amacrine neurons data.  `f9942g.txt` contains
the coordinates of the GCL neurons, and `f9942i.txt` contains the
coordinates of the INL neurons.

## kram

The chick cone data.  This spreadsheet was taken from PLOS One and one
field (DN4) extracted.

## w81s1

The cat beta retinal ganglion cell data.  The on-centre neurons are
stored in `w81s1on.txt` and off-centre neurons in `w81s1of.txt`.

