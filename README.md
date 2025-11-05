Matching between SciFi and emulsion data.

Needs sndsw
1. run_scifi_eventfilter.py

   filters events according to scifi hits. selection is so far
   hardcoded, hits in all 5 stations, few hits in stations 1 4 and 5,
   large n of hits in station 2. hits in both veto bars.

2. run_scifi_tracking.py

   runs hough tracking algorithm over output file of previous
   step. digihits are passed through (rawConv), and makes new tree for
   tracks rather than adding to rawConv, want to change this
   later. input filename hardcoded and takes no arguments, will
   improve later.

3. run_scifi_showering.py

   runs shower reconstruction using analysis tools over output file of
   previous step. so far only plots, expects both hits and tracks for
   plotting shower-track variables such as the angle between them.
