run in this order:

1. `make_histograms.py`
	
	run this in condor preferably. inputs are: 
	- brick (21 for now)
	- cellx, celly (i usually run 18 jobs, one job per cellx, in which i run all 18 celly)
	- binwidth
	+- V2 (new simulation)+ (v2 is assumed now)

	the output is a root file containing 2d histograms for each plate,
    for the sum of all plates, and a "full histogram" (hence the macro
    name) that spans the whole brick. the idea is that each output
    cell file will only fill a small portion of the full histogram,
    and upon hadding them all you get the full brick.
	
	due to the large number (18x18=324) of root files, these are
    stored in the subdirectory root/cells(_V2)
	
	after all cells are finished, these are supposed to be hadded
	
2. run `hadd brick<brick>_binning<binning>um_plates.root root/brick<brick>_cell_*_binning<binning>um.root`

3. `make_bin_tree.py`
   
   this macro compiles all the information contained in the th2f into
   a ttree. i know, kind of redundand. why not make a ttree to begin
   with? I probably should, at some point. but at least this doesn't
   take that long, and is supposed to be ran on the hadded file. it
   brings the filesize from ~93MB to ~5MB. 
   
   output appends "_bintree" 
   
4. `make_selection_emu.py`

	This takes as input the bintree file makes the clusterization, and
    (for now) selects the last 7 plates.

	the output filename of this is a bit of a mess, right now it's
	`brick21_binning100um_plates_bintree_exitshw_voxthres4_clus0.root`
	- here exitshw is because these are exit showers
	- voxthres4 is the voxel threshold for it to be considered in
      clusterization
	- clus0 is a different threshold upon the full cluster, but i have
      not implemented it so that's why it's at zero
