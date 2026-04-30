This directory contains the macros for matching the selections done previously using `run_filter_and_track.py` and the macros in `emu/`.


# Matching

- `match_xy_scifi_emu_cells.py`

  Main tool for Scifi/EMU matching. Does matching for all cells (or single row/column/cell using --cellx and/or --celly)

  It grabs all matched SciFi track projections (from hardcoded file containing TGraph objects (tgr_reco)). Then for each emulsion cell it grabs (from a hardcoded file) the background (TH2F, for aesthetics only, not used in matching) and a TGraph with the exit showers (tgr_emu).

  It removes the tgr_reco which fall outside the scope of the cell `is_sf_out_of_cell()`, then removes points (both emu and scifi)  which aren't close to any match. For each scifi point left, it performs a translation (x,y) to the nearest emu point (`move_to_closest()`). The average value of these translations for a single cell is a point in a histogram. The true translation for matching will show as a systematic bias from many individual cells.
  
- `mc_matching.py`

  It matches emulsion cluster showers to specific events in electronic detector cbmsim file. The cbmsim file contains the truth information, so this is useful for benchmarking the selections.

# Drawing macros

- `draw_mc_matching.py`

  Draws histograms from `mc_matching.py`

- `draw_xy_scifi_emu_cells.py`

  Draws roughly the same as `match_xy_scifi_emu_cells.py`, in case that was ran in no-drawing mode (--goff) and you want to draw again. But I'm actually not sure it's exactly the same drawing.

- `draw_xy_scifi_emu_reco.py`

  Draws TGraphs of scifi and emu for the entire brick. Likely broken.