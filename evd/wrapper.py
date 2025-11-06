# wrapper made with extensive help from AI (lumo)

import sys
import importlib.util
import pathlib
fname = "/afs/cern.ch/work/l/luarella/private/sf_emu_matching/scifi_eventfilter/scifi_selection_R4619_ev18860895.root"
geofile = "/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V4_2022.root"

# dont need path because fname is absolute path already
path_to_file = ""

path_to_script = "/afs/cern.ch/user/l/luarella/sw/rhel9_x86-64/sndsw/master-local2/shipLHC/scripts/2dEventDisplay.py"

fake_argv = [
    path_to_script,          # pretend this is the program name
    "-f", fname,
    "-g", geofile,
    "-p", path_to_file,
    "--extension","pdf",
    "-hf", "linearSlopeIntercept",
    "--rootbatch"
]

old_argv = sys.argv[:]
sys.argv = fake_argv

import ROOT
ROOT.gROOT.SetBatch(True)

script_path = pathlib.Path(path_to_script).resolve()
spec = importlib.util.spec_from_file_location("_orig_mod", script_path)
orig_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orig_mod)   # runs the top‑level code of myscript.py

sys.argv = old_argv

# event numbers in the root file, not the SND@LHC event numbers
events_to_print = [18,20,188,842,44]
orig_mod.loopEvents(events_to_print,auto=True,hitColour="q",withHoughTrack=2)
#orig_mod.loopEvents(events_to_print,auto=True)
