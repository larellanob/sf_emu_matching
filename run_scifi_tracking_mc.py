import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
from ROOT.snd import analysis_tools as ana

snd_root_path = os.environ["SNDSW_ROOT"]
geo_paths = f"{snd_root_path}/analysis/tools/geo_paths.csv"

import sys
ROOT.argv = sys.argv[:] # give ROOT its original list
# Now remove the flags that argparse shouldn't see.
clean_argv = [a for a in sys.argv if a != '-b']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--run" , dest="run"      , default=0    , type=int)
parser.add_argument("--mc"  , dest="mcfile"   , default="2025")
parser.add_argument("--frac", dest="run_frac" , default=1.0  , type=float)
parser.add_argument("--mcgeofile", dest="mcgeofile" , default="")
parser.add_argument("--jobid", dest="jobid", required=True, type=int)

args = parser.parse_args(clean_argv[1:])

# mcfile aliases
mcdir_base = "/eos/user/l/luarella/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_digi"

with open("success_PassingMuSimEmu.3083711.txt") as txtf:
    l_jobid_good = [ int(line) for line in txtf ]
print(l_jobid_good)

# only digitized jobs above 300
# (energy correlates with jobid)
l_jobid_good = [ job for job in l_jobid_good if job >= 299 ]
print(l_jobid_good)

if args.jobid not in l_jobid_good:
    print("ERROR: jobid not in list of good jobs. Exiting")
    exit()

ch = ROOT.TChain("cbmsim")
mcpath = f"{mcdir_base}/{args.jobid+1}"
mcfilename = f"{mcpath}/sndLHC.Ntuple-TGeant4_digCPP.root"
geo_file   = f"{mcpath}/geofile_full.Ntuple-TGeant4.root"
ch.Add(mcfilename)
print(f"added file {mcfilename}")
    
SciFiDet, MuFilter = ana.GetGeometry(geo_file)
chfile = ch.GetFile()
print(chfile)

f = chfile
# tasks
ioman   = ROOT.FairRootManager.Instance()
ioman.SetTreeName("cbmsim")

fairrun     = ROOT.FairRunAna()
source  = ROOT.FairFileSource(f)
outFile = ROOT.TMemFile('dummy',"create")
sink    = ROOT.FairRootFileSink(outFile)

ioman.SetSource(source)
ioman.SetSink(sink)
#avoid some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

# Add Hough Transform tracking
import SndlhcMuonReco
ht_sf_task = SndlhcMuonReco.MuonReco()
fairrun.AddTask(ht_sf_task)
# non-default tracking parameters
# mask scifi stations with nhits > 10
ht_sf_task.SetParFile(f"/eos/home-l/luarella/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_digi/TrackingParams_mask_showers_scifi2.xml")
ht_sf_task.SetHoughSpaceFormat("linearSlopeIntercept")
ht_sf_task.SetTrackingCase('passing_mu_Sf')
ht_sf_task.SetStandalone()

fairrun.Init()

# output file
out_filename = f"{mcdir_base}/{args.jobid+1}/sndLHC.Ntuple-TGeant4_digCPP_MuonTracks.root"
if args.run_frac != 1.0:
    out_filename = f"{mcdir_base}/{args.jobid+1}/sndLHC.Ntuple-TGeant4_digCPP_MuonTracks_frac{args.run_frac}.root"

out_file   = ROOT.TFile(out_filename,"recreate")

# copy other stuff from source file
bl     = f.Get("BranchList")
bl      .Write("BranchList", ROOT.TObject.kSingleKey)
tbl    = f.Get("TimeBasedBranchList")
tbl     .Write("TimeBasedBranchList", ROOT.TObject.kSingleKey)
ffh    = f.Get("FileHeader")
ffh     .Write("FileHeader")
cbmout = f.Get("cbmout")
cbmout  .Write("cbmout")

# track container
myTracks = ht_sf_task.kalman_tracks

# output ttree. we add a branch of reco muon tracks.
cbmsim = ch.CloneTree(0)
cbmsim.Branch("Reco_MuonTracks",myTracks)

# benchmarking
from time import process_time
time_start = process_time()
time_total = process_time()

# event loop
passing_events = 0
passing_fraction = 0
speed = 0
eta = 99999.
total_tracks = 0
nentries = ioman.GetInTree().GetEntries()
to_loop = -1
if to_loop == -1:
    to_loop = nentries
import math
print("about to enter event loop")
for ievent, event in enumerate(ioman.GetInTree()):
    ch.GetEntry(ievent) # need to loop this one as well for the cloned tree
    if ievent > to_loop: break
    time_end = process_time()
    time_delta = time_end-time_start
    time_total += time_delta
    remaining = to_loop-ievent
    if ievent != 0:
        speed = ievent/float(time_total)
        eta = remaining/speed
        passing_fraction = float(passing_events)/float(ievent)
    estimate_final_passing = passing_fraction*to_loop

    if ievent % 5000 == 0:
        # for stdout
        #print(f"Event {ievent}/{to_loop}. Tracks reco'd: {total_tracks}. Speed: {speed:.2f} (events/s). ETA: {eta/60.:.2f} (min)",end="\r")
        # for condor
        print(f"Event {ievent}/{to_loop}. Tracks reco'd: {total_tracks}. Speed: {speed:.2f} (events/s). ETA: {eta/60.:.2f} (min)")
    time_start = time_end
    
    SciFiDet.InitEvent(event.EventHeader)
    MuFilter.InitEvent(event.EventHeader)

    # muon tracking
    ht_sf_task.kalman_tracks.Clear()
    ht_sf_task.Exec(0)
    myTracks = ht_sf_task.kalman_tracks
    for i in range(len(myTracks)):
        # keep count of tracks
        total_tracks +=1
    cbmsim.Fill()
    passing_events+=1
cbmsim .Write("cbmsim")
print("seen events",passing_events)
print("total tracks:",total_tracks)
