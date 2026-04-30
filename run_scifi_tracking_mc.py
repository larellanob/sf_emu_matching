import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
from ROOT.snd import analysis_tools as ana

import sys
ROOT.argv = sys.argv[:] # give ROOT its original list
# Now remove the flags that argparse shouldn't see.
clean_argv = [a for a in sys.argv if a != '-b']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--input" , dest="file_in" , type=str  , required=True)
parser.add_argument("--output", dest="file_out", type=str                 )
parser.add_argument("--geo"   , dest="geofile" , type=str                 )
parser.add_argument("--nevents"  , dest="nevents", default=0, type=int  )
parser.add_argument("--startevent"  , dest="startevent", default=0, type=int  )
parser.add_argument("--workdir", dest="workdir", type=str, default="/afs/cern.ch/user/l/luarella/work-luarella/private/sf_emu_matching/")

args = parser.parse_args(clean_argv[1:])

# mcfile aliases
#base_dir = f"/eos/experiment/sndlhc/users/luarella/PassingMu/{pass_mu_sample}"

ch = ROOT.TChain("cbmsim")
input_filename = args.file_in
base_dir = "/".join(input_filename.split("/")[:-1])
if args.geofile == None:
    geo_file = f'{base_dir}/geofile_full.Ntuple-TGeant4.root'
else:
    geo_file = args.geofile
print("Geo file:",geo_file)
ch.Add(input_filename)
print(f"Input filename: {input_filename}")

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
ht_sf_task.SetParFile(f"/{args.workdir}/scifi_tracking/tracking_parameters/TrackingParams_mask_showers_scifi2.xml")
ht_sf_task.SetHoughSpaceFormat("linearSlopeIntercept")
ht_sf_task.SetTrackingCase('passing_mu_Sf')
ht_sf_task.SetStandalone()

fairrun.Init()

# output file
#out_tree_filename = f"{base_dir}/sndLHC.Ntuple-TGeant4_digCPP_MuonTracks.root"
out_tree_filename = input_filename.replace(".root","_MuonTracks.root")
to_loop = ioman.GetInTree().GetEntries()
if args.nevents != 0:
    to_loop = args.nevents
    if args.startevent != 0:
        to_loop = range(args.startevent,args.startevent+args.nevents)
    out_tree_filename = input_filename.replace(".root",f"_MuonTracks_ev{to_loop}.root")
print("OUT FILENAME",out_tree_filename)
out_file   = ROOT.TFile(out_tree_filename,"recreate")

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

    if ievent % (to_loop//100) == 0:
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
