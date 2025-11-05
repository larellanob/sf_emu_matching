import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
from ROOT.snd import analysis_tools as ana
import sys
ROOT.argv = sys.argv[:]                # give ROOT its original list
# Now remove the flags that argparse shouldn't see.
clean_argv = [a for a in sys.argv if a != '-b']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--rmin", dest="run_min", default=4575, type=int)
parser.add_argument("--rmax", dest="run_max", default=4575, type=int)
args = parser.parse_args(clean_argv[1:])

if args.run_max < args.run_min:
    args.run_max = args.run_min

# scifi
geo_file = "/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V4_2022.root"
if args.run_min > 5412:
    geo_file = "/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V3_2023.root"
SciFiDet, MuFilter = ROOT.snd.analysis_tools.GetGeometry(geo_file)

# 4575 - 4854
ch  = ROOT.TChain("rawConv")
ch2 = ROOT.TChain("rawConv")

#run_min = 4575
#run_min = 4619
#run_max = 4854
#run_max = 4619
nfiles = 100
for r in range(args.run_min,args.run_max+1):
    run_dir = f'/eos/experiment/sndlhc/convertedData/physics/2022/run_00{r}'
    if args.run_min > 5412:
        run_dir = f'/eos/experiment/sndlhc/convertedData/physics/2023/run_00{r}'
    if args.run_min == args.run_max:
        if not os.path.exists(run_dir):
            print("No directory for run",r,"EXIT")
            exit()
    if not os.path.exists(run_dir):
        continue
    for f in os.listdir(run_dir):
        current_file = run_dir+"/"+f
        ch.Add(current_file)
        ch2.Add(current_file)
        print("added file:",current_file)
chfile = ch.GetFile()

# this "goes to the end (??) of ch2, so we wasted it
# that's why we have ch and ch2
run_nevts = ch2.GetEntries()
if run_nevts == 0:
    print("WARNING: No events in run selection")
print("Added",run_nevts,"entries")

import cuts
# event loop
passing_events = 0
seen_events = 0
to_loop = run_nevts

# output file

if args.run_max > args.run_min:
    out_filename = f"scifi_eventfilter/scifi_selection_R{args.run_min}-{args.run_max}_ev{to_loop}.root"
elif args.run_max == args.run_min:
    out_filename = f"scifi_eventfilter/scifi_selection_R{args.run_min}_ev{to_loop}.root"
out_file   = ROOT.TFile(out_filename,"recreate")
# output tree of rawConv, our filtered events
out_raw = ch.CloneTree(0)


# copy other stuff from source file
bl     = chfile.Get("BranchList")
tbl    = chfile.Get("TimeBasedBranchList")
ffh    = chfile.Get("FileHeader")
cbmout = chfile.Get("cbmout")
bl    .Write("BranchList", ROOT.TObject.kSingleKey)
tbl   .Write("TimeBasedBranchList", ROOT.TObject.kSingleKey)
ffh   .Write("FileHeader")
cbmout.Write("cbmout")


for ievent, event in enumerate(ch):# event is an entry in the chain
    if ievent >= to_loop:
        print("breaking to loop")
        break

    SciFiDet.InitEvent(event.EventHeader)
    MuFilter.InitEvent(event.EventHeader)

    # veto cuts
    if not cuts.veto_cut(event.Digi_MuFilterHits):
        seen_events+=1
        continue
    elif not cuts.scifi_gt_n_hits(event.Digi_ScifiHits,160):
        seen_events+=1
        continue
    else:
        seen_events+=1
        passing_events+=1

    # grab horizontal and vertical hits
    h_hits = ROOT.std.vector('int')(5, 0)
    v_hits = ROOT.std.vector('int')(5, 0)
    ana.getSciFiHitsPerStation(event.Digi_ScifiHits,h_hits,v_hits)
    # will concatate them here
    hits=[]
    # hits in all scifi stations
    zero_hits_in_station = False
    for i in range(len(h_hits)):
        if h_hits[i] == 0 or v_hits[i] == 0:
            zero_hits_in_station = True
            break
        else:
            hits.append(h_hits[i]+v_hits[i])
    if zero_hits_in_station:
        passing_events-=1
        continue

    # hardcoded cuts
    # this shuold be a proper Task at some point
    if hits[1] < 100:
        passing_events-=1
        continue
    if hits[0] > 10 or hits[3] > 10 or hits[4] > 10:
        passing_events-=1
        continue

    pass_estimate = (float(passing_events)/seen_events)*run_nevts
    print(f"Last pass event: {ievent}/{run_nevts} ({100.*float(ievent)/run_nevts:.2f}% done). Pass so far: {passing_events}",
          end="\r")
    out_raw.Fill()

print("Writing rawConv")
out_file.Write()
out_file.Close()
    
print("passing_events:",passing_events,"/",run_nevts)
print("all events looped through vs tchain events:",seen_events,run_nevts)

