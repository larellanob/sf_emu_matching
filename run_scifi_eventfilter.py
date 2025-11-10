import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndTchainGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
from ROOT.snd import analysis_tools as ana

snd_root_path = os.environ["SNDSW_ROOT"]
run_paths = f"{snd_root_path}/analysis/tools/data_paths.csv"
geo_paths = f"{snd_root_path}/analysis/tools/geo_paths.csv"

import sys
ROOT.argv = sys.argv[:] # give ROOT its original list
# Now remove the flags that argparse shouldn't see.
clean_argv = [a for a in sys.argv if a != '-b']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--run" , dest="run"      , default=0    , type=int)
parser.add_argument("--mc"  , dest="mcfile"   , default="")
parser.add_argument("--frac", dest="run_frac" , default=1.0  , type=float)
args = parser.parse_args(clean_argv[1:])

if (args.run < 4575 or args.run > 4854) and args.mcfile == "":
    print("ERROR: Run number does not match emulsions Run 1 period and there is no MC file input.\n\tExiting.")
    exit()

import tools.run_database as rdb
# check our database for non-zero entries
# if not MC, use TchainGetter to get data file
if args.mcfile == "":
    ch = ana.GetTChain(run_paths,args.run)
    # geometry and detectors
    geo_file = ana.GetGeoPath(geo_paths,args.run)
    SciFiDet, MuFilter = ana.GetGeometry(geo_file)
    run_nevts = rdb.get_run_entries(args.run)
    print("run events",run_nevts)

else:
    mcpath = args.mcfile.split("/")[0:-1]
    mcpath = "/".join(mcpath)
    ch = ROOT.TChain("cbmsim")
    ch.Add(args.mcfile)
    ch2 = ROOT.TChain("cbmsim")
    ch2.Add(args.mcfile)
    geo_file = f"{mcpath}/geofile_full.Ntuple-TGeant4.root"
    SciFiDet, MuFilter = ana.GetGeometry(geo_file)
    run_nevts = ch2.GetEntries()
chfile = ch.GetFile()

# event loop
passing_events = 0
seen_events = 0

# output file
if args.mcfile == "":
    out_filename = f"scifi_eventfilter/scifi_selection_R{args.run}.root"
    if args.run_frac != 1.0:
        out_filename = f"scifi_eventfilter/scifi_selection_R{args.run}_frac{args.run_frac}.root"
else:
    mcfname = args.mcfile.split("/")[-1]
    mcfname = mcfname.replace(".root","")
    out_filename = f"scifi_eventfilter/scifi_selection_{mcfname}.root"
    if args.run_frac != 1.0:
        out_filename = f"scifi_eventfilter/scifi_selection_{mcfname}_frac{args.run_frac}.root"
        
out_file   = ROOT.TFile(out_filename,"recreate")
# output tree of rawConv, our filtered events
out_tree = ch.CloneTree(0)

# copy other stuff from source file
bl     = chfile.Get("BranchList")
tbl    = chfile.Get("TimeBasedBranchList")
ffh    = chfile.Get("FileHeader")
cbmout = chfile.Get("cbmout")
bl    .Write("BranchList", ROOT.TObject.kSingleKey)
tbl   .Write("TimeBasedBranchList", ROOT.TObject.kSingleKey)
ffh   .Write("FileHeader")
cbmout.Write("cbmout")

import tools.cuts as cuts

my_cutflow = cuts.Cutflow("Initial events:Both veto bars:Global number of hits:No zero hit stations:S2 hits > 100:Hits < 10 in S1, S3, or S4")
if args.mcfile == "":
    my_cutflow.set_data_name(f"Run {args.run}")
    if args.run_frac != 1.0:
        my_cutflow.set_data_name(f"Run {args.run} ({args.run_frac} fraction)")
else:
    my_cutflow.set_data_name(f"MC: {mcfname}")

for ievent, event in enumerate(ch):# event is an entry in the chain
    if ievent >= run_nevts*args.run_frac:
        print("Breaking event loop. Above looping fraction.")
        break

    SciFiDet.InitEvent(event.EventHeader)
    MuFilter.InitEvent(event.EventHeader)

    my_cutflow.Fill("Initial events")
    # veto cuts
    if not cuts.veto_cut(event.Digi_MuFilterHits):
        seen_events+=1
        continue
    my_cutflow.Fill("Both veto bars")

    if not cuts.scifi_gt_n_hits(event.Digi_ScifiHits,160):
        seen_events+=1
        continue
    my_cutflow.Fill("Global number of hits")

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
        seen_events+=1
        continue
    my_cutflow.Fill("No zero hit stations")

    # hardcoded cuts
    # this shuold be a proper Task at some point
    if hits[1] < 100:
        seen_events+=1
        continue
    my_cutflow.Fill("S2 hits > 100")
    if hits[0] > 10 or hits[3] > 10 or hits[4] > 10:
        seen_events+=1
        continue
    my_cutflow.Fill("Hits < 10 in S1, S3, or S4")
    # end of the cuts, event passed
    passing_events+=1
    pass_estimate = (float(passing_events)/seen_events)*run_nevts
    print(f"Last pass event: {ievent}/{run_nevts} ({100.*float(ievent)/run_nevts:.2f}% done). Pass so far: {passing_events}",
          end="\r")
    out_tree.Fill()

print("Writing tree")
out_file.Write()
out_file.Close()

h = my_cutflow.get_histogram()
c = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
histo_filename = out_filename.replace(".root","_cutflow.pdf")
c.Print(f"{histo_filename}[")
h.Draw()
c.Print(f"{histo_filename}")
c.SetLogy()
h.Draw()
c.Print(f"{histo_filename}")
c.Print(f"{histo_filename}]")
    
print("passing_events:",passing_events,"/",run_nevts)
if args.mcfile == "":
    print("Writing passing events (eventfilter) to run database")
    rdb.write_eventfilter(args.run,passing_events)
print("all events looped through vs tchain events:",seen_events,run_nevts)
