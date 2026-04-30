import ROOT
ROOT.gROOT.SetBatch(True)
import fedrarootlogon
import rootUtils as ut
from array import array

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--input", dest="input", required=True)
parser.add_argument("--ntracks", dest="ntracks", type=int, default=20)
parser.add_argument("--minplate", dest="minplate", type=int, default=1)  # start plate
parser.add_argument("--maxplate", dest="maxplate", type=int, default=40) # start plate
parser.add_argument("--endplate", dest="endplate", type=int, default=54) # end plate
parser.add_argument("--inputscifi", dest="inputscifi", default="")
args = parser.parse_args()

tfile_emu = ROOT.TFile(args.input)
print("input file",tfile_emu)
ttree_emu = tfile_emu.Get("emu_match")

d_tgraphs = {}
for x in range(1,19):
    for y in range(1,19):
        print(f"tgr_{x}_{y}")
        tgr = ROOT.TGraph()
        tgr.SetName(f"tgr_{x}_{y}")
        d_tgraphs[f"tgr_{x}_{y}"] = tgr

print("will loop through events")
for ievent, event in enumerate(ttree_emu):
    # cut on ntracks
    if event.shw_ntracks < 20:
        continue

    # cut on plate number (shower start)
    if event.shw_z_ini < args.minplate:
        continue
    if event.shw_z_ini > args.maxplate:
        continue

    # cut on min end plate
    if event.shw_z_end < args.endplate:
        continue

    d_tgraphs[f"tgr_{event.shw_cell_x}_{event.shw_cell_y}"].AddPoint(event.shw_x_end, event.shw_y_end)

print("finished loop through events")
fname_in  = args.input.split("/")[-1]
fname_out = f"root/{fname_in}"
fname_out = fname_out.replace(".root",f"_ntrk{args.ntracks}_zini{args.minplate}-{args.maxplate}_zend{args.endplate}.root")
print("Filename out",fname_out)
tfile_out = ROOT.TFile(fname_out,"recreate")
for x in range(1,19):
    for y in range(1,19):
        d_tgraphs[f"tgr_{x}_{y}"].Write()

if args.inputscifi != "":
    fname_in  = args.inputscifi.split("/")[-1]
    dir_in = args.inputscifi.replace(f'{fname_in}',"")
    fname_out = fname_in.replace(".root","_tgr.root")
    fname_out = f"{dir_in}/{fname_out}")
    
    tfile_out = ROOT.TFile(fname_out,"recreate")
    
    
