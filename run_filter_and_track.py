import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndTchainGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
from ROOT.snd import analysis_tools as ana
import shipunit as u
from array import array
import rootUtils as ut
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--input" , dest="file_in" , type=str, required=True)
parser.add_argument("--output", dest="file_out", type=str)
parser.add_argument("--geo"   , dest="geofile" , type=str)
parser.add_argument("--mcfile", dest="mcfile"  , type=int)
parser.add_argument("--data"  , action='store_true')
parser.add_argument("--run"   , dest="run"     , type=int)
parser.add_argument("--nevents", dest="nevents", type=int, default=0     )
parser.add_argument("--workdir", dest="workdir", type=str, default="/afs/cern.ch/user/l/luarella/work-luarella/private/sf_emu_matching/")
parser.add_argument("--s1"    , dest="s1cut"   , type=int, default=120)
parser.add_argument("--s2"    , dest="s2cut"   , type=int, default=600)
parser.add_argument("--s3"    , dest="s3cut"   , type=int, default=120)
parser.add_argument("--oob"    , action='store_true') # remove out of bounds (brick21)
parser.add_argument("--trackpars"  , dest="trackpars", default="nomasking", type=str)

args = parser.parse_args()

d_pars = {
    "supertight": 2,
    "tight"     : 10,
    "medium"    : 30,
    "loose"     : 50,
    "veryloose" : 75,
    "superloose":100,
    "nomasking" :200000
    }

# mcfile aliases

SZ = { 1:300.0, 2:313.0, 3:326.0, 4:339.0, 5:352.0 }
def extrapolation(pos,mom,z):
    # pos is one known point
    # mom is the slope (we have two)
    slope_xz = mom.X()/mom.Z()
    slope_yz = mom.Y()/mom.Z()
    # unknown point (at z)
    unk_z = z
    unk_x = slope_xz*(unk_z-pos.Z())+pos.X()
    unk_y = slope_yz*(unk_z-pos.Z())+pos.Y()
    return ROOT.Math.XYZPoint(unk_x, unk_y, unk_z)

snd_root_path = os.environ["SNDSW_ROOT"]
# data / geo paths
run_paths = f"{snd_root_path}/analysis/tools/data_paths.csv"
geo_paths = f"{snd_root_path}/analysis/tools/geo_paths.csv"
input_filename = args.file_in
out_tree_filename = input_filename.replace(".root",f"_{args.s1cut}-s1_{args.s2cut}-s2_{args.s3cut}-s3_{args.trackpars}-pars_MuonTracks.root")
if not args.data:
    pass_mu_sample = "LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_18cmB21_transpVeto_SWANBINNING_V2"
    base_input_dir = f"/eos/experiment/sndlhc/users/luarella/PassingMu/{pass_mu_sample}"
    base_out_dir = f"/eos/experiment/sndlhc/users/luarella/PassingMu/{pass_mu_sample}"
    if args.geofile == None:
        geo_file = f'{base_input_dir}/{args.mcfile}/geofile_full.Ntuple-TGeant4.root'
    else:
        geo_file = args.geofile
    full_input_filename = f"{base_input_dir}/{args.mcfile}/{input_filename}"
    full_out_tree_filename = f"{base_out_dir}/{args.mcfile}/{out_tree_filename}"
    
if args.data:
    base_out_dir = "/eos/experiment/sndlhc/users/luarella/data/2022"
    base_input_dir = f"/eos/experiment/sndlhc/convertedData/physics/2022"
    if args.run == None:
        print("ERROR: is data but no run specified!")
        exit()
    geo_file = ana.GetGeoPath(args.run, geo_paths)
    full_input_filename = f"{base_input_dir}/run_{str(args.run).zfill(6)}/{input_filename}"
    full_out_tree_filename = f"{base_out_dir}/run_{str(args.run).zfill(6)}/{out_tree_filename}"

print("Geo file:",geo_file)


print("(Full path) input filename", full_input_filename)
print("(Full path) output filename", full_out_tree_filename)
ch = ROOT.TChain("cbmsim")
ch.Add(full_input_filename)

SciFiDet, MuFilter = ana.GetGeometry(geo_file)
myopt = ROOT.snd.Configuration.Option.ti18_2022_2023
mySNDConfiguration = ROOT.snd.Configuration(myopt,
                                            SciFiDet,
                                            MuFilter)
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
parfile = f"/{args.workdir}/scifi_tracking/tracking_parameters/TrackingParams_{args.trackpars}.xml"
print(f"Using parameter file: {parfile}")
ht_sf_task.SetParFile(f"{parfile}")
ht_sf_task.SetHoughSpaceFormat("linearSlopeIntercept")
ht_sf_task.SetTrackingCase('passing_mu_Sf')
ht_sf_task.SetStandalone()

fairrun.Init()

# output file
to_loop = ioman.GetInTree().GetEntries()
if args.nevents != 0:
    to_loop = args.nevents
print("OUT FILENAME",full_out_tree_filename)
out_file   = ROOT.TFile(full_out_tree_filename,"recreate")

# copy other stuff from source file
# actually dont think we need them
'''
bl     = f.Get("BranchList")
bl      .Write("BranchList", ROOT.TObject.kSingleKey)
tbl    = f.Get("TimeBasedBranchList")
tbl     .Write("TimeBasedBranchList", ROOT.TObject.kSingleKey)
ffh    = f.Get("FileHeader")
ffh     .Write("FileHeader")
cbmout = f.Get("cbmout")
cbmout  .Write("cbmout")
'''

# track container
myTracks = ht_sf_task.kalman_tracks

# output ttree. we add a branch of reco muon tracks.
cbmsim = ch.CloneTree(0)
cbmsim.Branch("Reco_MuonTracks",myTracks)

eMCEvt = array('i',[0])
s2_reco_x = array('f',[0.0])
s2_reco_y = array('f',[0.0])
s2_reco_z = array('f',[0.0])
s2_builtin_x = array('f',[0.0])
s2_builtin_y = array('f',[0.0])
s2_builtin_z = array('f',[0.0])
s2_true_x = array('f',[0.0])
s2_true_y = array('f',[0.0])
s2_true_z = array('f',[0.0])
nhits_scifi = array('i',[-1,-1,-1,-1,-1])
n_masked = array('i',[-1])
n_scifi_tracking = array('i', [-1,-1,-1,-1,-1])

cbmsim.Branch("eMCEvt",eMCEvt,"eMCEvt/I")
cbmsim.Branch("s2_reco_x",s2_reco_x,"s2_reco_x/F")
cbmsim.Branch("s2_reco_y",s2_reco_y,"s2_reco_y/F")
cbmsim.Branch("s2_reco_z",s2_reco_z,"s2_reco_z/F")
cbmsim.Branch("s2_builtin_x",s2_builtin_x,"s2_builtin_x/F")
cbmsim.Branch("s2_builtin_y",s2_builtin_y,"s2_builtin_y/F")
cbmsim.Branch("s2_builtin_z",s2_builtin_z,"s2_builtin_z/F")
cbmsim.Branch("s2_true_x",s2_true_x,"s2_true_x/F")
cbmsim.Branch("s2_true_y",s2_true_y,"s2_true_y/F")
cbmsim.Branch("s2_true_z",s2_true_z,"s2_true_z/F")
cbmsim.Branch("nhits_scifi",nhits_scifi,"nhits_scifi[5]/I")
cbmsim.Branch("n_masked",n_masked,"n_masked/I")
cbmsim.Branch("n_scifi_tracking",n_scifi_tracking,"n_scifi_tracking[5]/I")

# histograms
h1 = {}
ut.bookHist(h1,"h_energy_truth",
            f"args.mcfile {args.mcfile};Muon energy (GeV);All events",
            500, 0.0, 5000.0)
nbins_hits = 2000
ut.bookHist(h1,"h_nhits_s1",
            f"args.mcfile{args.mcfile};Nhits S1;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s2",
            f"args.mcfile{args.mcfile};Nhits S2;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s3",
            f"args.mcfile{args.mcfile};Nhits S3;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s4",
            f"args.mcfile{args.mcfile};Nhits S4;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s5",
            f"args.mcfile{args.mcfile};Nhits S5;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s12345",
            f"args.mcfile{args.mcfile};Nhits all stations;All events",
            nbins_hits*5, 0, nbins_hits*5)
ut.bookHist(h1,"h_nhits_s1_corr",
            f"args.mcfile{args.mcfile};Nhits S1;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s2_corr",
            f"args.mcfile{args.mcfile};Nhits S2;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s3_corr",
            f"args.mcfile{args.mcfile};Nhits S3;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s4_corr",
            f"args.mcfile{args.mcfile};Nhits S4;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s5_corr",
            f"args.mcfile{args.mcfile};Nhits S5;All events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s12345_corr",
            f"args.mcfile{args.mcfile};Nhits all stations;All events",
            nbins_hits*5, 0, nbins_hits*5)
# nhits in s2 cut
ut.bookHist(h1,"h_energy_truth_cut1",
            f"args.mcfile {args.mcfile}, N hits selection S2;Muon energy (GeV);Selected events",
            500, 0.0, 5000.0)
ut.bookHist(h1,"h_nhits_s1_cut1",
            f"N hits selection S2;Nhits S1;Selected events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s2_cut1",
            f"N hits selection S2;Nhits S2;Selected events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s3_cut1",
            f"N hits selection S2;Nhits S3;Selected events",
            nbins_hits, 0, nbins_hits)
# nhits in s3 cut
ut.bookHist(h1,"h_energy_truth_cut2",
            f"args.mcfile {args.mcfile}, N hits selection S3;Muon energy (GeV);Selected events",
            500, 0.0, 5000.0)
ut.bookHist(h1,"h_nhits_s1_cut2",
            f"N hits selection S2, S3;Nhits S1;Selected events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s2_cut2",
            f"N hits selection S2, S3;Nhits S2;Selected events",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s3_cut2",
            f"N hits selection S2, S3;Nhits S3;Selected events",
            nbins_hits, 0, nbins_hits)
# track reconstructed
ut.bookHist(h1,"h_energy_truth_cut3",
            f"args.mcfile {args.mcfile}, N hits selection S3;Energy (GeV);Selected tracks",
            500, 0.0, 5000.0)
ut.bookHist(h1,"h_angle_xz_cut3",
            f"N hits selection S2, S3;#theta_{{XZ}} (mrad);Reconstructed SciFi tracks",
            100 , -40, 40)
ut.bookHist(h1,"h_angle_yz_cut3",
            f"N hits selection S2, S3;#theta_{{YZ}} (mrad);Reconstructed SciFi tracks",
            100 , -40, 40)
ut.bookHist(h1,"h_nhits_s1_cut3",
            f"N hits selection S2, S3;Nhits S1;Reconstructed SciFi tracks",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s2_cut3",
            f"N hits selection S2, S3;Nhits S2;Reconstructed SciFi tracks",
            nbins_hits, 0, nbins_hits)
ut.bookHist(h1,"h_nhits_s3_cut3",
            f"N hits selection S2, S3;Nhits S3;Reconstructed SciFi tracks",
            nbins_hits, 0, nbins_hits)

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
for ievent, event in enumerate(ioman.GetInChain()):
    ch.GetEntry(ievent) # need to loop this one as well for the cloned tree
    if args.mcfile != None:
        eMCEvt[0] = 100000*(args.mcfile-1)+ievent
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

    if to_loop > 100:
        if ievent % (to_loop//100) == 0:
            # for condor
            print(f"Event {ievent}/{to_loop}. Tracks reco'd: {total_tracks}. Speed: {speed:.2f} (events/s). ETA: {eta/60.:.2f} (min)")
    time_start = time_end
    
    SciFiDet.InitEvent(event.EventHeader)
    MuFilter.InitEvent(event.EventHeader)

    # large number of scifihits on station 2
    h_hits = ROOT.std.vector('int')(5, 0)
    v_hits = ROOT.std.vector('int')(5, 0)
    ana.getSciFiHitsPerStation(event.Digi_ScifiHits,h_hits,v_hits)

    # the hits in h_hits and v_hits have no timing correction
    # will get time-corrected hits in h_hits_corr and v_hits_corr
    h_hits_corr = ROOT.std.vector('int')(5, 0)
    v_hits_corr = ROOT.std.vector('int')(5, 0)
    for s in range(1,6):
        if h_hits[s-1] != 0:
            h_f_hits = ana.selectScifiHits(event.Digi_ScifiHits,
                                           s,
                                           False,
                                           26,
                                           0.0,
                                           52.0,
                                           1E9/(2*u.snd_freq/u.hertz),
                                           2E9/(u.snd_freq/u.hertz),
                                           True,
                                           not args.data
                                           )
            h_hits_corr[s-1] = len(h_f_hits)
        else:
            h_hits_corr[s-1] = 0
        if v_hits[s-1] != 0:
            v_f_hits = ana.selectScifiHits(event.Digi_ScifiHits,
                                           s,
                                           True,
                                           26,
                                           0.0,
                                           52.0,
                                           1E9/(2*u.snd_freq/u.hertz),
                                           2E9/(u.snd_freq/u.hertz),
                                           True,
                                           not args.data
                                           )
            v_hits_corr[s-1] = len(v_f_hits)
        else:
            v_hits_corr[s-1] = 0

    muon_energy = 0.0
    for s in range(0,5):
        nhits_scifi[s] = -1
        n_scifi_tracking[s] = -1
    n_masked[0] = -1
    if not args.data:
        for mc in event.MCTrack:
            if mc.GetMotherId()==-1:
                muon_energy = mc.GetEnergy()
                break

    # print out large events
    if event.Digi_ScifiHits.GetEntries() > 1000:
        print("----------")
        print(f"EVENT {ievent}")
        print("scifihits entries",event.Digi_ScifiHits.GetEntries())
        for s in range(0,5):
            print(f"S    {s+1}: {h_hits[s]+v_hits[s]}")
            print(f"Scorr{s+1}: {h_hits_corr[s]+v_hits_corr[s]}")

    # start filling histograms
    h1["h_energy_truth"].Fill(muon_energy)
    sum_corr_hits = 0
    for s in range(1,6):
        h1[f"h_nhits_s{s}"].Fill(h_hits[s-1]+v_hits[s-1])
        h1[f"h_nhits_s{s}_corr"].Fill(h_hits_corr[s-1]+v_hits_corr[s-1])
        sum_corr_hits += h_hits_corr[s-1]+v_hits_corr[s-1]
    h1["h_nhits_s12345"].Fill(event.Digi_ScifiHits.GetEntries())
    h1["h_nhits_s12345_corr"].Fill(sum_corr_hits)

    # cuts on number of hits in s1 and s2
    if h_hits_corr[0]+v_hits_corr[0] > args.s1cut:
        continue
    if h_hits_corr[1]+v_hits_corr[1] < args.s2cut:
        continue
    h1["h_energy_truth_cut1"].Fill(muon_energy)
    h1["h_nhits_s1_cut1"].Fill(h_hits[0]+v_hits[0])
    h1["h_nhits_s2_cut1"].Fill(h_hits[1]+v_hits[1])
    h1["h_nhits_s3_cut1"].Fill(h_hits[2]+v_hits[2])

    # cut on number of hits in s3
    if h_hits_corr[2]+v_hits_corr[2] > args.s3cut:
        continue

    # passing events
    h1["h_energy_truth_cut2"].Fill(muon_energy)
    h1["h_nhits_s1_cut2"].Fill(h_hits[0]+v_hits[0])
    h1["h_nhits_s2_cut2"].Fill(h_hits[1]+v_hits[1])
    h1["h_nhits_s3_cut2"].Fill(h_hits[2]+v_hits[2])

    n_to_mask = d_pars[f"{args.trackpars}"]
    n_masked_ = 0
    for s in range(0,5):
        nhits_scifi[s] = h_hits[s]+v_hits[s]
        if h_hits[s]+v_hits[s] > n_to_mask:
            n_scifi_tracking[s] = 0
            n_masked_ +=1 
        else:
            n_scifi_tracking[s] = 1
    n_masked[0] = n_masked_
    
    # perform muon tracking on passing events
    ht_sf_task.kalman_tracks.Clear()
    ht_sf_task.Exec(0)
    myTracks = ht_sf_task.kalman_tracks

    s2_reco_x[0] = 0
    s2_reco_y[0] = 0
    s2_reco_z[0] = 0
    s2_builtin_x[0] = 0
    s2_builtin_y[0] = 0
    s2_builtin_z[0] = 0
    s2_true_x[0] = 0
    s2_true_y[0] = 0
    s2_true_z[0] = 0

    if len(myTracks) == 0:
        print(f"No track reconstructed!! (event {ievent})")
        continue
    else:
        total_tracks +=1
        h1["h_energy_truth_cut3"].Fill(muon_energy)
        h1["h_nhits_s1_cut3"].Fill(h_hits[0]+v_hits[0])
        h1["h_nhits_s2_cut3"].Fill(h_hits[1]+v_hits[1])
        h1["h_nhits_s3_cut3"].Fill(h_hits[2]+v_hits[2])
        h1["h_angle_xz_cut3"].Fill(myTracks[0].getAngleXZ()*1000.)
        h1["h_angle_yz_cut3"].Fill(myTracks[0].getAngleYZ()*1000.)

        # my own extrapolation
        rtrk_pos = myTracks[0].getStart()
        rtrk_mom = myTracks[0].getTrackMom()
        rtrk_xyz = extrapolation(rtrk_pos,rtrk_mom, SZ[2])
        s2_reco_x[0] = rtrk_xyz.X()
        s2_reco_y[0] = rtrk_xyz.Y()
        s2_reco_z[0] = rtrk_xyz.Z()

        # builtin extrapolation
        rtrk_xyz_builtin = ROOT.TVector3()
        rtrk_xyz_builtin = myTracks[0].extrapolateToPlaneAtZ(SZ[2])
        s2_builtin_x[0] = rtrk_xyz_builtin.X()
        s2_builtin_y[0] = rtrk_xyz_builtin.Y()
        s2_builtin_z[0] = rtrk_xyz_builtin.Z()

        # my own extrapolation and builtin extrapolation are very similar

    if args.oob:
        # these tracks fall outside brick 21
        # discard them using the --oob argument
        if s2_reco_x[0] < -25.9 or s2_reco_y[0] > 33.65:
            print(f"RECO track out of bounds event {ievent}")
            total_tracks -=1
            continue
        if s2_builtin_x[0] < -25.9 or s2_builtin_y[0] > 33.65:
            print(f"BUILTIN track out of bounds event {ievent}")
            total_tracks -=1
            continue
        if s2_true_x[0] < -25.9 or s2_true_y[0] > 33.65:
            print(f"TRUE track out of bounds event {ievent}")
            total_tracks -=1
            continue

    ttrk_min_x = 0.0
    ttrk_min_y = 0.0
    ttrk_min_z = 500.0
    
    # loop over scifipoints (sfp) when using MC
    # we will fill a TGraph with the most upstream muon sfp (lowest Z) in S2
    if not args.data:
        for sfp in event.ScifiPoint:
            if sfp.station() != 2:
                continue
            if sfp.GetTrackID()==0:
                if sfp.GetZ() < ttrk_min_z:
                    ttrk_min_x = sfp.GetX()
                    ttrk_min_y = sfp.GetY()
                    ttrk_min_z = sfp.GetZ()
        s2_true_x[0] = ttrk_min_x
        s2_true_y[0] = ttrk_min_y
        s2_true_z[0] = ttrk_min_z

    # fill the output ttree
    cbmsim.Fill()
    passing_events+=1
# write the output ttree
cbmsim .Write("cbmsim")
print("seen events",passing_events)
print("total tracks:",total_tracks)


# write histograms to _hist file
hist_filename = full_out_tree_filename.replace(".root","_hist.root")
hist_file = ROOT.TFile(hist_filename,"recreate")
for h in h1:
    h1[h].Write()

# we now loop over the tracks we just created to make the tgraph and save that too.
print("Now looping over muon tracks and projecting to S2")

tgr_reco = ROOT.TGraph()
tgr_builtin = ROOT.TGraph()
tgr_true = ROOT.TGraph()
tgr_reco.SetName("tgr_reco")
tgr_builtin.SetName("tgr_builtin")
tgr_true.SetName("tgr_true")

for ievent, event in enumerate(cbmsim):
    tgr_reco.AddPoint(event.s2_reco_x,event.s2_reco_y)
    tgr_builtin.AddPoint(event.s2_builtin_x,event.s2_builtin_y)
    tgr_true.AddPoint(event.s2_true_x,event.s2_true_y)

tgr_reco.Write()
tgr_builtin.Write()
tgr_true.Write()


hist_file.Close()
