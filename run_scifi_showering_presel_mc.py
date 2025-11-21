import os
import ROOT
#ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gROOT.SetBatch(True)

ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
ROOT.gInterpreter.ProcessLine('#include "ScifiPoint.h"')
ROOT.gSystem.Load("libshipLHC.so")
from ROOT.snd import analysis_tools as ana

import rootUtils as ut

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--jobid", dest="jobid", required=True, type=int)
args = parser.parse_args()

with open("success_PassingMuSimEmu.3083711.txt") as txtf:
    l_jobid_good = [ int(line) for line in txtf ]

if args.jobid not in l_jobid_good:
    print("ERROR: jobid not in list of good jobs. Exiting")
    exit()
    
# tracking straight from MC file (no eventfilter)
input_filename = f"/eos/home-l/luarella/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_digi/{args.jobid+1}/sndLHC.Ntuple-TGeant4_digCPP_100000_MuonTracks.root"
geo_file = f"/eos/home-l/luarella/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_digi/{args.jobid+1}/geofile_full.Ntuple-TGeant4.root"
if args.jobid+1 != 300:
    print("ERROR: I have only processed the one file, use --jobid 299. Exiting.")
    exit()
    
S2z = 313.0
f = ROOT.TFile(input_filename)
print("Input filename:",input_filename)
tree = f.Get("cbmsim")

# scifi
SciFiDet, MuFilter = ana.GetGeometry(geo_file)
myopt = ROOT.snd.Configuration.Option.ti18_2022_2023
mySNDConfiguration = ROOT.snd.Configuration(myopt,
                                            SciFiDet,
                                            MuFilter)
# cut on SciFi station 2 (reco) hits
nhits_sf2_cut = 400

# histograms
h1 = {}
ut.bookHist(h1,"depo_energy_sh",
            f"Scifi 2, showers only;Energy deposited (GeV);Entries",
            100,0,1)
ut.bookHist(h1,"depo_energy_mu",
            f"Scifi 2, muon only;Energy deposited (GeV);Entries",
            100,0,1)


h2 = {}
ut.bookHist(h2,"rtrk_sfp_sh",
            f"S2, showers only;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-20,20,
            100,-20,20)
ut.bookHist(h2,"rtrk_sfp_mu",
            f"S2, from muon;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-20,20,
            100,-20,20)
ut.bookHist(h2,"rtrk_sfp_sh_zoom1",
            f"S2, showers only;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-2,2,
            100,-2,2)
ut.bookHist(h2,"rtrk_sfp_mu_zoom1",
            f"S2, from muon;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-2,2,
            100,-2,2)
ut.bookHist(h2,"rtrk_sfp_sh_zoom2",
            f"S2, showers only;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-0.2,0.2,
            100,-0.2,0.2)
ut.bookHist(h2,"rtrk_sfp_mu_zoom2",
            f"S2, from muon;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-0.2,0.2,
            100,-0.2,0.2)
ut.bookHist(h2,"rtrk_sfp",
            f"S2, from all;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-50,50,
            100,-50,50)
ut.bookHist(h2,"rtrk_sfp_sh_wt",
            f"S2, showers only, weighted;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-20,20,
            100,-20,20)
ut.bookHist(h2,"rtrk_sfp_sh_wt_zoom1",
            f"S2, showers only, weighted;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-2,2,
            100,-2,2)
ut.bookHist(h2,"rtrk_sfp_sh_wt_zoom2",
            f"S2, showers only, weighted;X_{{Reco track}} - X_{{Scifi point}} (cm);Y_{{Reco track}} - Y_{{Scifi point}} (cm)",
            100,-0.2,0.2,
            100,-0.2,0.2)

# sf centroid histograms
# reco track vs sf centroid
ut.bookHist(h2,"rtrk_sf",
            f"S2 reconstructed, weighted;X_{{Reco track}} - X_{{Scifi centroid}} (cm);Y_{{Reco track}} - Y_{{Scifi centroid}} (cm)",
            100,-20,20,
            100,-20,20)
ut.bookHist(h2,"rtrk_sf_zoom1",
            f"S2 reconstructed, weighted;X_{{Reco track}} - X_{{Scifi centroid}} (cm);Y_{{Reco track}} - Y_{{Scifi centroid}} (cm)",
            100,-2,2,
            100,-2,2)
ut.bookHist(h2,"rtrk_sf_zoom2",
            f"S2 reconstructed, weighted;X_{{Reco track}} - X_{{Scifi centroid}} (cm);Y_{{Reco track}} - Y_{{Scifi centroid}} (cm)",
            100,-0.2,0.2,
            100,-0.2,0.2)
# weighted sfpoints vs sf centroid 
ut.bookHist(h2,"sfp_sf",
            f"S2, weighted;X_{{Scifi point}} - X_{{Scifi centroid}} (cm);Y_{{Scifi point}} - Y_{{Scifi centroid}} (cm)",
            100,-20,20,
            100,-20,20)
ut.bookHist(h2,"sfp_sf_zoom1",
            f"S2, weighted;X_{{Scifi point}} - X_{{Scifi centroid}} (cm);Y_{{Scifi point}} - Y_{{Scifi centroid}} (cm)",
            100,-2,2,
            100,-2,2)
ut.bookHist(h2,"sfp_sf_zoom2",
            f"S2, weighted;X_{{Scifi point}} - X_{{Scifi centroid}} (cm);Y_{{Scifi point}} - Y_{{Scifi centroid}} (cm)",
            100,-0.2,0.2,
            100,-0.2,0.2)

if nhits_sf2_cut > 0:
    title_sf_cut = f" for events with N_{{hits}}^{{Scifi reco}} > {nhits_sf2_cut}"

ut.bookHist(h2,"toy_evd",
            f"Energy deposited in S2 (truth){title_sf_cut};x (cm);y (cm)",
            120,-40.9,-39.7,
            120, 19.2, 20.4)
ut.bookHist(h2,"toy_evd_mu",
            f"Energy deposited in S2 (truth){title_sf_cut}, from muon;x (cm);y (cm)",
            120,-40.9,-39.7,
            120, 19.2, 20.4)
ut.bookHist(h2,"toy_evd_sh",
            f"Energy deposited in S2 (truth){title_sf_cut}, from showers;x (cm);y (cm)",
            120,-40.9,-39.7,
            120, 19.2, 20.4)

ut.bookHist(h1,"mctrack_energy",
            f"Muon track energy;Energy (GeV);Entries",
            100,0,500)

def emu2scifi(xemu,yemu):
    xsf = (xemu - 473000)*1e-4
    ysf = (yemu + 158000)*1e-4
    return xsf, ysf

def scifi2emu(xsf,ysf):
    xem = (xsf*1E+4) + 473000
    yem = (ysf*1E+4) - 158000
    return xem, yem

def get_cell_borders(cell):
    width=10000
    extra_border = 1000
    cell_cx = (cell%18)*width+width
    cell_xi = cell_cx-width/2-extra_border
    cell_xf = cell_cx+width/2+extra_border
    cell_cy = (cell//18)*width+width
    cell_yi = cell_cy-width/2-extra_border
    cell_yf = cell_cy+width/2+extra_border
    return cell_xi, cell_xf, cell_yi, cell_yf

# emulsion cells
# we will move this to a later script, this takes too long and it's
# better to make some preselection first
cell_width=1.0 # all in cm
for cell in range(1,101):
    cell_xi, cell_xf, cell_yi, cell_yf = get_cell_borders(cell)
    # low x and y
    xsfi, ysfi = emu2scifi(cell_xi,cell_yi)
    # high x and y
    xsff, ysff = emu2scifi(cell_xf,cell_yf)

    ut.bookHist(h2,f"toy_evd_{cell}",
                f"Scifi S2 cell {cell} - truth dep. energy{title_sf_cut};x (cm);y (cm)",
                120,xsfi, xsff,
                120,ysfi, ysff)
    ut.bookHist(h2,f"toy_evd_mu_{cell}",
                f"Scifi S2 cell {cell} - truth dep. energy{title_sf_cut}, muon only;x (cm);y (cm)",
                120,xsfi, xsff,
                120,ysfi, ysff)
    ut.bookHist(h2,f"toy_evd_sh_{cell}",
                f"Scifi S2 cell {cell} - truth dep. energy{title_sf_cut}, shower only;x (cm);y (cm)",
                120,xsfi, xsff,
                120,ysfi, ysff)
    
# full canvas
xsfi, ysfi = emu2scifi(0.0,0.0)
#xsff, ysff = emu2scifi(180000,180000) # proper full canvas goes to 400000,400000
xsff, ysff = emu2scifi(400000,400000) # proper full canvas goes to 400000,400000
ut.bookHist(h2,f"full_canvas_em",
            f"Emulsion W2 - truth energy depositions{title_sf_cut};x (cm);y (cm)",
            2160,xsfi,xsff,
            2160,ysfi,ysff)
ut.bookHist(h2,f"full_canvas_sf",
            f"Scifi S2 - truth energy depositions{title_sf_cut};x (cm);y (cm)",
            2160,xsfi,xsff,
            2160,ysfi,ysff)

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

l_rtrks_x = []
l_rtrks_y = []

l_mctrks_x = []
l_mctrks_y = []

# out tree
make_output_root_file = False
if make_output_root_file:
    out_tree_filename = f"scifi_showering/emu_showering_sndLHC.Ntuple-TGeant4-160urad_100e6pp_FlukaEcut10_{args.jobid+1}_digCPP_100000_MuonTracks_presel.root"
    out_file = ROOT.TFile(out_tree_filename,"recreate")
    out_tree = tree.CloneTree(0)

# benchmarking
from time import process_time
time_start = process_time()
time_total = process_time()

# event loop
passing_events = 0
passing_fraction = 0
speed = 0
eta = 99999.
to_loop = tree.GetEntries()
total_pairs = 0

for ievent,event in enumerate(tree):
    if ievent > to_loop:
        break
    time_end = process_time()
    time_delta = time_end-time_start
    time_total += time_delta
    remaining = to_loop-ievent
    if ievent != 0:
        speed = ievent/float(time_total)
        eta = remaining/speed
        passing_fraction = float(passing_events)/float(ievent)
    estimate_final_passing = passing_fraction*to_loop
    print(f"Event {ievent}/{to_loop}. Passing events: {total_pairs}. Speed: {speed:.2f} (events/s). ETA: {eta/60.:.2f} (min)",end="\r")
    time_start = time_end

    SciFiDet.InitEvent(event.EventHeader)
    MuFilter.InitEvent(event.EventHeader)

    if len(event.Digi_ScifiHits) < 1:
        continue
    if len(event.Reco_MuonTracks) < 1:
        continue
    if len(event.EmulsionDetPoint) < 1:
        continue
    if len(event.ScifiPoint) < 1:
        continue

    # large number of scifihits on station 2
    h_hits = ROOT.std.vector('int')(5, 0)
    v_hits = ROOT.std.vector('int')(5, 0)
    ana.getSciFiHitsPerStation(event.Digi_ScifiHits,h_hits,v_hits)
    if h_hits[1]+v_hits[1] < nhits_sf2_cut: continue

    # emulsion points in wall 2
    # from muon
    emp_w2_mu = [ emp for emp in event.EmulsionDetPoint if ( emp.GetDetectorID()//10000==2 and emp.GetTrackID() == 0 ) ]
    # not from muon (i.e. from showers)
    emp_w2_sh = [ emp for emp in event.EmulsionDetPoint if ( emp.GetDetectorID()//10000==2 and emp.GetTrackID() != 0 ) ]

    # scifi points in station 2
    # from muon
    sfp_s2_mu = [ sfp for sfp in event.ScifiPoint if ( sfp.station() == 2 and sfp.GetTrackID() == 0 ) ]
    # not from muon (i.e. from showers)
    sfp_s2_sh = [ sfp for sfp in event.ScifiPoint if ( sfp.station() == 2 and sfp.GetTrackID() != 0 ) ]



    deposited_energy_s2_sh = 0.0
    deposited_energy_s2_mu = 0.0
    for sfp in sfp_s2_sh:
        deposited_energy_s2_sh += sfp.GetEnergyLoss()
    for sfp in sfp_s2_mu:
        deposited_energy_s2_mu += sfp.GetEnergyLoss()

    h1["depo_energy_sh"].Fill(deposited_energy_s2_sh)
    h1["depo_energy_mu"].Fill(deposited_energy_s2_mu)
    # reconstructed

    # initialise scifi plane objects, find, and get centroid
    v_ScifiStations = ana.FillScifi(mySNDConfiguration,
                                    event.Digi_ScifiHits,
                                    SciFiDet)
    v_ScifiStations[1].FindCentroid()
    sf2_xyz = v_ScifiStations[1].GetCentroid()

    # reco track at S2 vs deposited points x,y
    # "reco" track is just the MC track for now
    track_loop = event.Reco_MuonTracks
    for rtrk in track_loop:
        rtrk_pos = rtrk.getStart()
        rtrk_mom = rtrk.getTrackMom()
        total_pairs +=1 # pairs of reco track + reco digihits
        rtrk_xyz = extrapolation(rtrk_pos,rtrk_mom, S2z)
        l_rtrks_x.append(rtrk_xyz.X())
        l_rtrks_y.append(rtrk_xyz.Y())

        wt_x_mu = 0.0
        wt_y_mu = 0.0
        dep_energy_mu = 0.0
        nmu_depositions = 0

        wt_x_sh = 0.0
        wt_y_sh = 0.0
        dep_energy_sh = 0.0
        nsh_depositions = 0

        # fill scifi histograms and filter events
        for sfp in event.ScifiPoint:
            if sfp.station()!=2:
                continue
            # fill histos
            h2["rtrk_sfp"].Fill(rtrk_xyz.X()-sfp.GetX(),rtrk_xyz.Y()-sfp.GetY())
            if sfp.GetTrackID()==0:
                # from muon
                nmu_depositions += 1
                h2["rtrk_sfp_mu"].Fill(rtrk_xyz.X()-sfp.GetX(),rtrk_xyz.Y()-sfp.GetY())
                h2["rtrk_sfp_mu_zoom1"].Fill(rtrk_xyz.X()-sfp.GetX(),rtrk_xyz.Y()-sfp.GetY())
                h2["rtrk_sfp_mu_zoom2"].Fill(rtrk_xyz.X()-sfp.GetX(),rtrk_xyz.Y()-sfp.GetY())
                # x and y positiions weighted using deposited energy
                wt_x_mu += sfp.GetX()*sfp.GetEnergyLoss()
                wt_y_mu += sfp.GetY()*sfp.GetEnergyLoss()
                dep_energy_mu += sfp.GetEnergyLoss()
                h2["toy_evd_mu"].Fill(sfp.GetX(),sfp.GetY(),sfp.GetEnergyLoss())
                l_mctrks_x.append(sfp.GetX())
                l_mctrks_y.append(sfp.GetY())
            else:
                # not from muon (i.e. from showers)
                nsh_depositions += 1
                h2["rtrk_sfp_sh"].Fill(rtrk_xyz.X()-sfp.GetX(),rtrk_xyz.Y()-sfp.GetY())
                h2["rtrk_sfp_sh_zoom1"].Fill(rtrk_xyz.X()-sfp.GetX(),rtrk_xyz.Y()-sfp.GetY())
                h2["rtrk_sfp_sh_zoom2"].Fill(rtrk_xyz.X()-sfp.GetX(),rtrk_xyz.Y()-sfp.GetY())
                wt_x_sh += sfp.GetX()*sfp.GetEnergyLoss()
                wt_y_sh += sfp.GetY()*sfp.GetEnergyLoss()
                dep_energy_sh += sfp.GetEnergyLoss()
                h2["toy_evd_sh"].Fill(sfp.GetX(),sfp.GetY(),sfp.GetEnergyLoss())
                if sfp.GetEnergyLoss() == 0:
                    wt_x_sh += 1e-10
                    wt_y_sh += 1e-10
                    dep_energy_sh += 1e-10
                    #print("we have a passing sfp, not from the muon, which has energy 0.0")
            '''
            # we will do this in a different script, after preselection
            for cell in range(1,100):
                h2[f"toy_evd_{cell}"].Fill(sfp.GetX(),sfp.GetY(),sfp.GetEnergyLoss())
                if sfp.GetTrackID()==0:
                    h2[f"toy_evd_mu_{cell}"].Fill(sfp.GetX(),sfp.GetY(),sfp.GetEnergyLoss())
                else:
                    h2[f"toy_evd_sh_{cell}"].Fill(sfp.GetX(),sfp.GetY(),sfp.GetEnergyLoss())
            '''
            
            #h2["toy_evd"].Fill(sfp.GetX(),sfp.GetY(),sfp.GetEnergyLoss())
            h2["full_canvas_sf"].Fill(sfp.GetX(),sfp.GetY(),sfp.GetEnergyLoss())
        
        # end looping over scifipoints sfp

        # loop over emulsion points
        for emp in event.EmulsionDetPoint:
            if emp.GetDetectorID()//10000 != 2:
                continue
            if emp.GetEnergyLoss() == 0.0:
                continue
            h2["full_canvas_em"].Fill(emp.GetX(),emp.GetY(),emp.GetEnergyLoss())
        
        if nsh_depositions != 0:
            h2["rtrk_sfp_sh_wt"].Fill(rtrk_xyz.X()-wt_x_sh/dep_energy_sh, rtrk_xyz.Y()-wt_y_sh/dep_energy_sh)
            h2["rtrk_sfp_sh_wt_zoom1"].Fill(rtrk_xyz.X()-wt_x_sh/dep_energy_sh, rtrk_xyz.Y()-wt_y_sh/dep_energy_sh)
            h2["rtrk_sfp_sh_wt_zoom2"].Fill(rtrk_xyz.X()-wt_x_sh/dep_energy_sh, rtrk_xyz.Y()-wt_y_sh/dep_energy_sh)

        # reco track / sf centroid
        h2["rtrk_sf"].Fill(rtrk_xyz.X()-sf2_xyz.X(), rtrk_xyz.Y()-sf2_xyz.Y())
        h2["rtrk_sf_zoom1"].Fill(rtrk_xyz.X()-sf2_xyz.X(), rtrk_xyz.Y()-sf2_xyz.Y())
        h2["rtrk_sf_zoom2"].Fill(rtrk_xyz.X()-sf2_xyz.X(), rtrk_xyz.Y()-sf2_xyz.Y())

        wt_x = wt_x_mu + wt_x_sh
        wt_y = wt_y_mu + wt_y_sh
        dep_energy = dep_energy_mu + dep_energy_sh
        
        # weighted sfpoints / sf centroid
        h2["sfp_sf"].Fill(wt_x/dep_energy-sf2_xyz.X(), wt_y/dep_energy-sf2_xyz.Y())
        h2["sfp_sf_zoom1"].Fill(wt_x/dep_energy-sf2_xyz.X(), wt_y/dep_energy-sf2_xyz.Y())
        h2["sfp_sf_zoom2"].Fill(wt_x/dep_energy-sf2_xyz.X(), wt_y/dep_energy-sf2_xyz.Y())

print("")
print(f"Finished. Total pairs {total_pairs}")

c = ROOT.TCanvas()
h1["mctrack_energy"].Draw()
c.Print(f"scifi_showering/mctrack_energy_{args.jobid}.pdf")

if total_pairs == 0:
    print("no events passed selection!")
    exit()
    
label_jobid = f"_jobid{args.jobid}"
label_nhits_sf2_cut = f"_digihitsgt{nhits_sf2_cut}"
            
c = ROOT.TCanvas()
out_fig_filename = f"scifi_showering/EMU_mc_comparison_sf2gt{label_jobid}{label_nhits_sf2_cut}.pdf"

out_fig_evd_all_filename = f"scifi_showering/EMU_toy_evd_all_sf2gt{label_jobid}{label_nhits_sf2_cut}.pdf"
out_fig_evd_mu_filename = f"scifi_showering/EMU_toy_evd_mu_sf2gt{label_jobid}{label_nhits_sf2_cut}.pdf"
out_fig_evd_sh_filename = f"scifi_showering/EMU_toy_evd_sh_sf2gt{label_jobid}{label_nhits_sf2_cut}.pdf"

c.Print(f"{out_fig_filename}[")
#c.Print(f"{out_fig_evd_all_filename}[")
#c.Print(f"{out_fig_evd_mu_filename}[")
#c.Print(f"{out_fig_evd_sh_filename}[")
import array
#print(x_muon)
#print(y_muon)
for h in h1:
    h1[h].Draw()
    c.Print(f"{out_fig_filename}")
for h in h2:
    if "toy_evd" not in h:
        h2[h].Draw("colz")
        if "full_canvas" in h:
            ROOT.gStyle.SetOptStat(0)
            tg = ROOT.TGraph(len(l_rtrks_x),
                             array.array('d',l_rtrks_x),
                             array.array('d',l_rtrks_y))
            tg.SetMarkerSize(3)
            tg.SetMarkerStyle(52)
            tg.SetMarkerColor(ROOT.kRed)
            tg.Draw("P same")

            tgmc = ROOT.TGraph(len(l_mctrks_x),
                               array.array('d',l_mctrks_x),
                               array.array('d',l_mctrks_y))
            tgmc.SetMarkerSize(3)
            tgmc.SetMarkerStyle(50)
            tgmc.SetMarkerColor(ROOT.kGreen)
            tgmc.Draw("P same")

            leg = ROOT.TLegend(0.35,0.15,0.6,0.27)
            leg.AddEntry(tg,"Reco. scifi track","p")
            leg.AddEntry(tgmc,"Truth muon","p")
            leg.Draw("same")
            
        c.Print(f"{out_fig_filename}")
    else:
        continue
        h2[h].Draw("colz")
        if "mu" in h:
            c.Print(f"{out_fig_evd_mu_filename}")
        elif "sh" in h:
            c.Print(f"{out_fig_evd_sh_filename}")
        else:
            c.Print(f"{out_fig_evd_all_filename}")
c.Print(f"{out_fig_filename}]")
#c.Print(f"{out_fig_evd_all_filename}]")
#c.Print(f"{out_fig_evd_mu_filename}]")
#c.Print(f"{out_fig_evd_sh_filename}]")

print(f"Finished. Total pairs {total_pairs}")
