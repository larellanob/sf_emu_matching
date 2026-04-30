import ROOT
ROOT.gROOT.SetBatch(True)
import fedrarootlogon
import rootUtils as ut
from array import array
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
ROOT.gInterpreter.ProcessLine('#include "ScifiPoint.h"')
from ROOT.snd import analysis_tools as ana
ROOT.gStyle.SetOptStat(0)

from argparse import ArgumentParser

# first we draw two separate planes, to see if things check out, then
# we have to make a "get coordinate conversion" script


basedir = "/eos/experiment/sndlhc/users/luarella/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_18cmB21_transpVeto_SWANBINNING"
tfile_scifi      = ROOT.TFile(f"{basedir}/sndLHC.Ntuple-TGeant4_SFPremoved_digCPP_MuonTracks.root")
#tfile_emu        = ROOT.TFile(f"/afs/cern.ch/work/l/luarella/private/sf_emu_matching/emu/root/brick21_binning100um_top.root")
tfile_emu = ROOT.TFile("/afs/cern.ch/work/l/luarella/private/sf_emu_matching/emu/root/brick21_binning100um_plates_bintree_exitshw.root")


tree_cbm = tfile_scifi.Get("cbmsim")

# histograms
h2 = {}

# translations in x and y
x_t = -45.9
y_t = 15.2
nbinsx = 500
nbinsy = 500

xmin, xmax = 19, 39
ymin, ymax =  0, 20
ut.bookHist(h2,"emu_xy","FEDRA segments;x (cm);y (cm)",
            nbinsx,xmin, xmax,
            nbinsy,ymin, ymax)
ut.bookHist(h2,"emu_xy_t","FEDRA segments;x (cm);y (cm)",
            nbinsx,xmin+x_t, xmax+x_t,
            nbinsy,ymin+y_t, ymax+y_t)
ut.bookHist(h2,"emu_zy","FEDRA segments;z (plate);y (cm)",
            60,-60,0,
            nbinsy,ymin, ymax)
ut.bookHist(h2,"emu_zx","FEDRA segments;z (plate);x (cm)",
            60,-60,0,
            nbinsx,xmin, xmax)
lines_h = []
lines_v = []
for l in range(0,19):
    tlh = ROOT.TLine(20+x_t  , 0.45+l+y_t,   38+x_t, 0.45+l+y_t)
    tlv = ROOT.TLine(20+l+x_t, 0.45+y_t  , 20+l+x_t, 18.45+y_t)
    lines_h.append(tlh)
    lines_v.append(tlv)


sfx1, sfx2 = -27, -7
sfy1, sfy2 =  15, 35
sfnbinsx_fine = (sfx2-sfx1)*10 # 1mm bins
sfnbinsy_fine = (sfy2-sfy1)*10 # 1mm bins
sfnbinsx_coarse = (sfx2-sfx1)*2 # 5mm bins
sfnbinsy_coarse = (sfy2-sfy1)*2 # 5mm bins

h_scifi_base = ROOT.TH2F("h_scifi_base",
                         f";x (cm);y (cm)",
                         sfnbinsx_coarse, sfx1, sfx2, sfnbinsy_coarse, sfy1, sfy2)


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

looking_for = []
muon = []

dict_tgr = {}
for s in range(1,6):
    dict_tgr[f"tgr_scifi_true_{s}"] = ROOT.TGraph()
    dict_tgr[f"tgr_scifi_reco_{s}"] = ROOT.TGraph()
    
tgr_scifi_reco= {}
tgr_scifi_true= {}


#emulsions histograms
xmin, xmax = 19, 39
ymin, ymax =  0, 20
ut.bookHist(h2,"emu_xy","FEDRA segments;x (cm);y (cm)",
            nbinsx,xmin, xmax,
            nbinsy,ymin, ymax)
ut.bookHist(h2,"emu_xy_t","FEDRA segments;x (cm);y (cm)",
            nbinsx,xmin+x_t, xmax+x_t,
            nbinsy,ymin+y_t, ymax+y_t)
ut.bookHist(h2,"emu_zy","FEDRA segments;z (plate);y (cm)",
            60,-60,0,
            nbinsy,ymin, ymax)
ut.bookHist(h2,"emu_zx","FEDRA segments;z (plate);x (cm)",
            60,-60,0,
            nbinsx,xmin, xmax)


SZ = { 1:300.0, 2:313.0, 3:326.0, 4:339.0, 5:352.0 }
print("looping over scifi events and filling tgrs")
for ievent, event in enumerate(tree_cbm):

    if ievent %100 == 0: print(ievent)
    #if ievent > 100: break
    #SciFiDet.InitEvent(event.EventHeader)
    #MuFilter.InitEvent(event.EventHeader)
    # extra cuts on S2 or S3 hits here
    
    ttrk_min_x = 0.0
    ttrk_min_y = 0.0
    ttrk_min_z = 500.0
    seg_counter = 0
    for rtrk in event.Reco_MuonTracks:
        rtrk_pos = rtrk.getStart()
        rtrk_mom = rtrk.getTrackMom()
        rtrk_xyz = {}
        for s in range(1,6):
            # extrapolation of reco track to stations 1--5
            rtrk_xyz[s] = extrapolation(rtrk_pos,rtrk_mom, SZ[s])
            dict_tgr[f"tgr_scifi_reco_{s}"].AddPoint(rtrk_xyz[s].X(),rtrk_xyz[s].Y())
    for sfp in event.ScifiPoint:
        if sfp.station() != 2:
            continue
        if sfp.GetTrackID()==0:
            # we will fill the TGraph with the most upstream scifipoint (lowest Z)
            if sfp.GetZ() < ttrk_min_z:
                ttrk_min_x = sfp.GetX()
                ttrk_min_y = sfp.GetY()
                ttrk_min_z = sfp.GetZ()
    for s in range(1,6):
        dict_tgr[f"tgr_scifi_true_{s}"].AddPoint(ttrk_min_x,ttrk_min_y)
print("finished looping over events")

out_pdf_fname = f"img/xy_scifi_emu_cells.pdf"
c = ROOT.TCanvas()
c.SetCanvasSize(1000,1000)
c.Print(f"{out_pdf_fname}[")

h_scifi_base.Draw()
dict_tgr[f"tgr_scifi_true_2"].Draw("same P")
dict_tgr[f"tgr_scifi_true_2"].SetMarkerColor(ROOT.kGreen)
dict_tgr[f"tgr_scifi_true_2"].SetMarkerStyle(2)
dict_tgr[f"tgr_scifi_reco_2"].Draw("same P")
dict_tgr[f"tgr_scifi_reco_2"].SetMarkerColor(ROOT.kRed)
dict_tgr[f"tgr_scifi_reco_2"].SetMarkerStyle(5)
c.Print(f"{out_pdf_fname}")


for x in range(1,19):
    for y in range(1,19):
        print(f"Printing cell {x}, {y}")

        h2_cell = tfile_emu.Get(f"h2_{x}_{y}")
        tgr_emu = tfile_emu.Get(f"tgr_{x}_{y}")
        
        h2_cell.Draw("colz")
        tgr_emu.Draw("same P")
        # scifi
        dict_tgr[f"tgr_scifi_true_2"].Draw("same P")
        dict_tgr[f"tgr_scifi_reco_2"].Draw("same P")
        c.Print(f"{out_pdf_fname}")

c.Print(f"{out_pdf_fname}]")
