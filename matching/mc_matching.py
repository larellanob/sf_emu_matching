import ROOT
ROOT.gROOT.SetBatch(True)
import rootUtils as ut
from array import array
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--emufile", dest="fname_emu", required=True, default="")
parser.add_argument("--sample", dest="sample", default="FLUKA")
args = parser.parse_args()


# mctrack / scifi file
if args.sample == "V2":
    fname_sfi = "/eos/experiment/sndlhc/users/luarella/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_18cmB21_transpVeto_SWANBINNING_V2/sndLHC.Ntuple-TGeant4_digCPP_shwEgt50.root"
elif args.sample == "FLUKA":
    fname_sfi = "/eos/experiment/sndlhc/users/luarella/PassingMu/ALL_lhc_ir1_coll_2023_exp001_fort_spectregen_SWANBINNING/sndLHC.Ntuple-TGeant4_digCPP_shwEgt50.root"
else:
    print("ERROR: unrecognized sample. Exiting")
    exit()
tfile_sfi = ROOT.TFile(fname_sfi)
ttree_sfi = tfile_sfi.Get("sfi_match")

# emu file
#fname_emu = f"/afs/cern.ch/work/l/luarella/private/sf_emu_matching/emu/root/clusters_bin100_minvox10_sorted.root"
tfile_emu = ROOT.TFile(args.fname_emu)
ttree_emu = tfile_emu.Get("emu_match")


#fname_sync = "/afs/cern.ch/work/l/luarella/private/sf_emu_matching/matching/root/sync.root"
fname_sync = args.fname_emu.replace("/emu/","/matching/")
fname_sync = fname_sync.replace(".root","_sync.root")
tfile_sync = ROOT.TFile(fname_sync,"recreate")
ttree_sync = ROOT.TTree("mc_match","Matching tree, from scifi cbmsim and emu fedra")

emu_x_ini   = ROOT.std.vector('float')()
emu_x_end   = ROOT.std.vector('float')()
emu_y_ini   = ROOT.std.vector('float')()
emu_y_end   = ROOT.std.vector('float')()
emu_z_ini   = ROOT.std.vector('int')()
emu_z_end   = ROOT.std.vector('int')()
emu_z_max   = ROOT.std.vector('int')()
emu_cell_x  = ROOT.std.vector('int')()
emu_cell_y  = ROOT.std.vector('int')()
emu_binning = ROOT.std.vector('int')()
emu_brick   = ROOT.std.vector('int')()
emu_nclusters = ROOT.std.vector('int')()
emu_ntracks = ROOT.std.vector('float')()
emu_purity  = ROOT.std.vector('float')()

emu_emcevt         = ROOT.std.vector(ROOT.std.vector('int'))()
emu_emcevt_ntracks = ROOT.std.vector(ROOT.std.vector('int'))()
clus_x = ROOT.std.vector(ROOT.std.vector('float'))()
clus_y = ROOT.std.vector(ROOT.std.vector('float'))()
clus_z = ROOT.std.vector(ROOT.std.vector('float'))()
clus_ntracks = ROOT.std.vector(ROOT.std.vector('float'))()

ttree_sync.Branch("shw_nclusters", emu_nclusters)
ttree_sync.Branch("shw_ntracks", emu_ntracks)
ttree_sync.Branch("shw_x_ini", emu_x_ini)
ttree_sync.Branch("shw_y_ini", emu_y_ini)
ttree_sync.Branch("shw_z_ini", emu_z_ini)
ttree_sync.Branch("shw_x_end", emu_x_end)
ttree_sync.Branch("shw_y_end", emu_y_end)
ttree_sync.Branch("shw_z_end", emu_z_end)
ttree_sync.Branch("shw_z_max", emu_z_max)
ttree_sync.Branch("shw_cell_x", emu_cell_x)
ttree_sync.Branch("shw_cell_y", emu_cell_y)
ttree_sync.Branch("shw_binning", emu_binning)
ttree_sync.Branch("shw_brick", emu_brick)
ttree_sync.Branch("shw_purity",emu_purity)
ttree_sync.Branch("shw_eMCEvt",emu_emcevt)
ttree_sync.Branch("shw_eMCEvt_ntracks",emu_emcevt_ntracks)
ttree_sync.Branch("clus_x",clus_x)
ttree_sync.Branch("clus_y",clus_y)
ttree_sync.Branch("clus_z",clus_z)
ttree_sync.Branch("clus_ntracks",clus_ntracks)



# sfi variables
mc_shw_n = array('i',[0])
#mc_shw_pdg = array('i',[0])
mc_shw_pdg = ROOT.std.vector('float')()
mc_shw_x = ROOT.std.vector('float')()
mc_shw_y = ROOT.std.vector('float')()
mc_shw_z = ROOT.std.vector('float')()
mc_shw_e = ROOT.std.vector('float')()
mc_shw_anglexz = ROOT.std.vector('float')()
mc_shw_angleyz = ROOT.std.vector('float')()
mc_shw_projs2_x = ROOT.std.vector('float')()
mc_shw_projs2_y = ROOT.std.vector('float')()
mc_shw_projs2_z = ROOT.std.vector('float')()
ttree_sync.Branch("mc_shw_n",mc_shw_n,"mc_shw_n/I")
ttree_sync.Branch("mc_shw_pdg",mc_shw_pdg)
ttree_sync.Branch("mc_shw_e",mc_shw_e)
ttree_sync.Branch("mc_shw_x",mc_shw_x)
ttree_sync.Branch("mc_shw_y",mc_shw_y)
ttree_sync.Branch("mc_shw_z",mc_shw_z)
ttree_sync.Branch("mc_shw_anglexz",mc_shw_anglexz)
ttree_sync.Branch("mc_shw_angleyz",mc_shw_angleyz)
ttree_sync.Branch("mc_shw_projs2_x",mc_shw_projs2_x)
ttree_sync.Branch("mc_shw_projs2_y",mc_shw_projs2_y)
ttree_sync.Branch("mc_shw_projs2_z",mc_shw_projs2_z)
# scifi hits
nhits_scifi = array('i',[-1,-1,-1,-1,-1])
ttree_sync.Branch("nhits_scifi",nhits_scifi,"nhits_scifi[5]/I")

# reco scifi track
SZ = { 1:300.0, 2:313.0, 3:326.0, 4:339.0, 5:352.0 }
reco_s2_x = array('f',[0.0])
reco_s2_y = array('f',[0.0])
reco_s2_z = array('f',[0.0])
reco_s2_anglexz = array('f',[0.0])
reco_s2_angleyz = array('f',[0.0])
reco_s2_chi2 = array('f',[0.0])
reco_s2_ndf = array('i',[0])
true_s2_x = array('f',[0.0])
true_s2_y = array('f',[0.0])
true_s2_z = array('f',[0.0])
ttree_sync.Branch("reco_s2_x",reco_s2_x,"reco_s2_x/F")
ttree_sync.Branch("reco_s2_y",reco_s2_y,"reco_s2_y/F")
ttree_sync.Branch("reco_s2_z",reco_s2_z,"reco_s2_z/F")
ttree_sync.Branch("reco_s2_anglexz",reco_s2_anglexz,"reco_s2_anglexz/F")
ttree_sync.Branch("reco_s2_angleyz",reco_s2_angleyz,"reco_s2_angleyz/F")
ttree_sync.Branch("reco_s2_chi2",reco_s2_chi2,"reco_s2_chi2/F")
ttree_sync.Branch("reco_s2_ndf",reco_s2_ndf,"reco_s2_ndf/I")
ttree_sync.Branch("true_s2_x",true_s2_x,"true_s2_x/F")
ttree_sync.Branch("true_s2_y",true_s2_y,"true_s2_y/F")
ttree_sync.Branch("true_s2_z",true_s2_z,"true_s2_z/F")
eMCEvt = array('i',[0])
ttree_sync.Branch("eMCEvt",eMCEvt,"eMCEvt/I")

def clear_emu():
    for vec in [ emu_x_ini, emu_x_end, emu_y_ini, emu_y_end, emu_z_ini, emu_z_end, emu_z_max, emu_cell_x, emu_cell_y, emu_binning, emu_brick, emu_nclusters, emu_ntracks, emu_purity ]:
        vec.clear()
    for vecvec in [ emu_emcevt, emu_emcevt_ntracks, clus_x, clus_y, clus_z, clus_ntracks]:
        vecvec.clear()
        for vec in vecvec:
            vec.clear()
                
def clear_sfi():
    for arr in [ mc_shw_n, reco_s2_x, reco_s2_y, reco_s2_z, reco_s2_anglexz, reco_s2_angleyz, reco_s2_chi2,  reco_s2_ndf, true_s2_x, true_s2_y, true_s2_z]:
        arr[0] = 0
    for vec in [ mc_shw_pdg, mc_shw_x, mc_shw_y, mc_shw_z, mc_shw_e, mc_shw_anglexz, mc_shw_angleyz, mc_shw_projs2_x, mc_shw_projs2_y, mc_shw_projs2_z ]:
        vec.clear()
    for s in range(0,5):
        nhits_scifi[s] = -1

def fill_sfi(event):
    mc_shw_n[0] = event.mc_shw_n
    for i in range(event.mc_shw_n):
        mc_shw_pdg.push_back(event.mc_shw_pdg[i])
        mc_shw_x.push_back(event.mc_shw_x[i])
        mc_shw_y.push_back(event.mc_shw_y[i])
        mc_shw_z.push_back(event.mc_shw_z[i])
        mc_shw_e.push_back(event.mc_shw_e[i])
        mc_shw_anglexz.push_back(event.mc_shw_anglexz[i])
        mc_shw_angleyz.push_back(event.mc_shw_angleyz[i])
        mc_shw_projs2_x.push_back(event.mc_shw_projs2_x[i])
        mc_shw_projs2_y.push_back(event.mc_shw_projs2_y[i])
        mc_shw_projs2_z.push_back(event.mc_shw_projs2_z[i])
    for s in range(5):
        nhits_scifi[s] = event.nhits_scifi[s]
    reco_s2_x[0] = event.reco_s2_x
    reco_s2_y[0] = event.reco_s2_y
    reco_s2_z[0] = event.reco_s2_z
    reco_s2_anglexz[0] = event.reco_s2_anglexz
    reco_s2_angleyz[0] = event.reco_s2_angleyz
    reco_s2_chi2[0] = event.reco_s2_chi2
    reco_s2_ndf[0] = event.reco_s2_ndf
    true_s2_x[0] = event.true_s2_x
    true_s2_y[0] = event.true_s2_y
    true_s2_z[0] = event.true_s2_z
    eMCEvt[0] = event.eMCEvt
    

def fill_emu(event):
    emu_x_ini.push_back(event.shw_x_ini)
    emu_y_ini.push_back(event.shw_y_ini)
    emu_z_ini.push_back(event.shw_z_ini)
    emu_x_end.push_back(event.shw_x_end)
    emu_y_end.push_back(event.shw_y_end)
    emu_z_end.push_back(event.shw_z_end)
    emu_z_max.push_back(event.shw_z_max)
    emu_cell_x.push_back(event.shw_cell_x)
    emu_cell_y.push_back(event.shw_cell_y)
    emu_binning.push_back(event.shw_binning)
    emu_brick.push_back(event.shw_brick)
    emu_nclusters.push_back(event.shw_nclusters)
    emu_ntracks.push_back(event.shw_ntracks)
    emu_purity.push_back(event.shw_purity)
    # vectors of vectors
    emu_emcevt.push_back(event.shw_eMCEvt)
    emu_emcevt_ntracks.push_back(event.shw_eMCEvt_ntracks)
    for i in range(event.shw_nclusters):
        clus_x.push_back(event.clus_x)
        clus_y.push_back(event.clus_y)
        clus_z.push_back(event.clus_z)
        clus_ntracks.push_back(event.clus_ntracks)


""" Event loop

We loop through scifi and emulsion using ievent and jevent
respectively
"""
jevent = 0
jentries = ttree_emu.GetEntries()
print("will loop through scifi tree")
for ievent, event in enumerate(ttree_sfi):
    if jevent >= jentries:
        # exhausted emulsion tree
        break
    mcevt_sfi = event.eMCEvt
    ttree_emu.GetEntry(jevent)
    #print(ievent,jevent)
    mcevt_emu = ttree_emu.shw_eMCEvt[0]

    # j > i:
    # scifi tree is ahead of emulsion tree
    # write i-0, advance i (i advances with main loop)
    if mcevt_emu > mcevt_sfi:
        clear_sfi()
        clear_emu()
        fill_sfi(event)
        ttree_sync.Fill()
        continue

    # j < i:
    # emulsion has event not in scifi
    # write 0-j, advance j
    if mcevt_emu < mcevt_sfi:
        clear_sfi()
        clear_emu()
        fill_emu(ttree_emu)
        jevent+=1
        ttree_emu.GetEntry(jevent)
        ttree_sync.Fill()
        continue
    
    # j = i:
    # mc match
    # write i-j, advance j
    filled_sfi = False
    while mcevt_emu == mcevt_sfi:
        if not filled_sfi:
            clear_sfi()
            clear_emu()
            fill_sfi(event)
            filled_sfi = True
        fill_emu(ttree_emu)
        jevent += 1
        if jevent >= jentries:
            break
        ttree_emu.GetEntry(jevent)
        mcevt_emu = ttree_emu.shw_eMCEvt[0]
    ttree_sync.Fill()
ttree_sync.Write()
