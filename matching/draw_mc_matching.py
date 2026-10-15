import ROOT
ROOT.gROOT.SetBatch(True)
import rootUtils as ut
from array import array
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--input", dest="input", required=True, default="")
parser.add_argument("--sample",dest="sample",default="")
parser.add_argument("--binning",dest="binning",type=float,default=100) # binning in xy (um)
args = parser.parse_args()
#ROOT.gStyle.SetOptStat(0)
import math

fname_no_path = args.input.split("/")[-1]
tfile_sync = ROOT.TFile(args.input)
ttree_sync = tfile_sync.Get("mc_match")

x_em_points = [  -28.,  -7.  ]
y_em_points = [   15.,  35.  ]
z_em_points = [ 303.8, 311.6 ]

def emu2sfi(emux, emuy, emuz):
    x_t = -45.9
    y_t = 15.2
    z_t = 303.8
    sfix = emux/10000. + x_t
    sfiy = emuy/10000. + y_t
    sfiz = (emuz-1)*0.1350 + z_t
    return sfix, sfiy, sfiz

def sfi2emu(sfix, sfiy, sfiz):
    x_t = -45.9
    y_t = 15.2
    z_t = 303.8
    emux = (sfix - x_t)*10000.
    emuy = (sfiy - y_t)*10000.
    emuz = int(((sfiz - z_t)/0.1350) + 1)
    return emux, emuy, emuz

h1 = {}
ut.bookHist(h1, "h1_mcx",";Truth initial shower x (cm)", 100, -50, 50)
ut.bookHist(h1, "h1_mcy",";Truth initial shower y (cm)", 100, -50, 50)
ut.bookHist(h1, "h1_mcz",";Truth initial shower z (cm)", 100, 250, 350)
ut.bookHist(h1, "h1_clusx",";Cluster initial shower x (cm)", 100, -50, 50)
ut.bookHist(h1, "h1_clusy",";Cluster initial shower y (cm)", 100, -50, 50)
ut.bookHist(h1, "h1_clusz",";Cluster initial shower z (cm)", 100, 250, 350)
ut.bookHist(h1, "h1_purity",";Cluster purity", 100, 0.0, 1.0)

shw_start_z_min = 305
shw_start_z_max = 309

bins_ene = array('f',[0, 100, 200, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1600, 2000, 3000, 5000])
bins_ene_y = array('f', [0.0, 1.1])
h1["h2_eff_base"] = ROOT.TH2F("h1_eff_base",
                              ";Truth energy (GeV);Efficiency",
                              len(bins_ene)-1, bins_ene,
                              len(bins_ene_y)-1, bins_ene_y)
h1["h1_ene_reco"] = ROOT.TH1F("h1_ene_reco",
                              "Reconstructed clusters;Truth energy (GeV)",
                              len(bins_ene)-1, bins_ene)
h1["h1_ene_true"] = ROOT.TH1F("h1_ene_true",
                              "True showers;Truth energy (GeV)",
                              len(bins_ene)-1, bins_ene)
h1["h1_ene_reco1"] = ROOT.TH1F(f"h1_ene_reco1",
                               f"Reconstructed clusters, z>{shw_start_z_min};Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)
h1["h1_ene_true1"] = ROOT.TH1F(f"h1_ene_true1",
                               f"True showers, z>{shw_start_z_min};Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)
h1["h1_ene_reco2"] = ROOT.TH1F(f"h1_ene_reco2",
                               f"Reconstructed clusters, {shw_start_z_min}<z<{shw_start_z_max};Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)
h1["h1_ene_true2"] = ROOT.TH1F(f"h1_ene_true2",
                               f"True showers, {shw_start_z_min}<z<{shw_start_z_max};Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)
h1["h1_ene_reco3"] = ROOT.TH1F(f"h1_ene_reco3",
                               f"Reconstructed clusters, {shw_start_z_max}<z;Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)
h1["h1_ene_true3"] = ROOT.TH1F(f"h1_ene_true3",
                               f"True showers, {shw_start_z_max}<z;Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)
h1["h1_ene_reco4"] = ROOT.TH1F(f"h1_ene_reco4",
                               f"Reconstructed clusters, {shw_start_z_max}<z;Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)
h1["h1_ene_true4"] = ROOT.TH1F(f"h1_ene_true4",
                               f"True showers, {shw_start_z_max}<z;Truth energy (GeV)",
                               len(bins_ene)-1, bins_ene)


# cluster selection cuts
# look at energy spectrum (cut low energy)
l_cut_ntracks = [ 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200 ]
for ntrks in l_cut_ntracks:
    h1[f"h1_ene_reco_cutntracks_{ntrks}"] = ROOT.TH1F(f"h1_ene_reco_cutntracks_{ntrks}",
                                                      f"Reconstructed clusters, ntracks > {ntrks};Truth energy (GeV)",
                                                      len(bins_ene)-1, bins_ene)
    h1[f"h1_ene_true_cutntracks_{ntrks}"] = ROOT.TH1F(f"h1_ene_true_cutntracks_{ntrks}",
                                                      f"True showers, ntracks > {ntrks};Truth energy (GeV)",
                                                      len(bins_ene)-1, bins_ene)


h2 = {}
xend = -6
xini = -28
xbin = int((xend-xini)*10000. / args.binning)
ut.bookHist(h2, "h2_x",
            ";Initial #it{x} truth (cm);Initial #it{x} cluster reco (cm)",
            xbin, xini, xend,
            xbin, xini, xend)
yend = 35
yini = 15
ybin = int((yend-yini)*10000. / args.binning)
ut.bookHist(h2, "h2_y",
            ";Initial #it{y} truth (cm);Initial #it{y} cluster reco (cm)",
            ybin, yini, yend,
            ybin, yini, yend)
ut.bookHist(h2, "h2_z",
            ";Initial #it{z} truth (cm);Initial #it{z} cluster reco (cm)",
            100, 302.0, 312.0,
            100, 302.0, 312.0)
ut.bookHist(h2, "h2_zp",
            ";Initial #it{z} truth (cm);Initial #it{z} cluster reco (plate)",
            100, 302.0, 312.0,
            56, 1, 57)
ut.bookHist(h2, "h2_e",
            ";Truth shower energy (GeV);Number of tracks",
            500, 0.0, 5000.0,
            210, 0, 2100)
ut.bookHist(h2, "h2_e2",
            f"{shw_start_z_min} < True z < {shw_start_z_max};Truth shower energy (GeV);Number of tracks",
            500, 0.0, 5000.0,
            210, 0, 2100)

# cluster vs reco track
ut.bookHist(h2, "h2_x_reco",
            ";Scifi reco #it{x} (cm);Exit #it{x} cluster reco (cm)",
            xbin, xini, xend,
            xbin, xini, xend)
ut.bookHist(h2, "h2_y_reco",
            ";Scifi reco #it{y} (cm);Exit #it{y} cluster reco (cm)",
            ybin, yini, yend,
            ybin, yini, yend)
ut.bookHist(h2, "h2_z_reco",
            ";Scifi reco #it{z} (cm);Exit #it{z} cluster reco (cm)",
            100, 302.0, 312.0,
            100, 302.0, 312.0)

exit_plate = 54
# exiting cluster selection
ut.bookHist(h2, "h2_x_reco_exit",
            f"Cluster end > plate {exit_plate};Scifi reco #it{{x}} (cm);Exit #it{{x}} cluster reco (cm)",
            xbin, xini, xend,
            xbin, xini, xend)
ut.bookHist(h2, "h2_y_reco_exit",
            f"Cluster end > plate {exit_plate};Scifi reco #it{{y}} (cm);Exit #it{{y}} cluster reco (cm)",
            ybin, yini, yend,
            ybin, yini, yend)
ut.bookHist(h2, "h2_z_reco_exit",
            f"Cluster end > plate {exit_plate};Scifi reco #it{{z}} (cm);Exit #it{{z}} cluster reco (cm)",
            100, 302.0, 312.0,
            100, 302.0, 312.0)


h1["h_distx"] = ROOT.TH1F("h_distx",
                          ";Distance x emu - scifi (cm);Entries",
                          100,0,1000)
h1["h_disty"] = ROOT.TH1F("h_disty",
                          ";Distance y emu - scifi (cm);Entries",
                          100,0,1000)
h1["h_dist"] = ROOT.TH1F("h_dist",
                         ";Distance x emu - scifi (cm);Entries",
                         100,0,1000)

c = ROOT.TCanvas()
c.SetCanvasSize(1000,1000)
out_pdf = f"img/{fname_no_path}".replace(".root",".pdf")
c.Print(f"{out_pdf}[")



for ievent, event in enumerate(ttree_sync):
    if len(event.mc_shw_x) == 0:
        if len(event.shw_purity) > 0:
            h1["h1_purity"].Fill(event.shw_purity[0])
        continue
    
    h1["h1_mcx"].Fill(event.mc_shw_x[0])
    h1["h1_mcy"].Fill(event.mc_shw_y[0])
    h1["h1_mcz"].Fill(event.mc_shw_z[0])
    h1["h1_ene_true"].Fill(event.mc_shw_e[0])
    if event.mc_shw_z[0] > shw_start_z_min:
        h1["h1_ene_true1"].Fill(event.mc_shw_e[0])
    if event.mc_shw_z[0] > shw_start_z_min and event.mc_shw_z[0] < shw_start_z_max:
        h1["h1_ene_true2"].Fill(event.mc_shw_e[0])
    if event.mc_shw_z[0] > shw_start_z_max:
        h1["h1_ene_true3"].Fill(event.mc_shw_e[0])
    
    reco_cluster = False
    for clus in range(len(event.shw_z_ini)):
        emux, emuy, emuz = emu2sfi(event.shw_x_ini[clus], event.shw_y_ini[clus], event.shw_z_ini[clus])
        h1["h1_clusx"].Fill(emux)
        h1["h1_clusy"].Fill(emuy)
        h1["h1_clusz"].Fill(emuz)
        
        if not reco_cluster:
            h1["h1_purity"].Fill(event.shw_purity[0])
            h2["h2_e"].Fill(event.mc_shw_e[0],event.shw_ntracks[0])
            h1["h1_ene_reco"].Fill(event.mc_shw_e[0])
            if event.mc_shw_z[0] > shw_start_z_min:
                h1["h1_ene_reco1"].Fill(event.mc_shw_e[0])
            if event.mc_shw_z[0] > shw_start_z_min and event.mc_shw_z[0] < shw_start_z_max:
                h1["h1_ene_reco2"].Fill(event.mc_shw_e[0])
                h2["h2_e2"].Fill(event.mc_shw_e[0],event.shw_ntracks[0])
            if event.mc_shw_z[0] > shw_start_z_max:
                h1["h1_ene_reco3"].Fill(event.mc_shw_e[0])
            reco_cluster = True

        h2["h2_x"].Fill(event.mc_shw_x[0], emux)
        h2["h2_y"].Fill(event.mc_shw_y[0], emuy)
        h2["h2_z"].Fill(event.mc_shw_z[0], emuz)
        h2["h2_zp"].Fill(event.mc_shw_z[0], event.shw_z_ini[clus])

        emuxf, emuyf, emuzf = emu2sfi(event.shw_x_end[clus], event.shw_y_end[clus], event.shw_z_end[clus])
        h2["h2_x_reco"].Fill(event.reco_s2_x, emuxf)
        h2["h2_y_reco"].Fill(event.reco_s2_y, emuyf)
        h2["h2_z_reco"].Fill(event.reco_s2_z, emuzf)

        #if event.shw_z_end[clus] > exit_plate and event.shw_z_ini[clus] > shw_start_z_min and event.shw_z_ini[clus] < shw_start_z_max:
        if event.shw_z_end[clus] > exit_plate:
            h2["h2_x_reco_exit"].Fill(event.reco_s2_x, emuxf)
            h2["h2_y_reco_exit"].Fill(event.reco_s2_y, emuyf)
            h2["h2_z_reco_exit"].Fill(event.reco_s2_z, emuzf)

    # signal vs background
    # we define match distance dist = sqrt((xsfi-xemu)^2+(ysfi-yemu)^2)
    # signal: dist < value (e.g. 100um)
    # bkg1: dist > value (misreconstruction)
    # bkg2: no sfi reconstructed (orphan emu)
    # bkg3: no emu reconstructed (orphan sfi) <- will not appear if selection based on emu
    for ntrks in l_cut_ntracks:
        if len(event.shw_x_end) < 1:
            # fill bkg2
            continue
        distx = abs(event.reco_s2_x-event.shw_x_end[0])
        disty = abs(event.reco_s2_y-event.shw_y_end[0])
        dist = math.sqrt(pow(distx,2)+ pow(disty,2) )

        h1["h_distx"].Fill(distx)
        h1["h_disty"].Fill(disty)
        h1["h_dist"].Fill(dist)
        
        if event.shw_ntracks[0] > ntrks:
            #print("filling ntrks",ntrks,event.mc_shw_e[0])
            h1[f"h1_ene_reco_cutntracks_{ntrks}"].Fill(event.mc_shw_e[0])

for h in h1:
    if h == "h1_purity":
        ROOT.gStyle.SetOptStat(0)
    h1[h].Draw()
    c.Print(f"{out_pdf}")

teff = ROOT.TEfficiency(h1["h1_ene_reco"],h1["h1_ene_true"])
teff.Draw()
c.Print(f"{out_pdf}")
for i in range(1,4):
    teff = ROOT.TEfficiency(h1[f"h1_ene_reco{i}"],h1[f"h1_ene_true{i}"])
    if i == 2:
        sample = f"{args.sample}: " if args.sample != "" else ""
        teff.SetTitle(f"{sample}Truth shower start {shw_start_z_min} < z < {shw_start_z_max} (cm);Truth shower energy (GeV);Reconstruction efficiency")
    h1["h2_eff_base"].Draw()
    teff.Draw("same")
    c.Print(f"{out_pdf}")
for ntrks in l_cut_ntracks:
    teff = ROOT.TEfficiency(h1[f"h1_ene_reco_cutntracks_{ntrks}"],h1[f"h1_ene_true"])
    h1["h2_eff_base"].SetTitle(f"Reco cluster efficiency, ntracks > {ntrks};Truth shower energy (GeV);Reconstruction efficiency")
    h1["h2_eff_base"].Draw()
    teff.Draw("same")
    c.Print(f"{out_pdf}")
    
for h in h2:
    h2[h].Draw("colz")
    c.Print(f"{out_pdf}")
c.Print(f"{out_pdf}]")

nentries_xexit = h2["h2_x_reco_exit"].GetEntries()
print(f"h2_reco_exit_x entries: {nentries_xexit}")
