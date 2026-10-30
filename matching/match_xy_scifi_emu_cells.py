import ROOT
ROOT.gROOT.SetBatch(True)
import fedrarootlogon
import rootUtils as ut
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadRightMargin(0.145)
ROOT.gStyle.SetPadLeftMargin(0.125)
ROOT.gStyle.SetPadTopMargin(0.125)
ROOT.gStyle.SetTitleOffset(0.9,"z")
ROOT.gStyle.SetLegendBorderSize(0)

import utils.cell_geometry as cell_geom

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cellx", dest="cellx",default=0,type=int, help="If default, runs all X in a loop")
parser.add_argument("--celly", dest="celly",default=0,type=int, help="If default, runs all Y in a loop")
parser.add_argument("--tx", dest="tx",default=0,type=float)
parser.add_argument("--ty", dest="ty",default=0,type=float)
parser.add_argument("--thres", dest="thres",default=0.06,type=float)
parser.add_argument("--goff", action='store_true')
parser.add_argument("--input_emu", required=True)
parser.add_argument("--input_scifi", required=True)
parser.add_argument("--data", action='store_true')
parser.add_argument("--brick", type=int, required=True)

args = parser.parse_args()

# input files
print("Scifi file\n",args.input_scifi)
tfile_scifi = ROOT.TFile(args.input_scifi)
    
print("Emu file\n",args.input_emu)
tfile_emu = ROOT.TFile(args.input_emu)

# capture TGraphs in files
tgr_sfi = tfile_scifi.Get("tgr_reco")
tgr_sfi.SetMarkerColor(ROOT.kRed)
tgr_sfi.SetMarkerStyle(52)
tgr_sfi.SetMarkerSize(3)
tgr_sfi_truth = tfile_scifi.Get("tgr_true")
tgr_sfi_truth.SetMarkerColor(ROOT.kGreen)
tgr_sfi_truth.SetMarkerStyle(50)
tgr_sfi_truth.SetMarkerSize(3)

def translate_scifi_inplace(tgr, tx, ty):
    # modifies the original TGraph in place
    for p in range(tgr.GetN()):
        tgr.SetPoint(p, tgr.GetPointX(p) + tx, tgr.GetPointY(p) + ty)

translate_scifi_inplace(tgr_sfi, args.tx, args.ty)

from math import sqrt

def dist(x1, y1, x2, y2):
    return sqrt( (x2-x1)**2 + (y2-y1)**2)

def any_nearby(x1, y1, tgr2, thres = 0.06):
    # compare tgr1 point with all points of tgr2
    # if any closer than thres, return true, else discard
    for x2, y2 in zip(tgr2.GetX(), tgr2.GetY()):
        if dist(x1,y1,x2,y2) < thres:
            return True
    return False

def is_sf_out_of_cell(cellx,celly, sfx, sfy,brick=21):
    cxi, cxf, cyi, cyf = cell_geom.get_cell_borders(cellx, celly, brick)
    cxi, cyi = cell_geom.emu_to_sf(cxi, cyi, brick)
    cxf, cyf = cell_geom.emu_to_sf(cxf, cyf, brick)
    
    if sfx < cxi or sfx > cxf:
        return True
    elif sfy < cyi or sfy > cyf:
        return True
    else:
        return False

def move_to_closest(x,y, tgr2):
    closest = -1
    min_dist = 99999.
    movex = 0.0
    movey = 0.0
    for j, (x2,y2) in enumerate(zip(tgr2.GetX(), tgr2.GetY())):
        d = dist(x,y,x2,y2)
        if d < min_dist:
            closest = j
            min_dist = d
            movex = x2-x
            movey = y2-y
    return closest, movex, movey, min_dist

c = ROOT.TCanvas()
c.SetCanvasSize(1100,1000)
out_pdf = f"img/match_xy_scifi_emu_cells_brick{args.brick}_th{args.thres}_tx{args.tx}_ty{args.ty}_s2only.pdf"
if args.data:
    out_pdf = out_pdf.replace(".pdf","_data.pdf")
if not args.goff:
    c.Print(f"{out_pdf}[")
h_global_x = ROOT.TH1F("h_global_x", ";x translation (cm);Cells", 50, -0.1, 0.1)
h_global_y = ROOT.TH1F("h_global_y", ";y translation (cm);Cells", 50, -0.1, 0.1)

h_global_x.GetXaxis().SetNdivisions(505)
h_global_y.GetXaxis().SetNdivisions(505)

h_global_x_avg = ROOT.TH1F("h_global_x_avg", ";x translation (cm) (per-cell average);Cells", 50, -0.1, 0.1)
h_global_y_avg = ROOT.TH1F("h_global_y_avg", ";y translation (cm) (per-cell average);Cells", 50, -0.1, 0.1)

h_global_x_avg.GetXaxis().SetNdivisions(505)
h_global_y_avg.GetXaxis().SetNdivisions(505)

h_initial_emu = ROOT.TH1F("h_initial_emu",";Number of initial emu points;Cells",30,0,30)
h_dropped_emu = ROOT.TH1F("h_dropped_emu",";Number of dropped emu points;Cells",20,0,20)
h_dropped_sfi = ROOT.TH1F("h_dropped_sfi",";Number of dropped scifi points;Cells",50,0,50)
h_close_matches = ROOT.TH1F("h_close_matches",";Number of close matches per cell;Cells",20,0,10)

for cellx in range(1,19):
    if args.cellx != 0 and cellx != args.cellx:
        continue
    for celly in range(1,19):
        if args.celly != 0 and celly != args.celly:
            continue
        h_cell = tfile_emu.Get(f"h2_{cellx}_{celly}")
        if h_cell == None:
            print(f"Could not find cell {cellx}, {celly}!")
            continue
        print(f"Cell {cellx}, {celly}")
        h_cell.GetXaxis().SetNdivisions(505)
        h_cell.GetZaxis().SetTitleOffset(1.35)
        tgr_emu = tfile_emu.Get(f"tgr_{cellx}_{celly}")

        #tgr_emu = translate_emulsion(tgr_emu, args.emutx, args.emuty)
        
        h_initial_emu.Fill(tgr_emu.GetN())
        h_cell.SetTitle(f"Cell {cellx}, {celly};x (cm);y (cm);First segment count")
        h_cell.Draw("colz")
        tgr_emu.Draw("same P")
        tgr_sfi.Draw("same P")
        if not args.data:
            tgr_sfi_truth.Draw("same P")
        myleg = ROOT.TLegend(0.25,0.89,0.8,0.94)
        myleg.SetNColumns(3)
        myleg.AddEntry(tgr_emu,"Emu exit showers","P")
        myleg.AddEntry(tgr_sfi,"Scifi reco track","P")
        if not args.data:
            myleg.AddEntry(tgr_sfi_truth,"Scifi truth muon","P")
        myleg.Draw("same P")
        if not args.goff:
            c.Print(f"{out_pdf}")
        
        # discard non-easily paired SF and EMU
        bad_sfi = []
        bad_emu = []
        good_sfi = []
        good_emu = []
        tgr_good_emu = ROOT.TGraph()
        tgr_good_emu.SetMarkerColor(tgr_emu.GetMarkerColor())
        tgr_good_emu.SetMarkerStyle(tgr_emu.GetMarkerStyle())
        tgr_good_emu.SetMarkerSize(tgr_emu.GetMarkerSize())
        tgr_good_sfi = ROOT.TGraph()
        tgr_good_sfi.SetMarkerColor(tgr_sfi.GetMarkerColor())
        tgr_good_sfi.SetMarkerStyle(tgr_sfi.GetMarkerStyle())
        tgr_good_sfi.SetMarkerSize(tgr_sfi.GetMarkerSize())
        for isfi, (xsfi,ysfi) in enumerate(zip(tgr_sfi.GetX(), tgr_sfi.GetY())):
            if is_sf_out_of_cell(cellx, celly, xsfi,ysfi, args.brick):
                #bad_sfi.append(i)
                #print(f"out of cell {isfi} ({cellx},{celly}), scifi ({xsfi},{ysfi})")
                continue
            elif any_nearby(xsfi,ysfi,tgr_emu, args.thres):
                good_sfi.append(isfi)
                tgr_good_sfi.AddPoint(xsfi,ysfi)
            else:
                #print("not nearby: scifi:",xsfi,ysfi)
                bad_sfi.append(isfi)
        for iemu,(xemu,yemu) in enumerate(zip(tgr_emu.GetX(), tgr_emu.GetY())):
            if any_nearby(xemu,yemu,tgr_sfi, args.thres):
                good_emu.append(iemu)
                tgr_good_emu.AddPoint(xemu,yemu)
            else:
                bad_emu.append(iemu)
        h_cell.SetTitle(f"Cell {cellx}, {celly}, close matches")
        h_cell.Draw("colz")
        tgr_good_emu.Draw("same P")
        tgr_good_sfi.Draw("same P")
        myleg.Clear()
        myleg.SetNColumns(2)
        myleg.AddEntry(tgr_good_emu,"Emu exit showers","P")
        myleg.AddEntry(tgr_good_sfi,"Scifi reco track","P")
        myleg.Draw("same P")
        if not args.goff:
            c.Print(f"{out_pdf}")
        h_dropped_emu.Fill(len(bad_emu))
        h_dropped_sfi.Fill(len(bad_sfi))
        print("Bad emu, bad sfi:", len(bad_emu), len(bad_sfi))
        h_close_matches.Fill( (len(good_emu)+len(good_sfi))/2.0)

        if tgr_good_sfi.GetN() == 0 or tgr_good_emu.GetN() == 0:
            print(f"WARNING: cell {cellx}, {celly} empty after cleaning! (good sfi {tgr_good_sfi.GetN()}, good emu {tgr_good_emu.GetN()})")
            continue
        
        ## move sf to closest emu and calculate cell diff
        l_cell_diff = []
        movex_avg = 0.0
        movey_avg = 0.0
        for isfi, (xsfi,ysfi) in enumerate(zip(tgr_good_sfi.GetX(), tgr_good_sfi.GetY())):
            # we move a particular scifi point i to the closest emu
            dist_sum = 0.0
            closest, movex, movey, min_dist = move_to_closest(xsfi,ysfi,tgr_good_emu)
            h_global_x.Fill(movex)
            h_global_y.Fill(movey)
            print(f"point {isfi}, ({xsfi},{ysfi}), movex {movex}, movey {movey}")
            tgr_sf_move = ROOT.TGraph()
            tgr_sf_move.SetMarkerColor(tgr_good_sfi.GetMarkerColor())
            tgr_sf_move.SetMarkerStyle(tgr_good_sfi.GetMarkerStyle())
            tgr_sf_move.SetMarkerSize(tgr_good_sfi.GetMarkerSize())
            # now calculate a dist moving all other scifi to this position and sum it to total
            movex_avg += movex
            movey_avg += movey
            #print(f"Doing i {i}")
            for jsfi, (xjsfi, yjsfi) in enumerate(zip(tgr_good_sfi.GetX(), tgr_good_sfi.GetY())):
                cloj, movxj, movyj, min_dist_j = move_to_closest(xjsfi+movex,yjsfi+movey,tgr_good_emu)
                dist_sum += min_dist_j
                #print(f"\tj {j} mindist {min_dist_j}, dist_sum update {dist_sum}")
                tgr_sf_move.AddPoint(xjsfi+movex,yjsfi+movey)
                
            
            #print(f"Final dist_sum for i: {dist_sum}")
            l_cell_diff.append((isfi,dist_sum))
            if cellx == 7 and celly == 12:
                h_cell.SetTitle(f"Cell {cellx}, {celly}, move scifi {isfi}, sum of dist/N {dist_sum/tgr_sf_move.GetN():.3f}, movex {movex:.3f}, movey {movey:.3f}")
                h_cell.Draw("colz")
                tgr_good_emu.Draw("same P")
                tgr_sf_move.Draw("same P")
                if not args.goff:
                    c.Print(f"{out_pdf}")
        movex_avg /= tgr_good_sfi.GetN()
        movey_avg /= tgr_good_sfi.GetN()
        h_global_x_avg.Fill(movex_avg)
        h_global_y_avg.Fill(movey_avg)



#ROOT.gStyle.SetOptStat(1)
if not args.goff:
    h_global_x.Draw()
    c.Print(f"{out_pdf}")
    h_global_y.Draw()
    c.Print(f"{out_pdf}")
    h_global_x_avg.Draw()
    c.Print(f"{out_pdf}")
    h_global_y_avg.Draw()
    c.Print(f"{out_pdf}")
print(f"Mean global displacement")
print(f"x: {h_global_x.GetMean(1):.5f} ± {h_global_x.GetMeanError(1):.5f}")
print(f"y: {h_global_y.GetMean(1):.5f} ± {h_global_y.GetMeanError(1):.5f}")
print(f"Mean global displacement (per-cell average)")
print(f"x: {h_global_x_avg.GetMean(1):.5f} ± {h_global_x_avg.GetMeanError(1):.5f}")
print(f"y: {h_global_y_avg.GetMean(1):.5f} ± {h_global_y_avg.GetMeanError(1):.5f}")

def double_gaus_fit(h):
    expr = ("[0]*exp(-0.5*((x-[1])/[2])**2)"   # narrow Gaussian
            "+[3]*exp(-0.5*((x-0.0)/[4])**2)") # wide Gaussian, center locked at 0.0
    gaus = ROOT.TF1("gaus", expr, -0.1, 0.1)
    gaus.SetParameters( h.GetMaximum()/4.0, 0.0, 0.01,
                        h.GetMaximum()/2.0, 0.2)
    fit = h.Fit("gaus","","",-10.0,10.0)
    gaus_sig = ROOT.TF1("gaus_sig","gaus",-0.1,0.1)
    gaus_sig.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1), gaus.GetParameter(2))
    gaus_bkg = ROOT.TF1("gaus_bkg","gaus",-0.1,0.1)
    gaus_bkg.SetParameters(gaus.GetParameter(3), 0.0, gaus.GetParameter(4))
    return gaus, gaus_sig, gaus_bkg

gaus_x, gaus_sig_x, gaus_bkg_x = double_gaus_fit(h_global_x)
gaus_y, gaus_sig_y, gaus_bkg_y = double_gaus_fit(h_global_y)

gaus_x.SetLineColor(ROOT.kBlue)
gaus_y.SetLineColor(ROOT.kBlue)
gaus_sig_x.SetLineColor(ROOT.kGreen)
gaus_sig_y.SetLineColor(ROOT.kGreen)


if not args.goff:
    myleg = ROOT.TLegend(0.25,0.89,0.8,0.94)
    myleg.SetNColumns(3)
    h_global_x.Draw()
    gaus_x.Draw("same")
    gaus_sig_x.Draw("same")
    gaus_bkg_x.Draw("same")
    myleg.AddEntry(gaus_sig_x,f"{gaus_sig_x.GetParameter(0):.2f} exp((x-{gaus_sig_x.GetParameter(1):.3f})/{gaus_sig_x.GetParameter(2):.3f})")
    myleg.AddEntry(gaus_bkg_x,f"{gaus_bkg_x.GetParameter(0):.2f} exp((x-{gaus_bkg_x.GetParameter(1):.3f})/{gaus_bkg_x.GetParameter(2):.3f})")
    myleg.AddEntry(gaus_x,"Sum")
    myleg.Draw("same P")
    c.Print(f"{out_pdf}")


    h_global_y.Draw()
    gaus_y.Draw("same")
    gaus_sig_y.Draw("same")
    gaus_bkg_y.Draw("same")
    myleg.Clear()
    myleg.AddEntry(gaus_sig_y,f"{gaus_sig_y.GetParameter(0):.2f} exp((x-{gaus_sig_y.GetParameter(1):.3f})/{gaus_sig_y.GetParameter(2):.3f})")
    myleg.AddEntry(gaus_bkg_y,f"{gaus_bkg_y.GetParameter(0):.2f} exp((x-{gaus_bkg_y.GetParameter(1):.3f})/{gaus_bkg_y.GetParameter(2):.3f})")
    myleg.AddEntry(gaus_y,"Sum")
    myleg.Draw("same P")
    c.Print(f"{out_pdf}")


    h_global_x_avg.Draw()
    gaus_x.Draw("same")
    gaus_sig_x.Draw("same")
    gaus_bkg_x.Draw("same")
    c.Print(f"{out_pdf}")


    h_global_y_avg.Draw()
    gaus_y.Draw("same")
    gaus_sig_y.Draw("same")
    gaus_bkg_y.Draw("same")
    c.Print(f"{out_pdf}")

# numbers
ROOT.gStyle.SetOptStat(1)
if not args.goff:
    h_initial_emu.Draw()        
    c.Print(f"{out_pdf}")
    h_dropped_emu.Draw()        
    c.Print(f"{out_pdf}")
    h_dropped_sfi.Draw()        
    c.Print(f"{out_pdf}")
    h_close_matches.Draw()        
    c.Print(f"{out_pdf}")


    c.Print(f"{out_pdf}]")
        
