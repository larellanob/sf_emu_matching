import ROOT
ROOT.gROOT.SetBatch(True)
import fedrarootlogon
import rootUtils as ut
from array import array

from macro.cell_geometry import *

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--workdir", default=".")
parser.add_argument("--cellx", dest="cellx", type=int, default=0)
parser.add_argument("--celly", dest="celly", type=int, default=0)
parser.add_argument("--voxthres", dest="vox_thres", type=int, default=4) # voxel
parser.add_argument("--clusthres", dest="clus_thres", type=int, default=0)
parser.add_argument("--binwidth",dest="binwidth",type=int, default=100)
parser.add_argument("--sample", help="data or MC sample tag: V2, FLUKA, ...")
parser.add_argument("--run"  , type=int)
parser.add_argument("--brick", type=int)
parser.add_argument("--debug", action="store_true")
args = parser.parse_args()
ROOT.gStyle.SetOptStat(0)

tfilename = f"{args.workdir}/emu/root/brick{args.brick}_binning{args.binwidth}um_plates_{args.sample}_bintree.root"
tfile = ROOT.TFile(f"{tfilename}")
tree  = tfile.Get("bincontents")

# initialise a th2 for the cell x,y in scifi coordinates
def ini_th2_cell_sf(cellx, celly, title=""):
    if title == "":
        title=f"Cell {cellx}, {celly} Brick {args.brick}"
    # coordinates in emu
    xi, xf, yi, yf = get_cell_borders(cellx, celly, args.brick)
    # coordinates in sf
    xi, yi = emu_to_sf(xi, yi, args.brick)
    xf, yf = emu_to_sf(xf, yf, args.brick)
    if args.debug:
        print(f"emu to sf: (xi, yi): ({xi}, {yi}). (xf, yf): ({xf}, {yf})")
    nbins = int(10000/args.binwidth)
    h2 = ROOT.TH2F(f"h2_{cellx}_{celly}",
                   f"{title}",
                   nbins, xi, xf,
                   nbins, yi, yf)
    #h2.SetDirectory(0)
    return h2
    
def p_range_sum(pi, pf):
    # returns a string to make a selection in TTree.Draw() based on
    # the sum of the bin contents of plates i through plate f
    pc = pi
    sel = f"bp[{pi}]"
    while pc < pf:
        pc += 1
        sel = sel+f"+ bp[{pc}]"
    return sel


def cluster_xy(h2):
    if args.debug:
        print(f"cluster_xy(): start. h2.GetEntries(): {h2.GetEntries()}")
    ax = h2.GetXaxis()
    ay = h2.GetYaxis()
    d_tup = {} 
    d_tup_in_clus = {}
    for bx in range(1,ax.GetNbins()+1):
        for by in range(1,ay.GetNbins()+1):
            bxy = h2.GetBin(bx,by)
            bc = h2.GetBinContent(bxy)
            if bc > 0 and args.debug:
                print(f"bx by bc: {bx} {by} {bc}")
            d_tup[(bx,by)] = bc
            d_tup_in_clus[(bx,by)] = False
    l_tup = list(d_tup) # list of keys
    if args.debug:
        print(f"cluster_xy(): len(l_tup): {len(l_tup)}")

    # check for adjacency between tuples ti, tj
    def is_adjacent(ti,tj):
        if abs(tj[0]-ti[0]) < 2 and abs(tj[1]-ti[1]) < 2:
            return True
        else:
            return False

    # grow the cluster recursively
    def grow_cluster(ti, cluster):
        for tj in l_tup:
            if d_tup_in_clus[tj]:
                continue
            if d_tup[tj] == 0:
                # nothing to cluster
                continue
            if is_adjacent(ti, tj):
                d_tup_in_clus[tj] = True
                cluster.append(tj)
                # recursion
                grow_cluster(tj,cluster)

    l_all_clusters = []
    for ti in l_tup:
        if d_tup_in_clus[ti]:
            # this tuple already in a cluster
            continue
        if d_tup[ti] == 0:
            # nothing to cluster
            continue
        # if not part of existing cluster, start a new one
        cluster = [ti]
        d_tup_in_clus[ti] = True
        # grow the cluster recursively
        grow_cluster(ti, cluster)
        l_all_clusters.append(cluster)

    # transform clusters from index to real coordinates, and put
    # bincontent as 3rd index in a tuple
    l_xy_clusters = []
    for c in l_all_clusters:
        xy_cluster = []
        for ci in c:
            tup = (ax.GetBinCenter(ci[0]), ay.GetBinCenter(ci[1]), d_tup[ci])
            xy_cluster.append(tup)
        l_xy_clusters.append(xy_cluster)
    if args.debug:
        print(f"cluster_xy(): len(l_xy_clusters): {len(l_xy_clusters)}")
    return l_xy_clusters

def tgr_from_clusters(xy_clusters, cluster_thres=args.clus_thres):
    # now l_all_clusters contains all clusters
    tgr = ROOT.TGraph2D()
    for c in xy_clusters:
        cx, cy = 0.0, 0.0
        sum_bc = 0.0
        for ci in c:
            cx += ci[0]*ci[2]
            cy += ci[1]*ci[2]
            sum_bc += ci[2]
        if sum_bc > cluster_thres or cluster_thres == 0:
            tgr.AddPoint(cx/sum_bc, cy/sum_bc, sum_bc)
    return tgr
    
def draw_cell_slices(cellx,
                     celly,
                     c,
                     si = 0,
                     sf = 56,
                     thres = 57*args.vox_thres,
                     cluster=False):
    h_cell = ini_th2_cell_sf(cellx, celly)
    h_cell.Reset()
    h_cell.Draw("colz")
    name = h_cell.GetName()
    if cluster:
        h_cell.SetName(f"{name}_exit")
        name = h_cell.GetName()
    if si == 0 and sf == 56:
        ndraw = tree.Draw(f"sfy:sfx>>+{name}",f"(cellx == {cellx} && celly == {celly})*bc", "colz same")
        if args.debug:
            print(f"draw_cell_slices() (si==0 and sf==56): ndraw: {ndraw}")
    else:
        if args.clus_thres != 0:
            h_cell.SetTitle(f"Cell {cellx}, {celly} Plates {si+1} to {sf+1}: bins above {thres}, clusters above {args.clus_thres};x (cm);y (cm)")
        else:
            h_cell.SetTitle(f"Cell {cellx}, {celly} Plates {si+1} to {sf+1}: bins above {thres};x (cm);y (cm)")
        ndraw = tree.Draw(f"sfy:sfx>>+{name}",
                          f"(cellx == {cellx} && celly == {celly} && ({p_range_sum(si,sf)}>{thres}) )*({p_range_sum(si,sf)})",
                          "colz same")
        if args.debug:
            print(f"draw_cell_slices(): ndraw: {ndraw}")
    tgr = ROOT.TGraph2D()
    if cluster:
        xy_clusters = cluster_xy(h_cell)
        tgr = tgr_from_clusters(xy_clusters)
        tgr.SetName(f"tgr2d_{cellx}_{celly}")
        tgr.SetMarkerStyle(89)
        tgr.SetMarkerSize(4)
        tgr.SetMarkerColor(ROOT.kRed)
        #tgr.Draw("same P")
    #c.Print(f"{out_pdf}")
    #ROOT.gDirectory.Delete(name)
    #del h_cell
    if args.debug:
        print(f"draw_cell_slices(): h_cell.GetEntries(): {h_cell.GetEntries()}")
    return c,tgr, h_cell
    

c = ROOT.TCanvas()
c.SetCanvasSize(1000,1000)
out_pdf = f"{args.workdir}/emu/img/bin_tree_{args.sample}_b{args.run}{args.brick}_thvox{args.vox_thres}_thcl{args.clus_thres}.pdf"
if args.cellx != 0 or args.celly!= 0:
    out_pdf = f"{args.workdir}/emu/img/bin_tree_{args.sample}_b{args.run}{args.brick}_cell{args.cellx}.{args.celly}_thvox{args.vox_thres}_thcl{args.clus_thres}.pdf"

c.Print(f"{out_pdf}[")

'''
# tests
tree.Draw("bc")
c.Print(f"{out_pdf}")
tree.Draw("sfy:sfx","bc > 400")
c.Print(f"{out_pdf}")
tree.Draw("by:bx","bc > 400")
c.Print(f"{out_pdf}")
'''

# integer divisors of 56 (ignoring first plate)
# 1,56
# 2,28
# 4,14
# 7,8
# 8,7
# 14,4
# 28,2
# 56,1
slice_width=7
nslices=8
slice_thres = slice_width*args.vox_thres

'''
# if you want all slices:
for i in range(0,nslices):
    do_cluster = False
    if i == nslices-1:
        print(f"will cluster slice {i}, si {i*slice_width+1}, sf {i*slice_width+slice_width}")
        do_cluster=True
    else:
        continue
    c, tgr = draw_cell_slices(5,15,c,out_pdf, si=i*slice_width+1, sf = i*slice_width+slice_width, thres=thres, cluster=do_cluster)
'''

tfile_tgr = ROOT.TFile(f"{tfilename}".replace(".root",f"_exitshw_voxthres{args.vox_thres}_clus{args.clus_thres}.root"), "recreate")

# last slice for all cells
for x in range(1,19):
    if args.cellx != 0 and x != args.cellx:
        print(f"Skipping x = {x}")
        continue
    for y in range(1,19):
        if args.celly != 0 and y != args.celly:
            print(f"\tSkipping y = {y}")
            continue

        print(f"Doing cell {x}, {y}")
        c, tgr2d, hcell_exit = draw_cell_slices(x,y,c, si=56-slice_width, sf = 56, thres=slice_thres, cluster=True)
        tgx = tgr2d.GetX()
        tgy = tgr2d.GetY()
        tgn = tgr2d.GetN()
        if tgn == 0:
            print(f"CELL {x}, {y} has no valid points!")
            continue
        tgr = ROOT.TGraph(tgn,tgx,tgy)
        tgr.SetMarkerColor(tgr2d.GetMarkerColor())
        tgr.SetMarkerStyle(tgr2d.GetMarkerStyle())
        tgr.SetMarkerSize(tgr2d.GetMarkerSize())
        tgr.Draw("same P")
        tgr.SetName(tgr2d.GetName().replace("tgr2d","tgr"))
        tgr.Write()
        tgr2d.Write()
        hcell_exit.Write()

        c.Print(f"{out_pdf}")
        c, tgr2, hcell_full = draw_cell_slices(x,y,c) # draw whole cell
        tgr.Draw("same P")
        hcell_full.Write()
        c.Print(f"{out_pdf}")

c.Print(f"{out_pdf}]")
