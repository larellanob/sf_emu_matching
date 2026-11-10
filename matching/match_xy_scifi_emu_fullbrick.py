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
import numpy as np
from scipy.spatial import cKDTree

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cellx", dest="cellx",default=0,type=int, help="If default, runs all X in a loop")
parser.add_argument("--celly", dest="celly",default=0,type=int, help="If default, runs all Y in a loop")
parser.add_argument("--xmin", dest="xmin",default=0,type=float)
parser.add_argument("--xmax", dest="xmax",default=0,type=float)
parser.add_argument("--ymin", dest="ymin",default=0,type=float)
parser.add_argument("--ymax", dest="ymax",default=0,type=float)
parser.add_argument("--step", dest="step",default=0,type=float)
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

from math import sqrt
def dist(x1, y1, x2, y2):
    return sqrt( (x2-x1)**2 + (y2-y1)**2 )

def make_point_tree(tgr):
    # make tree using scifi (the one who is static)
    points = np.column_stack((tgr.GetX(), tgr.GetY()))
    return cKDTree(points)

def how_many_nearby(tgr1, kdtree, thres=0.06):
    points1 = np.column_stack((tgr1.GetX(), tgr1.GetY()))

    # numbe of points in kdtree within threshold of each point in tgr1
    nearby = kdtree.query_ball_point(points1,thres)

    # count points in tgr1 that have at least one point nearby in kdtree
    nmatches = sum(bool(matches) for matches in nearby)

    pairs = [ (i,j) for i, matches in enumerate(nearby) for j in matches ]

    return nmatches, pairs

ROOT.gStyle.SetPadLeftMargin(0.150)

c = ROOT.TCanvas()
c.SetCanvasSize(1100,1000)
out_pdf = f"img/match_xy_scifi_emu_brick{args.brick}_th{args.thres}_step{args.step}.pdf"
if args.data:
    out_pdf = out_pdf.replace(".pdf","_data.pdf")
if not args.goff:
    c.Print(f"{out_pdf}[")

# combine all tgraphs in a full brick tgraph
tgr_brick = ROOT.TGraph()
for cellx in range(1,19):
    for celly in range(1,19):
        tgr_emu = tfile_emu.Get(f"tgr_{cellx}_{celly}")
        if tgr_emu == None:
            print(f"Could not find cell {cellx}, {celly}!")
            continue
        for point in range(tgr_emu.GetN()):
            pointx = tgr_emu.GetPointX(point)
            pointy = tgr_emu.GetPointY(point)
            tgr_brick.AddPoint(pointx, pointy)

# arrays of original x,y positions of tgr_brick
x_orig = np.array(tgr_brick.GetX(), copy=True)
y_orig = np.array(tgr_brick.GetY(), copy=True)

def translate_tgr(tgr, x_orig, y_orig, tx, ty):
    for p in range(tgr.GetN()):
        tgr.SetPoint(p, x_orig[p] + tx, y_orig[p] + ty)

print(f"full brick points: {tgr_brick.GetN()}")

# kd-tree on scifi
kdtree = make_point_tree(tgr_sfi)

# loop over translations, fill th2f with n coincidences
nbinsx = int(round((args.xmax-args.xmin)/args.step,5))+1
nbinsy = int(round((args.ymax-args.ymin)/args.step,5))+1
halfstep = 0.5*args.step


th2_coincidences = ROOT.TH2F("th2_coincidences",
                             f"#splitline{{Number of coincidences between Scifi-emulsion}}{{Step size: {args.step} (cm), threshold: {args.thres} (cm)}};Emulsion displacement in #it{{x}} (cm);Emulsion displacement in #it{{y}} (cm)",
                             nbinsx, args.xmin-halfstep, args.xmax+halfstep,
                             nbinsy, args.ymin-halfstep, args.ymax+halfstep)

th2_coincidences.GetXaxis().SetNdivisions(505)

# we save the coincidences
nmatches_max = 0
nmatches_max_x = -999.0
nmatches_max_y = -999.0
nmatches_max_pairs = []

x = args.xmin
y = args.ymin
print(f"point 28 (initial): {tgr_brick.GetPointX(28):.3f}, {tgr_brick.GetPointY(28):.3f}")
while x <= args.xmax:
    while y <= args.ymax:
        translate_tgr(tgr_brick, x_orig, y_orig, x, y)
        # count coincidences here
        n, pairs = how_many_nearby(tgr_brick, kdtree, args.thres)
        if n > nmatches_max:
            nmatches_max = n
            nmatches_max_x = x
            nmatches_max_y = y
            nmatches_max_pairs = pairs
        th2_coincidences.Fill(x,y,n)
        # finish counting coincidences
        print(f"translate loop {x:.4f} {y:.4f}. point 28: {tgr_brick.GetPointX(28):.4f}, {tgr_brick.GetPointY(28):.4f}, matches: {n}")
        y = round(y + args.step,5)
    x = round(x + args.step,5)
    y = args.ymin

th2_coincidences.Draw("colz")
c.Print(f"{out_pdf}")
c.Print(f"{out_pdf}]")

print(f"Max number of matches {nmatches_max} at displacement ({nmatches_max_x:.4f}, {nmatches_max_y:.4f})")
print(f"matching pairs, {nmatches_max_pairs}")
