import ROOT
ROOT.gROOT.SetBatch(True)
import fedrarootlogon
import rootUtils as ut
import array


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--brick", dest="brick", type=int, default=21)
parser.add_argument("--run", dest="run", type=int, default=1)
parser.add_argument("--binwidth",dest="binwidth",type=int, default=100)
parser.add_argument("--cellx", dest="cellx", type=int, default=9)
parser.add_argument("--celly", dest="celly", type=int, default=12)
parser.add_argument("--data", action='store_true')
parser.add_argument("--fluka", action='store_true')
args = parser.parse_args()
ROOT.gStyle.SetOptStat(0)

if args.cellx > 18 or args.celly > 18 or args.cellx < 1 or args.celly < 1:
    print(f"ERROR, invalid cell {args.cellx} {args.celly}")
    exit()
cell0 = f"{args.cellx}0_{args.celly}0"
celldot = f"{args.cellx}.{args.celly}"

brick = args.brick
base_dir = f"/eos/experiment/sndlhc/users/dancc/FEDRA/muon_Euniform_RUN{args.run}_V2/cell_reco"
if args.fluka:
    base_dir = f"/eos/experiment/sndlhc/users/dancc/FEDRA/muon_Euniform_RUN{args.run}_FLUKA25/b0000{brick}/cell_reco"
if args.data:
    base_dir = f"/eos/experiment/sndlhc/emulsionData/emureco_Napoli/RUN{args.run}/b000{args.run}{brick}/cells/"
print("BASE DIR:",base_dir)
work_dir = f"/afs/cern.ch/user/l/luarella/work-luarella/public/sf_emu_matching/emu"

def get_cell_borders(x,y,brick, data):
    width=10000
    extra=0
    off_x, off_y = 0, -500
    if brick%10 == 3 or brick%10 == 1:
        off_x = 195000
    if brick%10 == 3 or brick%10 == 4:
        off_y = 195000
    if data:
        off_x = 0
        off_y = 0
        
    xi = off_x + x*width - width/2 - extra
    xf = off_x + x*width + width/2 + extra
    yi = off_y + y*width - width/2 - extra
    yf = off_y + y*width + width/2 + extra
    return xi, xf, yi, yf

xi, null1, yi, null2 = get_cell_borders(1,1,brick, args.data)
print(xi,null1,yi,null2)
null1,  xf, null2, yf = get_cell_borders(18,18,brick, args.data)
print(null1,xf,null2,yf)
print(xi,yi,xf,yf)

xi_cell, xf_cell, yi_cell, yf_cell = get_cell_borders(args.cellx, args.celly, brick, args.data)

um = 10000*18
binwidth=args.binwidth #um
nbins=int(um/binwidth)
nbins_cell = int((um/18)/binwidth)

print("Full histogram borders x:", xi, xf)
print("Full histogram borders y:", yi, yf)
h2_full = ROOT.TH2F("h2_full",
                    f"Brick {brick} ({binwidth}x{binwidth} #mum bins);x (#mum);y (#mum)",
                    nbins, xi, xf,
                    nbins, yi, yf
                    )
h2_cell = ROOT.TH2F(f"h2_{args.cellx}_{args.celly}",
                    f"Cell {args.cellx}, {args.celly} Brick {brick} ({binwidth}x{binwidth} #mum bins);x (#mum);y (#mum)",
                    nbins_cell, xi_cell, xf_cell,
                    nbins_cell, yi_cell, yf_cell
                    )
# h2_cell by depth
h2 = {}
for p in range(1,58):
    h2[f"h2_{args.cellx}_{args.celly}_p{p}"] = ROOT.TH2F(
        f"h2_{args.cellx}_{args.celly}_p{p}",
        f"Plate {p} Cell {args.cellx}, {args.celly} Brick {brick} ({binwidth}x{binwidth} #mum bins);x (#mum);y (#mum)",
        nbins_cell, xi_cell, xf_cell,
        nbins_cell, yi_cell, yf_cell
    )

out_filename = f"{work_dir}/root/cells_V2/brick{brick}_cell_{args.cellx}.{args.celly}_binning{binwidth}um.root"
if args.fluka:
    out_filename = f"{work_dir}/root/cells_fluka/brick{brick}_cell_{args.cellx}.{args.celly}_binning{binwidth}um.root"
if args.data:
    out_filename = f"{work_dir}/root/data/brick{brick}_cell_{args.cellx}.{args.celly}_binning{binwidth}um.root"
print("FILE OUT:",out_filename)
out_file = ROOT.TFile(out_filename,"recreate")
        

if args.data:
    file_in = ROOT.TFile(f"{base_dir}/cell_{cell0}/b000{args.run}{brick}/b000{args.run}{brick}.0.0.0.trk.root")
else:
    if args.fluka:
        file_in = ROOT.TFile(f"{base_dir}/cell_{cell0}/b0000{brick}/b0000{brick}.0.{celldot}.trk.root")
    else:
        file_in = ROOT.TFile(f"{base_dir}/cell_{cell0}/b0000{brick}/b0000{brick}.0.{celldot}.trk.root")
print("FILE IN:",file_in)

tree = file_in.Get("tracks")
entries = 0
for ievent, event in enumerate(tree):
    entries += 1
    h2_full.Fill(event.sf[0].X(),event.sf[0].Y())
    h2_cell.Fill(event.sf[0].X(),event.sf[0].Y())
    for p in range(1,58):
        if event.s[0].Plate() == p:
            h2[f"h2_{args.cellx}_{args.celly}_p{p}"].Fill(event.sf[0].X(),event.sf[0].Y())

print(f"Cell {args.cellx},{args.celly} added {entries} entries")

out_file.cd()
h2_full.Write()

ROOT.gDirectory.mkdir(f"cells")
ROOT.gDirectory.cd(f"cells")
ROOT.gDirectory.mkdir(f"c{args.cellx}_{args.celly}")
ROOT.gDirectory.cd(f"c{args.cellx}_{args.celly}")
h2_cell.Write()
for h in h2:
    h2[h].Write()
