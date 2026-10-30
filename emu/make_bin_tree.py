import ROOT
ROOT.gROOT.SetBatch(True)
import fedrarootlogon
import rootUtils as ut
from array import array

from utils.cell_geometry import *

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--brick", dest="brick", type=int, default=21)
parser.add_argument("--binwidth",dest="binwidth",type=int, default=100)
parser.add_argument("--input",dest="input", required=True)
args = parser.parse_args()


tfilename = args.input
'''
tfilename = f"root/brick{args.brick}_binning{args.binwidth}um_plates.root"
if args.sample == "V2":
    tfilename = f"root/brick{args.brick}_V2_binning{args.binwidth}um_plates.root"
elif args.sample == "FLUKA":
    tfilename = f"root/brick{args.brick}_FLUKA_binning{args.binwidth}um.root"
if args.data:
    tfilename = f"root/brick{args.brick}_binning{args.binwidth}um_plates_data.root"
'''
tfile = ROOT.TFile(f"{tfilename}")

toutfile = ROOT.TFile(tfilename.replace(".root","_bintree.root"),"recreate")
out_tree = ROOT.TTree("bincontents",f"Brick {args.brick} bin contents (binning {args.binwidth} um)")

f_cellx        = array('i',[0])
f_celly        = array('i',[0])
f_bincontent   = array('f',[0.0])
f_bx           = array('f',[0.0])
f_by           = array('f',[0.0])
f_scifix       = array('f',[0.0])
f_scifiy       = array('f',[0.0])
f_platecontent = array('f',[0.0]*57)

out_tree.Branch("cellx",f_cellx,"cellx/I")
out_tree.Branch("celly",f_celly,"celly/I")
out_tree.Branch("bc",f_bincontent,"bc/F")
out_tree.Branch("sfx",f_scifix,"sfx/F") # x-position in scifi coordinates
out_tree.Branch("sfy",f_scifiy,"sfy/F") # y-position in scifi coordinates
out_tree.Branch("bp",f_platecontent,"bp[57]/F")
out_tree.Branch("bx",f_bx,"bx/F")
out_tree.Branch("by",f_by,"by/F")


out_tree_shw = ROOT.TTree("showers",f"Brick {args.brick} showers")
f_sh_startx = array('i',[0])
f_sh_starty = array('i',[0])
f_sh_startz = array('i',[0])
f_sh_endx = array('i',[0])
f_sh_endy = array('i',[0])
f_sh_endz = array('i',[0])
f_sh_content = array('i'[0])

class fedra_voxel:
    # voxel centers
    x = 0.0
    y = 0.0
    z = 0.0
    # widths
    wx = 0.0
    wy = 0.0
    wz = 0.0
    # content
    vc = 0.0

    def __init__(self,x,y,z,wx,wy,wz,vc):
        self.x = x
        self.y = y
        self.z = z
        self.wx = wx
        self.wy = wy
        self.wz = wz
        self.vc = vc

    def get_vc(self):
        return self.vc

    def get_x(self):
        return self.x
    def get_y(self):
        return self.y
    def get_z(self):
        return self.z
        

class fedra_cluster:
    voxels = []
    startx = 0.0
    starty = 0.0
    startz = 0.0
    endx = 0.0
    endy = 0.0
    endz = 0.0
    # simple cone angle alpha
    aperture = 0.0
    sum_n = 0.0

    def __init__(self, voxels):
        self.voxels = voxels
        self.set_total_s0()
        self.set_start_end()

    def set_start_end(self):
        z_sorted = self.voxels.sort(key=lambda v: v.get_z())
        self.startx = z_sorted[0].get_x()
        self.starty = z_sorted[0].get_y()
        self.startz = z_sorted[0].get_z()
        self.endz = z_sorted[-1].get_z()

    def set_total_s0(self):
        for v in voxels:
            self.sum_n += v.get_vc()

    def add_voxel(self,v):
        self.voxels.append(v)


# clusterization
clus_thres = 2 # min number of sf[0] to define voxel
clus_sepx  = 1 # maximum separation in x
clus_sepy  = 1 # maximum separation in y
clus_sepz  = 1 # maximum separation in z

print("Will loop through histograms")
#for cellx in range(11,12):
for cellx in range(1,19):
    f_cellx[0] = cellx
#    for celly in range(3,4):
    for celly in range(1,19):
        f_celly[0] = celly
        h2_cell = tfile.Get(f"cells/c{cellx}_{celly}/h2_{cellx}_{celly}")
        if h2_cell == None:
            print(f"Missing cell {cellx}, {celly}!")
            continue
        else:
            print(f"Working over cell {cellx} {celly}")
        h2_plates = {}
        for p in range(1,58):
            h2_plates[p] = tfile.Get(f"cells/c{cellx}_{celly}/h2_{cellx}_{celly}_p{p}")

        # ++++n muons?+++++ (don't have this info) (first plate is good indicator, actually!)
        # bincontent pi (array length 57)
        # scifix, scifiy
        voxels = []
        ax = h2_cell.GetXaxis()
        ay = h2_cell.GetYaxis()
        for bx in range(1,ax.GetNbins()+1):
            f_bx[0] = ax.GetBinCenter(bx)

            center_x = ax.GetBinCenter(bx)
            for by in range(1,ay.GetNbins()+1):
                f_by[0] = ay.GetBinCenter(by)
                center_y = ay.GetBinCenter(by)
                f_scifix[0], f_scifiy[0] = emu_to_sf(center_x, center_y, args.brick)
                bxy = h2_cell.GetBin(bx,by)
                f_bincontent[0] = h2_cell.GetBinContent(bxy)
                for p in range(1,58):
                    p_content = h2_plates[p].GetBinContent(bxy)
                    f_platecontent[p-1] = p_content

                    '''
                    if p > 2 and p_content > clus_thres:
                        myvox = fedra_voxel(ax.GetBinCenter(bx),
                                            ay.GetBinCenter(by),
                                            p,
                                            args.binwidth,
                                            args.binwidth,
                                            1.0,
                                            p_content
                                            )
                        voxels.append(myvox)
                    '''
                out_tree.Fill()
        
print("Finished looping over histograms")
out_tree.Write()
