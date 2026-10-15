import ROOT
ROOT.gROOT.SetBatch(True)
import fedrarootlogon
import rootUtils as ut
from array import array

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--brick", dest="brick", type=int, default=21)
parser.add_argument("--binwidth",dest="binwidth",type=int, default=100)
parser.add_argument("--cellx", dest="cellx", type=int, default=9)
parser.add_argument("--celly", dest="celly", type=int, default=15)
parser.add_argument("--minclusters", dest="min_clusters", type=int, default=1)
parser.add_argument("--sample", dest="sample", default="FLUKA") # FLUKA, V2, data
parser.add_argument("--threshold", dest="threshold", type=int, default=4) # number of new tracks per bin for clustering
parser.add_argument("--maxgap", dest="maxgap", type=int, default=1) # max gap for clustering in xy and z
args = parser.parse_args()

nbins_cell = int(10000/args.binwidth)
def get_cell_borders(x,y,brick,mc=True):
    width=10000
    extra=0
    off_x, off_y = 0, -500
    if brick%10 == 3 or brick%10 == 1:
        off_x = 195000
    if brick%10 == 3 or brick%10 == 4:
        off_y = 195000
    if not mc:
        off_x = 0
        off_y = 0
    xi = off_x + x*width - width/2 - extra
    xf = off_x + x*width + width/2 + extra
    yi = off_y + y*width - width/2 - extra
    yf = off_y + y*width + width/2 + extra
    return xi, xf, yi, yf
is_mc = True
if args.sample == "data":
    is_mc = False
xi_cell, xf_cell, yi_cell, yf_cell = get_cell_borders(args.cellx, args.celly, args.brick, is_mc)

# work dir
eos_dir = f"/eos/experiment/sndlhc/users/luarella"
work_dir = f"/afs/cern.ch/user/l/luarella/work-luarella/public/sf_emu_matching"

# fedra file
if args.sample == "V2":
    fedra_dir   = f"/eos/experiment/sndlhc/users/dancc/FEDRA/muon_Euniform_RUN1_V2/cell_reco"
    histo_dir = f"{eos_dir}/emu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_spectregen_18cmB21_transpVeto_SWANBINNING_V2"
    output_dir = f"{work_dir}/emu/root/clusters_V2"
elif args.sample == "FLUKA":
    fedra_dir   = f"/eos/experiment/sndlhc/users/dancc/FEDRA/muon_Euniform_RUN1_FLUKA25/b0000{args.brick}/cell_reco_retrk"
    histo_dir = f"{eos_dir}/emu/ALL_lhc_ir1_coll_2023_exp001_fort_spectregen_SWANBINNING"
    output_dir = f"{work_dir}/emu/root/fluka"
elif args.sample == "data":
    fedra_dir = f"/eos/experiment/sndlhc/emulsionData/emureco_Napoli/RUN1/b000121/cells"
    histo_dir = f"{eos_dir}/emu/data"
    output_dir = f"{work_dir}/emu/root/clusters_data"

fedra_fname = f"{fedra_dir}/cell_{args.cellx}0_{args.celly}0/b0000{args.brick}/b0000{args.brick}.0.{args.cellx}.{args.celly}.trk.root"
if args.sample == "data":
    fedra_fname = f"{fedra_dir}/cell_{args.cellx}0_{args.celly}0/b0000{args.brick}/b0000{args.brick}.0.0.0.trk.root"

# histogram file
histo_fname = f"{histo_dir}/brick21_cell_{args.cellx}.{args.celly}_binning{args.binwidth}um.root"
try:
    histo_tfile = ROOT.TFile(f"{histo_fname}")
    print("Histo file opened OK: ", histo_tfile)
except:
    print("WARNING: histo tfile does not exist. Will make the histograms (~2.5 min)")
    # if file doesn't exist, make histograms yourself
    h2 = {}
    for p in range(1,58):
        h2[f"h2_{args.cellx}_{args.celly}_p{p}"] = ROOT.TH2F(
            f"h2_{args.cellx}_{args.celly}_p{p}",
            f"Plate {p} Cell {args.cellx}, {args.celly} Brick {args.brick} ({args.binwidth}x{args.binwidth} #mum bins);x (#mum);y (#mum)",
            nbins_cell, xi_cell, xf_cell,
            nbins_cell, yi_cell, yf_cell
        )

    # fill histograms
    print("Will fill plate histograms")
    histo_tfile = ROOT.TFile(f"{histo_fname}","recreate")
    fedra_tfile = ROOT.TFile(f"{fedra_fname}")
    fedra_ttree = fedra_tfile.Get("tracks")
    for ievent, event in enumerate(fedra_ttree):
        p = event.s[0].Plate()
        h2[f"h2_{args.cellx}_{args.celly}_p{p}"].Fill(event.sf[0].X(),event.sf[0].Y())
    print("finished filling plate histograms")
    histo_tfile.cd()

    cells_dir = ROOT.gDirectory.mkdir(f"cells")
    if not cells_dir:
        print("ERROR: failed to create cells dir")
        exit()
    cells_subdir = cells_dir.mkdir(f"c{args.cellx}_{args.celly}")
    if not cells_subdir:
        print("ERROR: failed to create cells SUBdir")
        exit()
    cells_subdir.cd()
    for h in h2:
        h2[h].Write()
    print("Finished filling plate histograms")

print("Opening FEDRA file")    
fedra_tfile = ROOT.TFile(f"{fedra_fname}")
fedra_ttree = fedra_tfile.Get("tracks")
print("Fedra file: ", fedra_tfile)

# clusterize
def cluster_xy(h2, plate, threshold, maxgap=0):
    ax = h2.GetXaxis()
    ay = h2.GetYaxis()
    d_tup = {} 
    d_tup_in_clus = {}
    for bx in range(1,ax.GetNbins()+1):
        for by in range(1,ay.GetNbins()+1):
            bxy = h2.GetBin(bx,by)
            bc = h2.GetBinContent(bxy)
            d_tup[(bx,by,plate)] = bc
            d_tup_in_clus[(bx,by,plate)] = False
    l_tup = list(d_tup) # list of keys

    # check for adjacency between tuples ti, tj
    def is_adjacent(ti,tj,maxgap):
        if abs(tj[0]-ti[0]) < 2+maxgap and abs(tj[1]-ti[1]) < 2+maxgap:
            return True
        else:
            return False

    # grow the cluster recursively
    def grow_cluster(ti, cluster):
        for tj in l_tup:
            if d_tup_in_clus[tj]:
                continue
            if d_tup[tj] < threshold:
                # nothing to cluster
                continue
            if is_adjacent(ti, tj,maxgap):
                d_tup_in_clus[tj] = True
                cluster.append(tj)
                # recursion
                grow_cluster(tj,cluster)

    l_all_clusters = []
    for ti in l_tup:
        if d_tup_in_clus[ti]:
            # this tuple already in a cluster
            continue
        if d_tup[ti] < threshold:
            # nothing to cluster
            continue
        # if not part of existing cluster, start a new one
        cluster = [ti]
        d_tup_in_clus[ti] = True
        # grow the cluster recursively
        grow_cluster(ti, cluster)
        l_all_clusters.append(cluster)

    # transform clusters from index to real coordinates, and put
    # bincontent as 4th index in a tuple
    l_xy_clusters = []
    for c in l_all_clusters:
        xy_cluster = []
        for ci in c:
            #tup = (ax.GetBinCenter(ci[0]), ay.GetBinCenter(ci[1]), plate, d_tup[ci])
            tup = (ci[0], ci[1], plate, d_tup[ci])
            xy_cluster.append(tup)
        l_xy_clusters.append(xy_cluster)
    return l_xy_clusters



def cluster_z(lcxyp1, lcxyp2,maxgap=0):
    '''clusterization in z. input is list of xy-clusters for two plates.

    the list of clusters 1 will "absorb" as many clusters in 2 as it
    can
    
    '''
    def is_adjacent_cluster(ci, cj, maxgap):
        for ti in ci:
            for tj in cj:
                if abs(ti[0]-tj[0]) < 2+maxgap and abs(ti[1]-tj[1]) < 2+maxgap and abs(ti[2]-tj[2]) < 2+maxgap:
                    return True
        return False

    # recursively combine adjacent clusters
    def combine_clusters(lc1,lc2):
        for i,ci in enumerate(lc1):
            for j, cj in enumerate(lc2):
                if is_adjacent_cluster(ci,cj,maxgap):
                    lc1[i] = ci+cj
                    lc2.pop(j)
                    # will keep combining until no more adjacent
                    combine_clusters(lc1,lc2)
        # exhsaused combinations. eat up whatever is left
        while len(lc2) > 0:
            for j, cj in enumerate(lc2):
                lc2.pop(j)
                lc1.append(cj)
        return lc1

    # shallow copies of the input lists, so they are not modified
    l1 = lcxyp1[:] or list(lcxyp1)
    l2 = lcxyp2[:] or list(lcxyp2)
    
    l_cluster_sum = []
    l_cluster_sum = combine_clusters(l1,l2)

    # finished
    return l_cluster_sum
                
        
# first, per-plate clusters in xy
xy_clus = {}
cell_xyz_clusters = []
print("Will make XY clusters")
for p in range(3,57): # FEDRA tracks cannot start in p=57, so excluded
    print(f"clustering plate {p}")
    h2_p = histo_tfile.Get(f"cells/c{args.cellx}_{args.celly}/h2_{args.cellx}_{args.celly}_p{p}")
    thres = args.threshold
    if h2_p.GetEntries() > 20000:
        thres+=1
    xy_clus[p] = cluster_xy(h2_p, p, thres,args.maxgap)
print("Will make Z clusters")
for p in range(3,57):
    print(f"Z-clustering up to plate {p}")
    cell_xyz_clusters = cluster_z(cell_xyz_clusters, xy_clus[p],args.maxgap)

print("Clusterized  XYZ")
# the clustering is done at this point, but i want to be fancy and
# make a class
class cluster_shower():
    is_mc = False
    cell = (0,0)
    brick = 21
    binning = 0
    # list of 3d-bins
    l_clusters = []
    # th2f for coordinate conversion
    th2f = ROOT.TH2F()
    # clusters in emulsion coordinates
    l_clus_emu = []
    # clusters in detector coordinates
    l_clus_det = []

    # total number of first segments (sf0)
    sum_sf0 = 0.0
    
    # initial and end plates
    p_ini = 0
    p_end = 0

    def __init__(self,clusters,x=0,y=0,binning=100):
        # clusters are one 3d bin. list of tuples with 4 elemtents
        # (bin x, bin y, bin z, content)
        # binx and biny depend on binning argument in constructor
        # binz is plate
        # content is the number of new tracks in that bin
        self.l_clusters = clusters
        self.set_total_sf0()
        self.set_start_end_z()

        self.cell = (x,y)
        self.binning = binning
       
    def set_start_end_z(self):
        z_sorted = sorted(self.l_clusters,key=lambda v: v[2])
        self.p_ini = z_sorted[0][2]
        self.p_end = z_sorted[-1][2]

    def set_total_sf0(self):
        for v in self.l_clusters:
            self.sum_sf0 += v[3]
        
    def get_start_plate(self):
        return self.p_ini
    def get_end_plate(self):
        return self.p_end
    def get_max_plate(self):
        d_per_plate = {}
        for p in range(self.p_ini, self.p_end+1):
            d_per_plate[p] = 0
        for v in self.l_clusters:
            d_per_plate[v[2]] += v[3]
        return max(d_per_plate, key=d_per_plate.get)

    def get_content(self):
        return self.sum_sf0
    def get_nclusters(self):
        return len(self.l_clusters)
    def get_cluster(self,i):
        return self.l_clusters[i]
    def get_emu_coord(self,cluster):
        x_index = cluster[0]
        y_index = cluster[1]
        z_index = cluster[2]
        xi, xf, yi, yf = get_cell_borders(self.cell[0], self.cell[1], self.brick)
        ndiv = int(10000/self.binning)
        stepx = (xf-xi)/ndiv
        stepy = (yf-yi)/ndiv
        edgesx = [ xi+i*stepx for i in range(ndiv+1) ]
        edgesy = [ yi+i*stepy for i in range(ndiv+1) ]
        # x index is actually a bin number, from 1 to nbins+1
        # edges is a python list, from 0 to nbins
        centerx = (edgesx[x_index-1]+edgesx[x_index])/2.0
        centery = (edgesy[y_index-1]+edgesy[y_index])/2.0
        centerz = 3040000 + z_index*1015
        return centerx, centery, centerz

    def get_xy_center_at_plate(self,p):
        # returns the x, y center (in emu coordinates) for the weighted clusters in plate p

        # first we select only the clusters in plate p
        l_clus_plate_p = [ c for c in self.l_clusters if c[2] == p ]

        x_wt = 0.0
        y_wt = 0.0
        sum_tracks = 0
        # we loop over them and get a center weighted by 
        for c in l_clus_plate_p:
            sum_tracks += c[3]
            x, y, z = self.get_emu_coord(c)
            x_wt += x*c[3]
            y_wt += y*c[3]
        x_wt /= float(sum_tracks)
        y_wt /= float(sum_tracks)
        return x_wt, y_wt

    def get_boxy_dict(self,clusterid,boxy_dict={}):
        """Boxy selection

        This returns a dictionary to make a quick selection of a
        bounding box around the shower. Circumscribed cuboid, with
        sides parallel to its respective coordinate axis.

        The idea is to use the dict in a ROOT selection over the full
        FEDRA track file, with ~270k tracks. The selection should
        return around 500 tracks. The processing time goes from a few
        minutes to less than a second.
        """
        minx = 999999999
        miny = 999999999
        maxx = -1
        maxy = -1
        for i in range(self.get_nclusters()):
            clus = self.get_cluster(i)
            x,y,z = self.get_emu_coord(clus)
            z = clus[2] # actually, for z use plate, not coordinate
            nbin_half = args.binwidth/2
            if x-nbin_half < minx:
                minx = x-nbin_half
            if y-nbin_half < miny:
                miny = y-nbin_half
            if x+nbin_half > maxx:
                maxx = x+nbin_half
            if y+nbin_half > maxy:
                maxy = y+nbin_half
        # value will be tuple (minx, miny, minz, maxx, maxy, maxz
        boxy_dict[clusterid] = (minx, miny, self.get_start_plate(), maxx, maxy, self.get_end_plate())
        return boxy_dict

    def get_fedra_selection(self, full_selection=""):
        """Fully grained cluster selection
        
        This function is supposed to be used like boxy selection, but
        fits the cluster perfectly. In practice ROOT cannot take such
        a large number of selection operators, so it fails for showers
        with 100 bins/clusters or so. The boxy selection gives enough
        of a time reduction in my opinion.

        """
        # intended to be used with a single cluster
        for i in range(self.get_nclusters()):
            c = self.get_cluster(i)
            x,y,z = self.get_emu_coord(c)
            z = c[2] # actually, for z use plate, not coordinate
            nbin_half = args.binwidth/2
            selx = f"sf[0].X() > {x-nbin_half} && sf[0].X() < {x+nbin_half}"
            sely = f"sf[0].Y() > {y-nbin_half} && sf[0].Y() < {y+nbin_half}"
            selz = f"s[0].Plate() == {z}"
            sel = f"{selx} && {sely} && {selz}"
            if full_selection != "":
                full_selection += f"|| ( {sel} ) "
            else:
                full_selection += f"({sel}) "
        return full_selection

    
    def is_xyz_in_cluster(self,x,y,z):
        xyz_is_in_cluster = False
        # x, y in emulsion coordinates
        # z in plate
        half_bin = self.binning / 2.0
        for i in range(self.get_nclusters()):
            c = self.get_cluster(i)
            if z == c[2]:
                xi,yi,zi = self.get_emu_coord(c)
                if x > xi-half_bin and x < xi+half_bin and y > yi-half_bin and y < yi+half_bin:
                    xyz_is_in_cluster = True
                    break
        return xyz_is_in_cluster
        

#xi, xf, yi, yf = get_cell_borders(args.cellx, args.celly, args.brick)
h3 = ROOT.TH3F("h3",
               f"Cluster showers in cell {args.cellx},{args.celly};z (plate);x (#mum);y (#mum)",
               57, 0, 57,
               nbins_cell, xi_cell, xf_cell,
               nbins_cell, yi_cell, yf_cell
               )

print(f"Printing showers above {args.min_clusters} clusters")
print("Shower\t clusters\t startplate\t endplate\t content\t cluster 0\t\t c0 x,y,z")
d_cluster_sel = {}
n_sum_bins = 0
for i,c in enumerate(cell_xyz_clusters):
    """
    this loop does NOT make the clusters, it only makes a "boxy
    selection"
    """
    myshower = cluster_shower(c,args.cellx, args.celly, args.binwidth)
    if myshower.get_nclusters() >= args.min_clusters:
        c0 = myshower.get_cluster(0)
        x,y,z = myshower.get_emu_coord(c0)
        print(f"{i}\t {myshower.get_nclusters()}\t\t {myshower.get_start_plate()}\t\t {myshower.get_end_plate()}\t\t {myshower.get_content()}\t\t {c0}\t {x}, {y}, {z}")
        # fill selection string (to be used later)
        cluster_sel = ""
        cluster_sel = myshower.get_fedra_selection(cluster_sel)
        d_cluster_sel = myshower.get_boxy_dict(i,d_cluster_sel)
        n_sum_bins += myshower.get_nclusters()
        for c in range(myshower.get_nclusters()):
            clus = myshower.get_cluster(c)
            x,y,z = myshower.get_emu_coord(clus)
            #h3.Fill(clus[2], clus[0], clus[1], clus[3])
            # fill 3d histo
            h3.Fill(clus[2], x, y, clus[3])
print("now we have the rough selection, we can copy the tree")
print("")
boxy_fname = f"{histo_dir}/boxy_{args.cellx}_{args.celly}_bin{args.binwidth}_minclus{args.min_clusters}_thres{args.threshold}_maxgap{args.maxgap}.root"

if is_mc:
    try:
        tfile_boxy = ROOT.TFile(boxy_fname)
        tree_sel = tfile_boxy.Get("tracks")
    except:
        print("Boxy tree does not exist. Creating one... (~6 minutes)")
        tfile_boxy = ROOT.TFile(boxy_fname,"recreate")
        tree_sel = fedra_ttree.CloneTree(0)
        nentries_fedra = fedra_ttree.GetEntries()
        for ievent, event in enumerate(fedra_ttree):
            save_event = False
            if ievent % 1000 == 0 :
                print(f"event {ievent}/{nentries_fedra}",end="\r")
            for cid in d_cluster_sel:
                if event.sf[0].X() < d_cluster_sel[cid][0] or event.sf[0].X() > d_cluster_sel[cid][3]:
                    continue
                if event.sf[0].Y() < d_cluster_sel[cid][1] or event.sf[0].Y() > d_cluster_sel[cid][4]:
                    continue
                if event.s[0].Plate() < d_cluster_sel[cid][2] or event.s[0].Plate() > d_cluster_sel[cid][5]:
                    continue
                save_event = True
                break
            if save_event:
                tree_sel.Fill()
        tree_sel.Write()

        # we have to run n_clusters times over a reduced ttree, hopefully will be faster
        print(f"TTree copied. Running over reduced tree, entries: {tree_sel.GetEntries()} / {fedra_ttree.GetEntries()}")

# output ttree
out_tfile = ROOT.TFile(f"{output_dir}/cluster_{args.cellx}_{args.celly}_bin{args.binwidth}_minclus{args.min_clusters}_thres{args.threshold}_maxgap{args.maxgap}.root","recreate")
out_ttree = ROOT.TTree("emu_match","Matching tree, from emu fedra")

#
emu_x_ini   = array('f',[0])
emu_x_end   = array('f',[0])
emu_y_ini   = array('f',[0])
emu_y_end   = array('f',[0])
emu_z_ini   = array('i',[0])
emu_z_end   = array('i',[0])
emu_z_max   = array('i',[0])
emu_cell_x  = array('i',[0])
emu_cell_y  = array('i',[0])
emu_binning = array('i',[0])
emu_brick   = array('i',[0])
emu_nclusters = array('i',[0])
emu_ntracks = array('f',[0])
emu_emcevt  = ROOT.std.vector('int')()
emu_emcevt_ntracks = ROOT.std.vector('int')()
emu_purity  = array('f',[0])

emu_cell_x[0]  = args.cellx
emu_cell_y[0]  = args.celly
emu_binning[0] = args.binwidth

clus_x = ROOT.std.vector('float')()
clus_y = ROOT.std.vector('float')()
clus_z = ROOT.std.vector('float')()
clus_ntracks = ROOT.std.vector('float')()

out_ttree.Branch("shw_nclusters", emu_nclusters, "emu_nclusters/I")
out_ttree.Branch("shw_ntracks", emu_ntracks, "emu_ntracks/F")
out_ttree.Branch("shw_x_ini", emu_x_ini, "emu_x_ini/F")
out_ttree.Branch("shw_x_end", emu_x_end, "emu_x_end/F")
out_ttree.Branch("shw_y_ini", emu_y_ini, "emu_y_ini/F")
out_ttree.Branch("shw_y_end", emu_y_end, "emu_y_end/F")
out_ttree.Branch("shw_z_ini", emu_z_ini, "emu_z_ini/I")
out_ttree.Branch("shw_z_end", emu_z_end, "emu_z_end/I")
out_ttree.Branch("shw_z_max", emu_z_max, "emu_z_max/I")
out_ttree.Branch("shw_cell_x", emu_cell_x, "emu_cell_x/I")
out_ttree.Branch("shw_cell_y", emu_cell_y, "emu_cell_y/I")
out_ttree.Branch("shw_binning", emu_binning, "emu_binning/I")
out_ttree.Branch("shw_brick", emu_brick, "emu_brick/I")
out_ttree.Branch("shw_eMCEvt",emu_emcevt)
out_ttree.Branch("shw_eMCEvt_ntracks",emu_emcevt_ntracks)
out_ttree.Branch("shw_purity",emu_purity, "emu_purity/F")
out_ttree.Branch("clus_x",clus_x)
out_ttree.Branch("clus_y",clus_y)
out_ttree.Branch("clus_z",clus_z)
out_ttree.Branch("clus_ntracks",clus_ntracks)

for i,c in enumerate(cell_xyz_clusters):

    emu_x_ini[0] = 0.0
    emu_x_end[0] = 0.0
    emu_y_ini[0] = 0.0
    emu_y_end[0] = 0.0
    emu_z_ini[0] = 0
    emu_z_end[0] = 0
    emu_z_max[0] = 0
    emu_nclusters[0] = 0
    emu_ntracks[0] = 0
    emu_emcevt.clear()
    emu_emcevt_ntracks.clear()
    emu_purity[0] = 0.0
    clus_x.clear()
    clus_y.clear()
    clus_z.clear()
    clus_ntracks.clear()
    
    myshower = cluster_shower(c,args.cellx, args.celly, args.binwidth)
    if myshower.get_nclusters() >= args.min_clusters:
        emu_nclusters[0] = myshower.get_nclusters()
        emu_ntracks[0] = myshower.get_content()
        emu_z_ini[0] = myshower.get_start_plate()
        emu_z_end[0] = myshower.get_end_plate()
        emu_z_max[0] = myshower.get_max_plate()
        emu_x_ini[0], emu_y_ini[0] = myshower.get_xy_center_at_plate(myshower.get_start_plate())
        emu_x_end[0], emu_y_end[0] = myshower.get_xy_center_at_plate(myshower.get_end_plate())
        
        for i in range(myshower.get_nclusters()):
            clus = myshower.get_cluster(i)
            clus_x.push_back(clus[0])
            clus_y.push_back(clus[1])
            clus_z.push_back(clus[2])
            clus_ntracks.push_back(clus[3])

        if is_mc:
            print(f"======================= SHOWER {i} ===========================")
            d_mcevt = {}
            l_trackcounter = 0
            print(f"Copying cluster {i} tree... (~?? minutes). Cluster selection:")
            print(f"TTree copied. Running over sel tree, entrie:s {tree_sel.GetEntries()} / {tree_sel.GetEntries()}")
            for ievent, event in enumerate(tree_sel):
                if myshower.is_xyz_in_cluster(event.sf[0].X(), event.sf[0].Y(), event.s[0].Plate()):
                    l_trackcounter += 1
                    if event.s[0].MCEvt() not in d_mcevt:
                        d_mcevt[event.s[0].MCEvt()] = 1
                    else:
                        d_mcevt[event.s[0].MCEvt()] += 1

            print("mc tracks in shower",l_trackcounter)
            print("\tnbins:",myshower.get_content())
            # sort d_mcevt from high to low value
            d_mcevt_sorted = sorted(d_mcevt.items(), key=lambda item: item[1], reverse=True)
            d_mcevt = dict(d_mcevt_sorted)
            print("mcevt dictionary:",d_mcevt)
            mcevt = max(d_mcevt, key=d_mcevt.get)
            print("event with most tracks:",mcevt)
            purity = d_mcevt[mcevt]/float(l_trackcounter)
            print("Purity",purity)
            emu_purity[0] = purity
            for mcevt in d_mcevt:
                emu_emcevt.push_back(mcevt)
                emu_emcevt_ntracks.push_back(d_mcevt[mcevt])
        out_ttree.Fill()
# finished looping over clusters
# write output root file
print("Finished looping over clusters. Writing output root file")
out_ttree.Write()
    
            
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas()
c.SetCanvasSize(1000,1000)
out_pdf = f"{work_dir}/emu/img/clusters_{args.cellx}_{args.celly}_bin{args.binwidth}_minclus{args.min_clusters}_thres{args.threshold}_maxgap{args.maxgap}.pdf"
c.Print(f"{out_pdf}[")

# 3d
h3.Draw("box2z")
c.Print(f"{out_pdf}")

# xz
c.SetTheta(90.0)
c.SetPhi(0.0)
c.Print(f"{out_pdf}")

# yz
c.SetTheta(0.0)
c.SetPhi(0.0)
c.Print(f"{out_pdf}")

c.Print(f"{out_pdf}]")
