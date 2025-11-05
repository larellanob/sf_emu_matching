import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
ROOT.gInterpreter.ProcessLine('#include "sndSciFiTools.h"')
from ROOT.snd import analysis_tools as ana


#inp_filename = "scifi_eventfilter/scifi_selection_R4619_ev18860895.root"
inp_filename = "scifi_eventfilter/scifi_selection_R6640_ev369364945.root"
inp_filename_nodir = inp_filename.split("/")[-1]
output_file_id = inp_filename_nodir.replace("scifi_selection_","")
output_file_id = output_file_id.replace(".root","")
print("Output file tag:",output_file_id)
runnumber = int(output_file_id.split("_")[0].replace("R",""))
print("Run number parsed:",runnumber)

#geo_file = "/afs/cern.ch/work/d/dannc/public/MuonMatching/files/geofile_sndlhc_TI18_2022_NagoyaEmu.root"
geo_file = "/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V4_2022.root"
if runnumber > 5412:
    geo_file = "/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V3_2023.root"
SciFiDet, MuFilter = ana.GetGeometry(geo_file)

f = ROOT.TFile(inp_filename)

# tasks
ioman   = ROOT.FairRootManager.Instance()
ioman.SetTreeName("rawConv")

run     = ROOT.FairRunAna()
source  = ROOT.FairFileSource(f)
run.SetSource(source)
outFile = ROOT.TFile(f'scifi_tracking/myfile_run_scifi_tracking_{output_file_id}.root','RECREATE')
sink    = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)
# Don't use FairRoot's default event header settings
#run.SetEventHeaderPersistence(False)

#ioman.SetSource(source)
#ioman.SetSink(sink)
#avoid some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

# Add Hough Transform tracking
import SndlhcMuonReco
ht_sf_task = SndlhcMuonReco.MuonReco()
run.AddTask(ht_sf_task)
#ht_sf_task.SetParFile("{}/python/TrackingParams_V1_28November2022.xml".format(os.environ["SNDSW_ROOT"]))
ht_sf_task.SetParFile("{}/python/TrackingParams.xml".format(os.environ["SNDSW_ROOT"]))
ht_sf_task.SetHoughSpaceFormat("linearSlopeIntercept")
ht_sf_task.SetTrackingCase('passing_mu_Sf')
#ht_sf_task.Passthrough()
ht_sf_task.SetStandalone()

run.Init()

#run.Run(1,10)
#exit()
#outFile.Write()
#outFile.Close()

# histograms
h2 = ROOT.TH2F("h2",
               "Hough tracks through SciFi station 2 (0.25x0.25 cm^{2});#it{x} (cm);#it{y} (cm)",
               180,-50,-5,
               180,10,55
               )
h2_1cm = ROOT.TH2F("h2_1cm",
                   "Hough tracks through SciFi station 2 (1x1 cm^{2});#it{x} (cm);#it{y} (cm)",
                   45,-50,-5,
                   45,10,55
                   )
h2_5cm = ROOT.TH2F("h2_5cm",
                   "Hough tracks through SciFi station 2 (5x5 cm^{2});#it{x} (cm);#it{y} (cm)",
                   9,-50,-5,
                   9,10,55
                   )
h_angxz = ROOT.TH1F("h_angxz",
                    "XZ angle of Hough tracks;#theta_{XZ} (rad);N. tracks",
                    400,-0.02,0.02
                    )
h_angyz = ROOT.TH1F("h_angyz",
                    "YZ angle of Hough tracks;#theta_{YZ} (rad);N. tracks",
                    400,-0.02,0.02
                    )
h_theta = ROOT.TH1F("h_theta",
                    "#theta angle of Hough tracks;#theta (rad);N. tracks",
                    400,-0.02,0.02
                    )


# out tree
out_filename = f"scifi_tracking/scifi_tracks_{output_file_id}.root"
print(f"Using output filename: {out_filename}")
out_file = ROOT.TFile(out_filename,"RECREATE")

# copy other stuff from source file
bl     = f.Get("BranchList")
bl      .Write("BranchList", ROOT.TObject.kSingleKey)
tbl    = f.Get("TimeBasedBranchList")
tbl     .Write("TimeBasedBranchList", ROOT.TObject.kSingleKey)
ffh    = f.Get("FileHeader")
ffh     .Write("FileHeader")
cbmout = f.Get("cbmout")
cbmout  .Write("cbmout")
rawConv = f.Get("rawConv").CloneTree()
rawConv .Write("rawConv")

out_tree = ROOT.TTree("ht_scifi","Hough tracks in SciFi")




from array import array
#o_ntracks = ROOT.Int_t()
o_ntracks = array('i',[0])
o_start = ROOT.std.vector('TVector3')()
o_stop  = ROOT.std.vector('TVector3')()
o_mom   = ROOT.std.vector('TVector3')()
#std_vec_tv3 = ROOT.std.vector('TVector3')
#o_trk_pts   = std_vec_tv3()
o_trk_pts   = ROOT.std.vector(ROOT.std.vector('TVector3'))()
o_trk_type  = ROOT.std.vector('int')()
o_trk_flag  = ROOT.std.vector('bool')()
o_ndf       = ROOT.std.vector('int')()
o_chi2      = ROOT.std.vector('float')()
o_chi2ndf   = ROOT.std.vector('float')()
o_slopexz   = ROOT.std.vector('float')()
o_slopeyz   = ROOT.std.vector('float')()
o_anglexz   = ROOT.std.vector('float')()
o_angleyz   = ROOT.std.vector('float')()
o_theta     = ROOT.std.vector('float')()

out_tree.Branch("ntracks",o_ntracks,"ntracks/I")
out_tree.Branch("start", o_start,32000,0)
out_tree.Branch("stop", o_stop)
out_tree.Branch("mom", o_mom)
out_tree.Branch("points", o_trk_pts)
out_tree.Branch("type", o_trk_type)
out_tree.Branch("flag", o_trk_flag)
out_tree.Branch("chi2", o_chi2)
out_tree.Branch("ndf", o_ndf)
out_tree.Branch("chi2ndf", o_chi2ndf)
out_tree.Branch("slope_xz", o_slopexz)
out_tree.Branch("slope_yz", o_slopeyz)
out_tree.Branch("angle_xz", o_anglexz)
out_tree.Branch("angle_yz", o_angleyz)
out_tree.Branch("theta", o_theta)

# matching functions
import match

# from daniele
def getOriginAndDims(nodepath):
    nav = ROOT.gGeoManager.GetCurrentNavigator()
    nav.cd(nodepath)
    O = {'X':0, 'Y':0, 'Z':0}
    D = {'X':0, 'Y':0, 'Z':0}
    N = nav.GetCurrentNode()
    S = N.GetVolume().GetShape()
    D['X'], D['Y'], D['Z'] = S.GetDX(),S.GetDY(),S.GetDZ()
    O['X'], O['Y'], O['Z'] = S.GetOrigin()[0],S.GetOrigin()[1],S.GetOrigin()[2]
    OriginArray = array('d', O.values())
    OriginTrans = array('d', [0, 0, 0])
    nav.LocalToMaster(OriginArray, OriginTrans)
    O['X'], O['Y'], O['Z'] = OriginTrans[0],OriginTrans[1],OriginTrans[2]
    return O, D

# extrapolate track to station
def track_extrapolation(track,station):
    ini = track.getStart()
    fin = track.getStop()
    slope_xz = track.getSlopeXZ()
    slope_yz = track.getSlopeYZ()
    Scifiplane_path = '/cave_1/Detector_0/volTarget_1/ScifiVolume{0}_{0}000000/'.format(station)
    scifi_or,scifi_dim = getOriginAndDims(Scifiplane_path)
    scifi_z = scifi_or['Z']
    scifi_x = slope_xz*(scifi_z-fin.Z())+fin.X()
    scifi_y = slope_yz*(scifi_z-fin.Z())+fin.Y()
    #ini.Print()
    #fin.Print()
    #print("--------")
    return scifi_x, scifi_y, scifi_z
    
# benchmarking
from time import process_time
time_start = process_time()
time_total = process_time()

# event loop
passing_events = 0
passing_fraction = 0
speed = 0
eta = 99999.
total_tracks = 0
nentries = ioman.GetInTree().GetEntries()
to_loop = nentries
if to_loop == -1:
    to_loop = nentries
import math
for ievent, event in enumerate(ioman.GetInTree()):
    
    if ievent > to_loop: break
    if passing_events > 5000:
        print("Stopping at 5000 passing events!")
        break
    time_end = process_time()
    time_delta = time_end-time_start
    time_total += time_delta
    remaining = to_loop-ievent
    if ievent != 0:
        #speed = 1./float(time_delta)
        speed = ievent/float(time_total)
        eta = remaining/speed
        passing_fraction = float(passing_events)/float(ievent)
    estimate_final_passing = passing_fraction*to_loop
    
    
    print(f"Event {ievent}/{to_loop}. Passing: {passing_events}. Speed: {speed:.2f} (events/s). ETA: {eta/60.:.2f} (min)",end="\r")
    time_start = time_end

    
    
    
    SciFiDet.InitEvent(event.EventHeader)
    MuFilter.InitEvent(event.EventHeader)

    ht_sf_task.kalman_tracks.Clear()
    ht_sf_task.Exec(0)
    if (len(ht_sf_task.kalman_tracks)) != 0:
        #print("found",len(ht_sf_task.kalman_tracks),"hough tracks",ievent)
        pass
    myTracks = ht_sf_task.kalman_tracks

    o_start   .clear()
    o_stop    .clear()
    o_mom     .clear()
    o_trk_type.clear()
    o_trk_flag.clear()
    o_ndf     .clear()
    o_chi2    .clear()
    o_chi2ndf .clear()
    o_slopexz .clear()
    o_slopeyz .clear()
    o_anglexz .clear()
    o_angleyz .clear()
    o_theta   .clear()
    
    o_ntracks[0] = len(myTracks)
    #if len(myTracks) == 0: continue
    for i in range(len(myTracks)):
        total_tracks +=1
        #slope_xz = myTracks[i].getSlopeXZ()
        #slope_yz = myTracks[i].getSlopeYZ()
        #print("\tslope xz: {}".format(slope_xz))
        #print("\tslope yz: {}".format(slope_yz))
        angle_xz = myTracks[i].getAngleXZ()
        angle_yz = myTracks[i].getAngleYZ()
        h_angxz.Fill(angle_xz)
        h_angyz.Fill(angle_yz)
        #print("\tangle xz: {}".format(angle_xz))
        #print("\tangle yz: {}".format(angle_yz))
        # track extrapolation to station 2
        trk_x, trk_y, trk_z = track_extrapolation(myTracks[i],2)

        start = myTracks[i].getStart()
        o_start.push_back(start)
        #o_start.SetXYZ(start.X(),start.Y(),start.Z())
        stop = myTracks[i].getStop()
        o_stop.push_back(stop)
        #o_stop.SetXYZ(stop.X(),stop.Y(),stop.Z())
        mom = myTracks[i].getTrackMom() # TVector3 object
        o_mom.push_back(mom)
        #o_mom.SetXYZ(mom.X(),mom.Y(),mom.Z())
        theta = math.atan(math.sqrt(mom.X()*mom.X()+mom.Y()*mom.Y())/mom.Z())

        h_theta.Fill(theta)
        h2.Fill(trk_x,trk_y)
        h2_1cm.Fill(trk_x,trk_y)
        h2_5cm.Fill(trk_x,trk_y)
        
        o_trk_pts.clear()
        trk_pts_temp = myTracks[i].getTrackPoints()
        o_trk_pts.push_back(trk_pts_temp)
        #for trk in trk_pts_temp:
            #print(f"Printing track (out of {len(trk_pts_temp)})")
            #trk.Print()
            #o_trk_pts.push_back(trk)
        o_trk_type.push_back(myTracks[i].getTrackType())
        o_trk_flag.push_back(myTracks[i].getTrackFlag())
        o_chi2.    push_back(myTracks[i].getChi2())
        o_ndf.     push_back(myTracks[i].getNdf())
        o_chi2ndf. push_back(myTracks[i].getChi2Ndf())
        o_slopexz. push_back(myTracks[i].getSlopeXZ())
        o_slopeyz. push_back(myTracks[i].getSlopeYZ())
        o_anglexz. push_back(angle_xz)
        o_angleyz. push_back(angle_yz)
        o_theta.   push_back(theta)

        
        continue
        #print("x,y,z:",trk_x,trk_y,trk_z)
        if not match.is_in_b1(trk_x,trk_y):
            #print("unfavorable track!")
            #print("---------------------------------------")
            continue
        # now we loop over emulsion tracks!
        cell = match.which_cell(trk_x,trk_y)
        #print(cell)
        #print(match.is_in_cell(trk_x,trk_y,cell))
        #print(match.get_cell_borders(cell))
    # finished looping over tracks in the event
    out_tree.Fill()
    passing_events+=1
#out_raw.Write()
out_tree.Write()
print("passing_events",passing_events)
c = ROOT.TCanvas("c","my canvas",3000,3000)
fig_filename = f"scifi_tracking/scifi_tracks_{output_file_id}.pdf"
c.Print(f"{fig_filename}[")
h2.Draw("colz")
c.Print(f"{fig_filename}")
h2_1cm.Draw("colz")
c.Print(f"{fig_filename}")
h2_5cm.Draw("colz")
c.Print(f"{fig_filename}")
h_angxz.Draw()
c.Print(f"{fig_filename}")
h_angyz.Draw()
c.Print(f"{fig_filename}")
h_theta.Draw()
c.Print(f"{fig_filename}")
c.Print(f"{fig_filename}]")
print("total tracks:",total_tracks)
