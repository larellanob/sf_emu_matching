import os
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gInterpreter.ProcessLine('#include "sndGeometryGetter.h"')
from ROOT.snd import analysis_tools as ana

geo_file = "/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V4_2022.root"
SciFiDet, MuFilter = ana.GetGeometry(geo_file)

#f = ROOT.TFile("scifi_selection/scifi_selection_R4619-4619_NE-1.root")
f = ROOT.TFile("scifi_tracking/scifi_tracks_R4619_ev18860895.root")
#f = ROOT.TFile("scifi_eventfilter/scifi_selection_R4619_ev18860895_bigshowers.root")
#tree = f.Get("rawConv")

# configurations for SND-plane analysis tools
myopt = ROOT.snd.Configuration.Option.ti18_2022_2023
mySNDConfiguration = ROOT.snd.Configuration(myopt,
                                            SciFiDet,
                                            MuFilter)

tree_hits = f.Get("rawConv").GetTree()
tree_trks = f.Get("ht_scifi").GetTree()

# extrapolate track to a z position
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


# histograms


h_sh_ini = ROOT.TH1F("h_sh_ini","Shower start;Station;n showers",6,-0.5,5.5)
h_sh_end = ROOT.TH1F("h_sh_end","Shower end;Station;n showers",6,-0.5,5.5)
h_sh_len = ROOT.TH1F("h_sh_len","Shower length;Station;n showers",6,-0.5,5.5)

h_sh_pos = ROOT.TH2F("h_sh_pos",
                  "Shower position in station 2;x (cm);y (cm)",
                  100, -50,-5,
                  100, 10, 55
                  )

h_diff = ROOT.TH2F("h_diff",
                   "(tr-sh) Position;x (cm);y (cm)",
                   100, -10,10,
                   100, -10,10
                   )
h_diff_zoom = ROOT.TH2F("h_diff_zoom",
                        "(tr-sh) Position (zoomed in);x (cm);y (cm)",
                        100, -1,1,
                        100, -1,1
                        )

h_ang = ROOT.TH1F("h_ang",
                  "Angle between track and shower;Angle (rad);track-shower pairs",
                  158,-0.1,1.57)

h_ang_1 = ROOT.TH1F("h_ang_1",
                     "Angle between track and shower (single station);Angle (rad);track-shower pairs",
                     158,-0.1,1.57)

h_ang_23 = ROOT.TH1F("h_ang_23",
                     "Angle between track and shower (nstations >1);Angle (rad);track-shower pairs",
                     158,-0.1,1.57)

h_ang_for = ROOT.TH1F("h_ang_for",
                      "Forward angle between track and shower;Angle (rad);track-shower pairs",
                      210,-0.01,0.2)
                  
h_off_ang_pos = ROOT.TH2F("h_off_ang_pos",
                          ";Angle (rad);Pos. diff in quadrature",
                          100,-0.1,1.0,
                          100, 0.0,3.0
                          )

h_sh_ang_xz_1 = ROOT.TH1F("h_sh_ang_xz_1",
                        "SciFi showers (single station);Angle xz (rad);N showers",
                        158,-1.57,1.57)
h_sh_ang_yz_1 = ROOT.TH1F("h_sh_ang_yz_1",
                        "SciFi showers (single station);Angle yz (rad);N showers",
                        158,-1.57,1.57)
h_sh_ang_xz_2 = ROOT.TH1F("h_sh_ang_xz_2",
                        "SciFi showers (multiple stations);Angle xz (rad);N showers",
                        158,-1.57,1.57)
h_sh_ang_yz_2 = ROOT.TH1F("h_sh_ang_yz_2",
                        "SciFi showers (multiple stations);Angle yz (rad);N showers",
                        158,-1.57,1.57)

h_centroid_1 = ROOT.TH2F("h_centroid_1",
                         "Centroid of single-station showers (station 2);x (cm);y (cm)",
                         100, -50,-5,
                         100, 10, 55)
h_centroid_2 = ROOT.TH2F("h_centroid_2",
                         "Centroid of multi-station showers (station 2);x (cm);y (cm)",
                         100, -50,-5,
                         100, 10, 55)
h_centroid_3 = ROOT.TH2F("h_centroid_3",
                         "Centroid of multi-station showers (station 3);x (cm);y (cm)",
                         100, -50,-5,
                         100, 10, 55)

h_centroid_x = ROOT.TH1F("h_centroid_x",
                         "Centroid of showers;x (cm);n showers",
                         100, -50,-5)
h_centroid_y = ROOT.TH1F("h_centroid_y",
                         "Centroid of showers;y (cm);n showers",
                         100, 10, 55)
h_centroid_z = ROOT.TH1F("h_centroid_z",
                         "Centroid of showers;z (cm);n showers",
                         120, 280,400)



n_shower_events = 0
b_test_only = False
for ievent, (ehits, etrks) in enumerate(zip(tree_hits,tree_trks)):
    if b_test_only and ievent > 20:
        break
    
    SciFiDet.InitEvent(ehits.EventHeader)
    MuFilter.InitEvent(ehits.EventHeader)


    # initialise scifi plane objects
    v_ScifiStations = ana.FillScifi(mySNDConfiguration,
                                    ehits.Digi_ScifiHits,
                                    SciFiDet)

    # need US for pos,mom of showers
    v_US = ana.FillUS(mySNDConfiguration,
                      ehits.Digi_MuFilterHits,
                      MuFilter)

    show_ini = ana.GetScifiShowerStart(v_ScifiStations)
    show_end = ana.GetScifiShowerEnd(v_ScifiStations)
    #if show_end <= 3: continue

    # the function to find center and direction uses the centroids. lets see how they look like
    # stations z
    stx = { 0:000.0 }
    sty = { 0:000.0 }
    stz = { 0:000.0 }
    #stz = { 0:000.0, 1:300.0, 2:313.00, 3:329.0, 4:350.0, 5:380.0 }
    for si in range(5): # si = station index
        v_ScifiStations[si].FindCentroid()
        stx[si+1] = v_ScifiStations[si].GetCentroid().X()
        sty[si+1] = v_ScifiStations[si].GetCentroid().Y()
        stz[si+1] = v_ScifiStations[si].GetCentroid().Z()
        pass
        v_US[si].FindCentroid()
        #ci = v_ScifiStations[si].GetCentroid()
        #print("scifi station",si+1,"centroid",f"({ci.x():.2f},{ci.y():.2f},{ci.z():.2f})")
        h_centroid_x.Fill(v_ScifiStations[si].GetCentroid().X())
        h_centroid_y.Fill(v_ScifiStations[si].GetCentroid().Y())
        h_centroid_z.Fill(v_ScifiStations[si].GetCentroid().Z())
    
        
    pos, mom = ana.GetShowerInterceptAndDirection(mySNDConfiguration,
                                                  v_ScifiStations,
                                                  v_US)
    '''
    # "proper" way used by simona
    sh_scifi_planes, sh_us_planes = ana.GetShoweringPlanes(v_ScifiStations, v_US)
    pos, mom = ana.GetShowerInterceptAndDirection(mySNDConfiguration,
                                                  sh_scifi_planes,
                                                  sh_us_planes)
    '''
    
    if (math.isnan(mom.x())):
        continue

    # the returned shower position is the intercept, i.e. at z=0
    # we consider the shower to be between its start and end planes
    # so we use the stz dict to get the post at the centroid of the
    # initial and end plane.

    sh_ini = extrapolation(pos,mom,stz[show_ini])
    sh_end = extrapolation(pos,mom,stz[show_end])

    h_sh_ini.Fill(show_ini)
    h_sh_end.Fill(show_end)
    h_sh_len.Fill(show_end-show_ini+1)
    
    n_shower_events += 1

    
    if (show_end-show_ini) == 0: # single station showers
        h_sh_ang_xz_1.Fill(math.atan(mom.X()/mom.Z()))
        h_sh_ang_yz_1.Fill(math.atan(mom.Y()/mom.Z()))
        h_centroid_1.Fill(stx[2],sty[2])
    elif (show_end-show_ini) > 0: # multiple station showers
        h_sh_ang_xz_2.Fill(math.atan(mom.X()/mom.Z()))
        h_sh_ang_yz_2.Fill(math.atan(mom.Y()/mom.Z()))
        h_centroid_2.Fill(stx[2],sty[2])
        h_centroid_3.Fill(stx[3],sty[3])
    
    if b_test_only:
        print("----------")
        print(ievent,"event")
        print("Stz",stz)
        print("scifi shower ini, end",show_ini, show_end)
        print(f"Shower momentum: {mom.x():.4f} {mom.y():.4f} {mom.z():.4f}")
        print(f"Shower origin  : {pos.x():.4f} {pos.y():.4f} {pos.z():.4f}")
        print(f"Shower ini ({show_ini}) : {sh_ini.x():.4f} {sh_ini.y():.4f} {sh_ini.z():.4f}")
        print(f"Shower end ({show_end}) : {sh_end.x():.4f} {sh_end.y():.4f} {sh_end.z():.4f}")
    for trk in range(etrks.ntracks):
        #tx, ty, tz = extrapolation(etrks.start[trk],etrks.mom[trk],pos.z())
        t_v = extrapolation(etrks.start[trk],etrks.mom[trk],stz[2])
        tx = t_v.x()
        ty = t_v.y()
        tz = t_v.z()
        if b_test_only:
            print(f"extrap ini: {tx:.2f} {ty:.2f} {tz:.2f}")
            print(f"Track origin   : {tx:.2f} {ty:.2f} {tz:.2f}")
            print(f"Track momentum : {etrks.mom[trk].x():.2f} {etrks.mom[trk].y():.2f} {etrks.mom[trk].z():.2f}")
        h_sh_pos.Fill(sh_ini.x(),sh_ini.y())
        h_diff.Fill(tx-sh_ini.x(), ty-sh_ini.y())
        h_diff_zoom.Fill(tx-sh_ini.x(), ty-sh_ini.y())

        # angle between them
        trk_mom_xyz = ROOT.Math.XYZVector(etrks.mom[trk].x(),etrks.mom[trk].y(),etrks.mom[trk].z())
        #ang = math.acos(etrks.mom[trk].Unit().Dot(mom.Unit()))
        ang = math.acos(trk_mom_xyz.Unit().Dot(mom.Unit()))
        h_ang.Fill(ang)
        if (show_end - show_ini) == 0:
            h_ang_1.Fill(ang)
        elif (show_end - show_ini) > 0:
            h_ang_23.Fill(ang)
        h_ang_for.Fill(ang)

        posdiff = math.sqrt((tx-sh_ini.x())**2+(ty-sh_ini.y())**2)
        h_off_ang_pos.Fill(ang,posdiff)
        
    
print("events:",ievent)    
print("number of shower events:",n_shower_events)



if not b_test_only:

    # we stack a couple of histograms
    ths_ang = ROOT.THStack("ths_ang","Angle between track and shower (stacked);Angle (rad);N shower-track pairs")
    ths_ang.Add(h_ang_1)
    ths_ang.Add(h_ang_23)
    h_ang_1.SetFillColor(3)
    h_ang_23.SetFillColor(4)
        
    
    fig_fname = f"scifi_showering/tracks_vs_showers.pdf"
    c = ROOT.TCanvas()
    c.Print(f"{fig_fname}[")

    leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
    
    # shower start and finish
    h_sh_ini.SetLineWidth(2)
    h_sh_end.SetLineWidth(2)
    h_sh_len.SetLineWidth(2)
    h_sh_ini.SetLineColor(1)
    h_sh_end.SetLineColor(2)
    h_sh_len.SetLineColor(3)

    h_sh_ini.Draw()
    h_sh_end.Draw("same")
    h_sh_len.Draw("same")
    leg.AddEntry(h_sh_ini)
    leg.AddEntry(h_sh_end)
    leg.AddEntry(h_sh_len)
    leg.Draw("same")
    c.Print(f"{fig_fname}")
    
    h_sh_pos.Draw("colz")
    c.Print(f"{fig_fname}")

    h_diff.Draw("colz")
    c.Print(f"{fig_fname}")

    h_diff_zoom.Draw("colz")
    c.Print(f"{fig_fname}")

    '''
    h_ang.Draw()
    c.Print(f"{fig_fname}")

    h_ang_1.Draw()
    c.Print(f"{fig_fname}")
    h_ang_23.Draw()
    c.Print(f"{fig_fname}")
    '''
    leg.Clear()
    leg.AddEntry(h_ang_1,"Single SciFi station shower")
    leg.AddEntry(h_ang_23,"Multiple SciFi station shower")
    ths_ang.Draw()
    leg.Draw("same")
    c.Print(f"{fig_fname}")

    h_ang_for.Draw()
    c.Print(f"{fig_fname}")

    h_off_ang_pos.Draw("colz")
    c.Print(f"{fig_fname}")

    h_sh_ang_xz_1.Draw()
    c.Print(f"{fig_fname}")
    h_sh_ang_yz_1.Draw()
    c.Print(f"{fig_fname}")

    h_sh_ang_xz_2.Draw()
    c.Print(f"{fig_fname}")
    h_sh_ang_yz_2.Draw()
    c.Print(f"{fig_fname}")

    h_centroid_1.Draw("colz")
    c.Print(f"{fig_fname}")
    h_centroid_2.Draw("colz")
    c.Print(f"{fig_fname}")
    h_centroid_3.Draw("colz")
    c.Print(f"{fig_fname}")

    h_centroid_x.Draw()
    c.Print(f"{fig_fname}")
    h_centroid_y.Draw()
    c.Print(f"{fig_fname}")
    h_centroid_z.Draw()
    c.Print(f"{fig_fname}")




    c.Print(f"{fig_fname}]")
