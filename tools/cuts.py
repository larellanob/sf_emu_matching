def veto_cut(Digi_MuFilterHits):
    Veto1BarHit = False
    Veto2BarHit = False
    for aHit in Digi_MuFilterHits:
        if Veto1BarHit and Veto2BarHit:
            break
        if not aHit.isValid():
            continue
        if aHit.GetSystem() == 1:
            if aHit.GetPlane() == 0:
                Veto1BarHit = True
            if aHit.GetPlane() == 1:
                Veto2BarHit = True
    if Veto1BarHit == False or Veto2BarHit == False:
        return False
    else:
        return True

def scifi_gt_n_hits(Digi_ScifiHits,n):
    # scifi cuts:
    ## single hit on station 1
    ## large variation in n_hits from station 1 to station 2
    if len(Digi_ScifiHits) < n:
        return False
    else:
        return True

