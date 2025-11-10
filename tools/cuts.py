import ROOT

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

class Cutflow:
    l_cuts = {} # "cut name": nentries
    dataname = ""

    def __init__(self,cuts_string_separated_by_colons):
        self.l_cuts = {}
        cut_split = cuts_string_separated_by_colons.split(":")
        for cut in cut_split:
            self.l_cuts[cut] = 0

    def set_data_name(self,dataname):
        self.dataname = dataname
    
    def add_cut(self,cutname):
        self.l_cuts[cutname] = 0

    def get_cut_names(self):
        return self.l_cuts.keys()

    def get_cut_entries(self,cut):
        return self.l_cuts[cut]
    
    def Fill(self,cutname):
        self.l_cuts[cutname] += 1
    
    def set_entries(self,cutname,entries):
        self.l_cuts[cutname] = entries

    def check_cutflow(self):
        is_consistent = True
        l_keys = list(self.l_cuts.keys())
        prev_entry = self.l_cuts[l_keys[0]]
        for key in l_keys:
            entry = self.l_cuts[key]
            if entry > prev_entry:
                is_consistent = False
                break
        return is_consistent
    
    def get_histogram(self):
        h = ROOT.TH1F("h",f"Cutflow {self.dataname};;Passing events",len(self.l_cuts),0,len(self.l_cuts))
        xax = h.GetXaxis()
        for i in range(1,h.GetNbinsX()+1):
            # root bin index start from 1, python from 0
            xax.SetBinLabel(i,list(self.l_cuts.keys())[i-1])
            h.SetBinContent(i,list(self.l_cuts.values())[i-1])
        
        # consistency check, should be decreasing!
        if self.check_cutflow():
            return h
        else:
            print(f"ERROR: Cutflow {self.dataname} inconsistent, check histogram!!!!")
            h.SetTitle("INCONSISTENT Cutflow;Selection criterion;N entries")
            return h
