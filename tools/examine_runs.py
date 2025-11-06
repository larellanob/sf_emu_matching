import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = False
ROOT.gInterpreter.ProcessLine('#include "sndTchainGetter.h"')
from ROOT.snd import analysis_tools as ana

import csv

runs_events = {}
rmin = 4575
rmax = 4854

#rmin = 4595
#rmax = 4610

snd_root_path = os.environ["SNDSW_ROOT"]
csv_file_path = f"{snd_root_path}/analysis/tools/data_paths.csv"

for run in range(rmin,rmax+1):
    ch = ana.GetTChain(csv_file_path,run)
    runs_events[run] = ch.GetEntries()

with open("run_database.csv",'w',newline='') as csvfile:
    writer = csv.writer(csvfile,delimiter=',')
    # header
    writer.writerow(["run","entries","eventfilter","tracks","showers","matched"])
    for run in range(rmin,rmax+1):
        print(f"{run}\t{runs_events[run]}")
        writer.writerow([run,runs_events[run],"","","",""])
