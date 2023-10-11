import pandas as pd


hits = pd.read_csv("/afs/cern.ch/user/r/rlindley/workspace/muon_reco/acc_dev/txtfiles/hits_pileup_allphi.txt",sep=",")
particles = pd.read_csv("/afs/cern.ch/user/r/rlindley/workspace/muon_reco/acc_dev/txtfiles/truthtracks_pileup_allphi.txt",sep=",")
merge = pd.merge(hits, particles, on="event")

merge.to_csv('merge_pileup_allphi.txt', sep=" ", index=False)
