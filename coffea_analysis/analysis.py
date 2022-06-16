import uproot
import awkward as ak
from coffea import processor #,hist
import coffea.hist as hist
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
import matplotlib.pyplot as plt
from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)
#import cairo

class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(int),
            "pt": hist.Hist(
               "Events",
               hist.Cat("dataset", "Dataset"),
               hist.Bin("pt", "$pT_{\mu\mu}$ [GeV]", 20, 20, 500),
            ),
        })
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        
        dataset = events.metadata['dataset']
        muons = ak.zip({
            "pt": events.Muon_pt,
            "eta": events.Muon_eta,
            "phi": events.Muon_phi,
            "mass": events.Muon_mass,
            "charge": events.Muon_charge,
        }, with_name="PtEtaPhiMCandidate")

        cut = (ak.num(muons) > 0) & (ak.sum((muons.pt > 50),axis=1)>=1)
        muons = muons[cut]
        print('muons =',len(muons))
        output['sumw'][dataset] += len(events) 
        output['pt'].fill(
            dataset=dataset,
            pt=ak.sum(muons.pt,axis=1)
        )

        return output

    def postprocess(self, accumulator):
        return accumulator

import uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema

uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource
filename = "~/Downloads/SUSY_T5tttt_CMSNANOAOD.root" #"root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/Run2012B_DoubleMuParked.root"
file = uproot.open(filename)
events = NanoEventsFactory.from_root(
    file,
    #entry_stop=10000,
    metadata={"dataset": "DoubleMuon"},
    schemaclass=BaseSchema,
).events()

print('events =',len(events))
p = MyProcessor()
out = p.process(events)

hist.plot1d(out['pt'])
plt.show()
