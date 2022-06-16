import uproot
import awkward as ak
from coffea import processor #,hist
import coffea.hist as hist
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
import matplotlib.pyplot as plt
# register our candidate behaviors
from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)
#import cairo

class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(float),
            "mass": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("mass", "$m_{\mu\mu}$ [GeV]", 60, 60, 120),
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

        cut = (ak.num(muons) == 2) & (ak.sum(muons.charge) == 0)
        # add first and second muon in every event together
        dimuon = muons[cut][:, 0] + muons[cut][:, 1]

        output["sumw"][dataset] += len(events)
        output["mass"].fill(
            dataset=dataset,
            mass=dimuon.mass,
        )

        return output

    def postprocess(self, accumulator):
        return accumulator

import uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema

# https://github.com/scikit-hep/uproot4/issues/122
uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

filename = "~/Downloads/SUSY_T5tttt_CMSNANOAOD.root" #"root://eospublic.cern.ch//eos/root-eos/cms_opendata_2012_nanoaod/Run2012B_DoubleMuParked.root"
file = uproot.open(filename)
events = NanoEventsFactory.from_root(
    file,
    entry_stop=10000,
    metadata={"dataset": "DoubleMuon"},
    schemaclass=BaseSchema,
).events()
p = MyProcessor()
out = p.process(events)
print(out)

muonPt = events.Muon_pt
muonEta = events.Muon_eta

# tree = file["Events"]

# # let's build the lepton arrays back into objects
# # in the future, some of this verbosity can be reduced
# arrays = {k.replace('Electron_', ''): v for k, v in tree.arrays(filter_name="Electron_*", how=dict).items()}
# electrons = ak.zip({'x': arrays.pop('Px'),
#                     'y': arrays.pop('Py'),
#                     'z': arrays.pop("Pz"),
#                     't': arrays.pop("E"),
#                     },
#                     with_name="LorentzVector"
# )


# arrays = {k.replace('Muon_', ''): v for k,v in tree.arrays(filter_name="Muon_*", how=dict).items()}
# muons = ak.zip({'x': arrays.pop('Px'),
#                 'y': arrays.pop('Py'),
#                 'z': arrays.pop("Pz"),
#                 't': arrays.pop("E"),
#                 },
#                 with_name="LorentzVector"
# )

# print("Avg. electrons/event:", ak.sum(ak.num(electrons))/tree.num_entries)
# print("Avg. muons/event:", ak.sum(ak.num(muons))/tree.num_entries)

lepton_kinematics = hist.Hist(
    "Events",
    hist.Cat("flavor", "Lepton flavor"),
    hist.Bin("pt", "$p_{T}$", 19, 10, 100),
    hist.Bin("eta", "$\eta$", [-2.5, -1.4, 0, 1.4, 2.5]),
)

## Pass keyword arguments to fill, all arrays must be flat numpy arrays
## User is responsible for ensuring all arrays have same jagged structure!
# lepton_kinematics.fill(
#     flavor="electron",
#     pt=ak.flatten(electrons.pt),
#     eta=ak.flatten(electrons.eta)
# )

lepton_kinematics.fill(
    flavor="muon",
    pt=ak.flatten(muonPt),
    eta=ak.flatten(muonEta)
)

## Now we can start to manipulate this single histogram to plot different views of the data
## here we look at lepton pt for all eta
lepton_pt = lepton_kinematics.integrate("eta")
x = [10,2,6,3,7,4,5,7,2,4,5]
print(ak.flatten(muonPt))
print(x)
plt.hist(x) #ak.flatten(muonPt))
plt.show()

#print(lepton_kinematics.pt)
#plt.hist(lepton_kinematics.pt)
#plt.show()

# ax = hist.plot1d(
#     lepton_pt,
#     overlay="flavor",
#     stack=True,
#     fill_opts={'alpha': .5, 'edgecolor': (0,0,0,0.3)}
# )

# ## all plot calls return the matplotlib axes object, from which
# ## you can edit features afterwards using matplotlib object-oriented syntax
# ## e.g. maybe you really miss '90s graphics...
# ax.get_legend().shadow = True
# plt.show()
