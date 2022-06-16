#!/usr/bin/env python

import sys
import os,re
import datetime
import numpy
import multiprocessing

import awkward as ak
from coffea import nanoevents
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import processor, hist
from coffea import lookup_tools
from coffea import analysis_tools

from typing import Iterable, Callable, Optional, List, Generator, Dict, Union
import collections

import ROOT


## IMPLEMENT SELECTION AND HISTOS
import Selection as selection
import Histograms as Histograms


class Analysis(processor.ProcessorABC):
    def __init__(self,dataset,pu_data): ####run,isData, pu_data, pu_mc):
        self.run = dataset.run
        self.year = dataset.run[:4]
        self.isData = dataset.isData

        if self.isData:
            parsePileUpJSON2(self.year)

        if not self.isData:
            self.pu_data = pu_data
            self.pu_mc   = dataset.getPileup()
            self.pu_weight = PileupWeight.PileupWeight(self.pu_data,self.pu_mc)

        self.lumimask = LumiMask.LumiMask(dataset)

        self.histo = Histograms.AnalysisHistograms(self.isData)

        self.counter = Counter.Counters()
        self.counter.book(dataset.histograms["skimCounter"])

        self.book_histograms()
        self.first = True

    def book_histograms(self):
        self.histograms = {}
        self.histograms.update(self.counter.get())
        self.histograms.update(self.histo.book())

        if not self.isData:
            self.addHistogram('pu_orig', 100, 1, 100)
            self.cloneHistogram('pu_orig', 'pu_corr')
            self.cloneHistogram('pu_orig', 'pu_data')
        print("Booked",len(self.histograms),"histograms")
        self._accumulator = processor.dict_accumulator(self.histograms)

    def addHistogram(self, name, nbins, binmin, binmax):
        self.histograms[name] = hist.Hist(
            "",
            hist.Bin("value", name, nbins, binmin, binmax)
        )

    def cloneHistogram(self, nameOrig, nameClone):
        self.histograms[nameClone] = self.histograms[nameOrig].copy()

    def getArrays(self, histo):
        x = []
        axis  = histo.axis()
        edges = axis.edges()
        for i in range(0,len(edges)-1):
            bincenter = edges[i] + 0.5*(edges[i+1]-edges[i])
            x.append(bincenter)
        y = histo.values()
        return ak.from_numpy(numpy.array(x)),y

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        out = self.accumulator.identity()

        self.counter.setAccumulatorIdentity(out)
        self.counter.setSkimCounter()

        # Weights
        eweight = analysis_tools.Weights(len(events),storeIndividual=True)
        if not self.isData:
            genw = events['genWeight']
            eweight.add('pileup',genw)

            pu = self.pu_weight.getWeight(events.Pileup.nTrueInt)
            eweight.add('pileup',pu)
        events["weight"] = eweight.weight()
        self.counter.increment('all events',events)

        # if self.isData:
        #     events = events[selection.lumimask(events,self.lumimask)]
        # self.counter.increment('JSON filter',events)

    
        ## Constructing variables

        
        # Plotting
        out = self.histo.fill(events,out)
        
        return out

    def postprocess(self, accumulator):
        self.counter.printCounters(accumulator)
        return accumulator

# def usage():
    # print
    # print( "### Usage:  ",os.path.basename(sys.argv[0]),"<multicrab skim>" )
    # print

def main():

    # multicrabdir = os.path.abspath(sys.argv[1])
    # if not os.path.exists(multicrabdir) or not os.path.isdir(multicrabdir):
    #     usage()
    #     sys.exit()
    # year = multicrabdatasets.getYear(multicrabdir)

    # starttime = datetime.datetime.now().strftime("%Y%m%dT%H%M")

    # datasets = multicrabdatasets.getDatasets(multicrabdir,whitelist=whitelist,blacklist=blacklist)
    # pileup_data = multicrabdatasets.getDataPileupMulticrab(multicrabdir)
    # lumi = multicrabdatasets.loadLuminosity(multicrabdir,datasets)

    # outputdir = os.path.basename(os.path.abspath(multicrabdir))+"_processed"+starttime+lepton
    # if not os.path.exists(outputdir):
    #     os.mkdir(outputdir)

    # import time
    # t0 = time.time()

    # print("Number of cores used",MAX_WORKERS)
    # if len(datasets) == 0:
    #     print("No datasets to be processed")
    #     print("  whitelist:",whitelist)
    #     print("  blacklist:",blacklist)
    #     sys.exit()

    # for i,d in enumerate(datasets):
    #     print("Dataset %s/%s %s"%(i+1,len(datasets),d.name))
    #     t00 = time.time()
    #     subdir = os.path.join(outputdir,d.name)
    #     if not os.path.exists(subdir):
    #         os.mkdir(subdir)
    #         os.mkdir(os.path.join(subdir,"results"))

    #     samples = {d.name: d.getFileNames()}

    #     #job_executor = processor.FuturesExecutor(workers = MAX_WORKERS)
    #     job_executor = processor.IterativeExecutor()
    #     run = processor.Runner(
    #         executor = job_executor,
    #         schema=nanoevents.NanoAODSchema,
    #         chunksize = CHUNKSIZE,
    #         maxchunks = MAXCHUNKS
    #     )
    #     result = run(samples, 'Events', Analysis(d,pileup_data))
    #     """
    #     result = processor.run_uproot_job(
    #         samples,
    #         "Events",
    #         Analysis(d,pileup_data), ####d.run,d.isData,pileup_data,d.getPileup()),
    #         processor.iterative_executor,
    #         {"schema": NanoAODSchema},
    #     )
    #     """
    #     t01 = time.time()
    #     dt0 = t01-t00
    #     print("Processing time %s min %s s"%(int(dt0/60),int(dt0%60)))

    #     fOUT = ROOT.TFile.Open(os.path.join(subdir,"results","histograms.root"),"RECREATE")
    #     fOUT.cd()
    #     fOUT.mkdir("configInfo")
    #     fOUT.cd("configInfo")

    #     days = ["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"]
    #     now = datetime.datetime.now()
    #     m = "produced: %s %s"%(days[now.weekday()],now)
    #     timestamp = ROOT.TNamed(m,"");
    #     timestamp.Write()

    #     gitCommit = aux.execute("git rev-parse HEAD")[0]
    #     gc = "git commit:"+gitCommit
    #     gitcommit = ROOT.TNamed(gc,gitCommit);
    #     gitcommit.Write()

    #     h_lumi = ROOT.TH1D("lumi","",1,0,1)
    #     h_lumi.SetBinContent(1,d.lumi)
    #     h_lumi.Write()

    #     h_isdata = ROOT.TH1D("isdata","",1,0,1)
    #     h_isdata.SetBinContent(1,int(d.isData))
    #     h_isdata.Write()

    #     fOUT.cd()
    #     fOUT.mkdir("analysis")

    #     for key in result.keys():
    #         #print("check keys",key)
    #         if 'counter' in key or key in ['pu_orig','pu_data','pu_corr']:
    #             fOUT.cd("configInfo")
    #         else:
    #             fOUT.cd("analysis")
    #         histo = hist2root.convert(result[key]).Clone(key)
    #         histo.Write()
    #         fOUT.cd()
    #     fOUT.Close()
    #     dt1 = time.time()-t01
    #     print("Converting hist2root time %s min %s s"%(int(dt1/60),int(dt1%60)))

    # dt = time.time()-t0

    # print("Total processing time %s min %s s"%(int(dt/60),int(dt%60)))
    # print("output in",outputdir)

if __name__ == "__main__":

    main()
    #os.system("ls -lt")
    #os.system("pwd")
