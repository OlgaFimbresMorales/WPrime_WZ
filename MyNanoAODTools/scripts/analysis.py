#!/usr/bin/env python3

import os
import sys
import ROOT
import math
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object


class LeptonAnalysis(Module):
    def __init__(self, outputFile):
        self.minLeptons = 3
        self.outputFile = outputFile
        self.histograms = {}

    def beginJob(self):
        self.outputFile = ROOT.TFile(self.outputFile, "RECREATE")
        self.histograms["electron_pt"] = ROOT.TH1F("electron_pt", "Electron pT; pT (GeV); Events", 50, 0, 200)
        self.histograms["muon_pt"] = ROOT.TH1F("muon_pt", "Muon pT; pT (GeV); Events", 50, 0, 200)
        self.histograms["invariant_mass"] = ROOT.TH1F("invariant_mass", "Invariant Mass of Leading Electrons; Mass (GeV); Events", 50, 0, 200)
        #self.histograms["invariant_mass"] = ROOT.TH1F("invariant_mass", "Invariant Mass of Leading Electrons; Mass (GeV); Events", 50, 0, 200)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("invariant_mass", "F")  # Define new branch

    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        
        
        #---------------ELECTRONES
        
        good_electrons = sorted([ele for ele in electrons if ele.pt > 10 and abs(ele.eta) < 2.5], key=lambda x: x.pt, reverse=True)
        
        # Ordenar electrones por pT de mayor a menor
        good_electrons.sort(key=lambda x: x.pt, reverse=True)
        
        
        # Asegurarse de que el primer electron tenga pT > 50 y el segundo > 10
        if len(good_electrons) >= 2:
            if good_electrons[0].pt > 50 and good_electrons[1].pt > 10:
                # Solo usamos los dos electrones mas importantes
                e1, e2 = good_electrons[0], good_electrons[1]
                # Calcular la masa invariante de estos dos electrones
                invariant_mass = self.computeInvariantMass(e1, e2)
                self.histograms["invariant_mass"].Fill(invariant_mass)
                
                # Llenar los histogramas de pT
                self.histograms["electron_pt"].Fill(e1.pt)
                self.histograms["electron_pt"].Fill(e2.pt)
                
        #-----------------MUONES
        
        
        good_muons = [mu for mu in muons if mu.pt > 20 and abs(mu.eta) < 2.4]
        
        #Ordenas muones por pT de mayor a menor
        good_muons.sort(key=lambda x: x.pt, reverse=True)
        
        # Asegurarse de que el primer muon tenga pT > 70 y el segundo > 20
        if len(good_muons) >= 2:
            if good_muons[0].pt > 70 and good_muons[1].pt > 20:
                # Solo usamos los dos muones mas importantes
                m1, m2 = good_muons[0], good_muons[1]
                # Calcular la masa invariante de estos dos electrones
                #invariant_mass = self.computeInvariantMass(e1, e2)
                #self.histograms["invariant_mass"].Fill(invariant_mass)
                
                # Llenar los histogramas de pT
                self.histograms["muon_pt"].Fill(m1.pt)
                self.histograms["muon_pt"].Fill(m2.pt)
                
        
        return len(good_electrons) + len(good_muons) >= self.minLeptons
    
    def computeInvariantMass(self, lepton1, lepton2):
        e1, px1, py1, pz1 = self.getLorentzVector(lepton1)
        e2, px2, py2, pz2 = self.getLorentzVector(lepton2)
        mass2 = (e1 + e2) ** 2 - (px1 + px2) ** 2 - (py1 + py2) ** 2 - (pz1 + pz2) ** 2
        return math.sqrt(mass2) if mass2 > 0 else 0
    
    def getLorentzVector(self, lepton):
        e = math.sqrt(lepton.pt**2 * math.cosh(lepton.eta)**2 + 0.000511**2)  # Electron mass ~0.511 MeV
        px = lepton.pt * math.cos(lepton.phi)
        py = lepton.pt * math.sin(lepton.phi)
        pz = lepton.pt * math.sinh(lepton.eta)
        return e, px, py, pz
    
    def endJob(self):
        self.outputFile.cd()
        for hist in self.histograms.values():
            hist.Write()
        self.outputFile.Close()
        
        


# Read input file and Condor job ID arguments
if len(sys.argv) < 5:
    print("Usage: filterNanoAOD.py <input.root> <cluster_id> <process_id> <out_folder>")
    sys.exit(1)

inputFile  = sys.argv[1]
out_folder  = sys.argv[2]
process    = sys.argv[3]
outputDir  = f"filteredNanoAOD/{outfolder}/{process}"  # Corrected
#outputDir = "filteredNanoAOD"

outputDir = out_folder
os.makedirs(outputDir, exist_ok=True)

# Use Condor job IDs for unique histogram filenames
#histOutputFile = os.path.join(outputDir, f"histograms_{cluster_id}_{process_id}.root")
branchSelFile = "branchsel.txt"


p = PostProcessor(
    outputDir, [inputFile],
    cut=None,
    branchsel=branchSelFile,  # Keeps only selected branches
    outputbranchsel=None,  # Ensures all new branches are included
    modules=[LeptonAnalysis()]#histOutputFile)],
    noOut=False,
    justcount=False
)

p.run()


