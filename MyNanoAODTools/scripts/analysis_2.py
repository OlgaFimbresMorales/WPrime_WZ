#!/usr/bin/env python3

import os
import sys
import ROOT
import math
import itertools
from array import array
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object


class LeptonAnalysis(Module):
    def __init__(self, outputFile):
        self.minLeptons = 3
        self.outputFile = outputFile
        self.histograms = {}
        self.cutflow = {}

    def beginJob(self):
        self.outputFile = ROOT.TFile(self.outputFile, "RECREATE")
        self.histograms["electron_pt"] = ROOT.TH1F("electron_pt", "Electron pT; pT (GeV); Events", 50, 0, 200)
        self.histograms["muon_pt"] = ROOT.TH1F("muon_pt", "Muon pT; pT (GeV); Events", 50, 0, 200)
        self.histograms["invariant_mass"] = ROOT.TH1F("invariant_mass", "Invariant Mass of Leading Electrons; Mass (GeV); Events", 50, 0, 200)
        
        self.cutflow = {
            "total_events": 0,
            #"trigger_pass": 0,
            "min_leptons_pass": 0,
            "z_candidate_found_A": 0,
            "z_candidate_found_B": 0,
            "z_candidate_found_C": 0,
            "z_candidate_found_D": 0,
            "final_events_A": 0,
            "final_events_B": 0
        }


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("invariant_mass", "F")  # Define new branch
        self.out.branch("nLeptons", "F")
        self.out.branch("nev_passA","I")
        self.out.branch("nev_passB","I")
        self.out.branch("nev_passC","I")
        self.out.branch("nev_passD","I")
        
        branches = ["Zmass", "Wmass", "Dr_Z", "Sum_pt", "Sum_mass", "nlep"]
        
        for channel in ['A', 'B', 'C', 'D']:
            self.out.branch(f"{channel}_Lep3W_pt", "F")
            for branch in branches:
                self.out.branch(f"{channel}_{branch}", "F")
            for i in range(1, 3):
                self.out.branch(f"{channel}_Lep{i}Z_pt", "F")
        
        
        


        inputTree.SetBranchStatus("event",1)
        inputTree.SetBranchStatus("Electron_pdgId",1)
        inputTree.SetBranchStatus("Electron_charge",1)
        inputTree.SetBranchStatus("Electron_pt",1)
        inputTree.SetBranchStatus("Electron_eta",1)
        inputTree.SetBranchStatus("Electron_phi",1)
        inputTree.SetBranchStatus("Muon_pt",1)
        inputTree.SetBranchStatus("Muon_eta",1)
        inputTree.SetBranchStatus("Muon_phi",1)
        inputTree.SetBranchStatus("Muon_charge",1)
        inputTree.SetBranchStatus("Muon_pdgId",1)
        inputTree.SetBranchStatus("Muon_mass",1)
        inputTree.SetBranchStatus("MET_phi",1)
        
        

    def analyze(self, event):
        
        self.cutflow["total_events"] += 1
        
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leptons = list(muons) + list(electrons)
        
        met_pt = event.MET_pt 
        met_phi = event.MET_phi
        
        #---------------ELECTRONES minimut cuts
        good_electrons = sorted([ele for ele in electrons if ele.pt > 10 and abs(ele.eta) < 2.5], key=lambda x: x.pt, reverse=True) 
                
        #---------------MUONES minimut cuts
        good_muons = sorted([mu for mu in muons if mu.pt > 20 and abs(mu.eta) < 2.4], key=lambda x: x.pt, reverse=True)
        
        
        good_leptons = good_electrons + good_muons

        self.out.fillBranch("nLeptons", len(good_leptons)) #number of leptons per event


        # Initialize all channel flags to zero
        self.out.fillBranch("nev_passA", 0)
        self.out.fillBranch("nev_passB", 0)
        self.out.fillBranch("nev_passC", 0)
        self.out.fillBranch("nev_passD", 0)
        
        if len(good_leptons) >= self.minLeptons: 

            if len(good_electrons) >= 3:  # Canal A  eee

                self.out.fillBranch("nev_passA", 1) #number of leptons per event
                
                best_pair, best_mass_Z = self.findBestZCandidate(good_electrons)
                self.out.fillBranch("A_Zmass", best_mass_Z)



            if len(good_electrons) >= 2 and len(good_muons)>=1:  # Canal B eemu

                self.out.fillBranch("nev_passB", 1) #number of leptons per event
                
                best_pair, best_mass_Z = self.findBestZCandidate(good_electrons)
                self.out.fillBranch("B_Zmass", best_mass_Z)



            if len(good_muons) >= 2 and len(good_electrons)>=1:  # Canal C   mumue

                self.out.fillBranch("nev_passC", 1) #number of leptons per event
                
                best_pair, best_mass_Z = self.findBestZCandidate(good_muons)
                self.out.fillBranch("C_Zmass", best_mass_Z)


            if len(good_muons) >= 3:  # Canal D  mumumu

                self.out.fillBranch("nev_passD", 1) #number of leptons per event
                
                best_pair, best_mass_Z = self.findBestZCandidate(good_muons)
                self.out.fillBranch("D_Zmass", best_mass_Z)

            
            

        return True      
        
#        if len(good_leptons) >= self.minLeptons:   # at least with 3 leptons

    #         #if passed_in_channel:
    #         #    return False
                
    #         passed_in_channel = True
                
    #         self.cutflow["z_candidate_found_A"] += 1

#                best_pair, best_mass_Z = self.findBestZCandidate(good_leptons)
#                if not best_pair:
#                    return False  # No se encontro un par Z
#                self.out.fillBranch("A_Zmass", best_mass_Z)
#              return True
#            return True
#       return True          
    #         lepton1, lepton2 = best_pair
    #         lepton1_pt = lepton1.pt
    #         lepton2_pt = lepton2.pt

    #         leptons = [lepton for lepton in good_leptons if lepton != lepton1 and lepton != lepton2]

    #         best_lep, best_mass_W = self.findBestWCandidate(leptons, met_pt, met_phi)
    #         if not best_lep:
    #             return False #no se encontro un lepton para W
            
    #         lepton3 = best_lep
    #         lepton3_pt = lepton3.pt

    #         dr = self.dr_l1l2_Z(best_pair)
    #         total_pt = lepton1_pt + lepton2_pt + lepton3_pt
    #         total_mass = best_mass_Z + best_mass_W
    #         nlep = len(good_leptons)

    #     # Almacenar las ramas para el canal A
    #          self.out.fillBranch("A_Zmass", best_mass_Z)
    #         self.out.fillBranch("A_Wmass", best_mass_W)
    #         self.out.fillBranch("A_Sum_mass", total_mass)
    #         self.out.fillBranch("A_Dr_Z", dr)
    #         self.out.fillBranch("A_Lep1Z_pt", lepton1_pt)
    #         self.out.fillBranch("A_Lep2Z_pt", lepton2_pt)
    #         self.out.fillBranch("A_Lep3W_pt", lepton3_pt)
    #         self.out.fillBranch("A_Sum_pt", total_pt)
    #         self.out.fillBranch("A_nlep", nlep)

    #         self.cutflow["final_events_A"] += 1
    #         return True  # Termina el analisis para el canal A

    # # **Canal B**
    #     elif len(good_electrons) >= 2 and len(good_muons) >= 1:  # Asegurate de que hay al menos 2 electrones y 1 muon
    #         #if passed_in_channel:
    #         #    return False
                
    #         passed_in_channel = True
                
    #         self.cutflow["z_candidate_found_B"] += 1

    #         best_pair, best_mass_Z = self.findBestZCandidate(good_leptons)
    #         if not best_pair:
    #             print("No se encontro par Z en Canal B")
    #             return False  # No se encontro un par Z

    #         lepton1, lepton2 = best_pair
    #         lepton1_pt = lepton1.pt
    #         lepton2_pt = lepton2.pt

    #         leptons = [lepton for lepton in good_leptons if lepton != lepton1 and lepton != lepton2]

    #         best_lep, best_mass_W = self.findBestWCandidate(leptons, met_pt, met_phi)
    #         if not best_lep:
    #             return False #no se encontro un lepton para W
    #         lepton3 = best_lep
    #         lepton3_pt = lepton3.pt

    #         dr = self.dr_l1l2_Z(best_pair)
    #         total_pt = lepton1_pt + lepton2_pt + lepton3_pt
    #         total_mass = best_mass_Z + best_mass_W
    #         nlep = len(good_leptons)

    #     # Almacenar las ramas para el canal B
    #         self.out.fillBranch("B_Zmass", best_mass_Z)
    #         self.out.fillBranch("B_Wmass", best_mass_W)
    #         self.out.fillBranch("B_Sum_mass", total_mass)
    #         self.out.fillBranch("B_Dr_Z", dr)
    #         self.out.fillBranch("B_Lep1Z_pt", lepton1_pt)
    #         self.out.fillBranch("B_Lep2Z_pt", lepton2_pt)
    #         self.out.fillBranch("B_Lep3W_pt", lepton3_pt)
    #         self.out.fillBranch("B_Sum_pt", total_pt)
    #         self.out.fillBranch("B_nlep", nlep)

    #         self.cutflow["final_events_B"] += 1
    #         return True  # Termina el analisis para el canal B

# Si el evento no pasa para ningun canal, no hacer nada o rechazarlo
#        print("Evento no paso por Canal A ni B")
#        return True
       
#    return True

    def etaphiplane(self, lepton1, lepton2):
        dr_etaphi = ()       
    
    def computeInvariantMass(self, lepton1, lepton2):
        #Obtener los 4-vectores de los dos leptones
        e1, px1, py1, pz1 = self.getLorentzVector(lepton1)
        e2, px2, py2, pz2 = self.getLorentzVector(lepton2)
        
        #Calcular la masa invariante (M^2 = E^2 - p^2)
        mass2 = (e1 + e2) ** 2 - (px1 + px2) ** 2 - (py1 + py2) ** 2 - (pz1 + pz2) ** 2
        return math.sqrt(mass2) if mass2 > 0 else 0
    
    def getLorentzVector(self, lepton):
        
        
        m_lepton = 0.000511 if abs(lepton.pdgId) == 11 else 0.105  # Electron: 0.511 MeV, Muon: 105 MeV
        e = math.sqrt(lepton.pt**2 * math.cosh(lepton.eta)**2 + 0.000511**2)  # Electron mass ~0.511 MeV
        
        #Componentes del momento
        px = lepton.pt * math.cos(lepton.phi)
        py = lepton.pt * math.sin(lepton.phi)
        pz = lepton.pt * math.sinh(lepton.eta)
        
        #Valor absoluto del momento
        #abs_p = lepton.pt * math.cosh(lepton.eta)
        return e, px, py, pz
        
    def findBestZCandidate(self, good_leptons):
        Z_MASS = 91.2  # Z boson mass in GeV
        #Finds the lepton pair with invariant mass closest to the Z boson mass"""
        best_pair = None
        best_mass_Z = float("inf")
        min_diff = float("inf")

        for lepton1, lepton2 in itertools.combinations(good_leptons, 2):
            if lepton1.charge + lepton2.charge != 0:  # Require opposite charge
                continue
            if abs(lepton1.pdgId) != abs(lepton2.pdgId): 
                continue
            
            mass = self.computeInvariantMass(lepton1, lepton2)
            diff = abs(mass - Z_MASS)

            if diff < min_diff:
                min_diff = diff
                best_mass_Z = mass
                best_pair = (lepton1, lepton2)
                
        if best_pair:        
           return best_pair, best_mass_Z  
        else:
           return None, 0.0
           
    def findBestWCandidate(self, leptons, met_pt, met_phi):
        W_MASS = 80.4 #W boson mass in GeV
        best_lep = None
        min_diff = float("inf")
        best_mass_W = float("inf")
        
        for lepton in leptons:
            mass = self.WMass(lepton, met_pt, met_phi)
            diff = abs(mass - W_MASS)
            
            if diff < min_diff:
                min_diff = diff
                best_mass_W = mass
                best_lep = lepton
                
        if best_lep:
            return best_lep, best_mass_W
        else:
            return None, 0.0
        
           
    def dr_l1l2_Z(self, best_pair):
        dr_max = 1.5
        if best_pair:
           lepton1, lepton2 = best_pair
           dr = math.sqrt((lepton1.phi - lepton2.phi)**2 + (lepton1.eta - lepton2.eta)**2) 
           return dr
        else:
           return 0
        
               
    def endJob(self):
        self.outputFile.Write()
        self.outputFile.Close()
        
        print("Cutflow:")
        for cut, count in self.cutflow.items():
            print(f"{cut}: {count}")
        
    def WMass(self, lepton, met_pt, met_phi):
        e, px, py, pz = self.getLorentzVector(lepton)
        
        #El MET puede tratarse como un neutrino con energia igual a MET_pt
        # y phi igual a MET_phi
        met_px = met_pt * math.cos(met_phi)
        met_py = met_pt * math.sin(met_phi)
        
        met_e = met_pt
        
        mass2 = (e + met_e) ** 2 - (px + met_px) ** 2 - (py + met_py) ** 2 - (pz) ** 2
        return math.sqrt(mass2) if mass2 > 0 else 0
        
    


# Read input file and Condor job ID arguments
if len(sys.argv) < 4:
    print("Usage: filterNanoAOD.py <input.root> <cluster_id> <process_id> <out_folder>")
    sys.exit(1)

inputFile  = sys.argv[1]
outfolder  = sys.argv[2]
process    = sys.argv[3]
outputDir  = f"filteredNanoAOD/{outfolder}/{process}"  # Corrected
#outputDir = "filteredNanoAOD"

#outputDir = outfolder
os.makedirs(outputDir, exist_ok=True)

# Use Condor job IDs for unique histogram filenames
histOutputFile = os.path.join(outputDir, f"histograms_{process}.root")
branchSelFile = "branchsel.txt"


p = PostProcessor(
    outputDir, [inputFile],
    cut=None,
    branchsel=branchSelFile,  # Keeps only selected branches
    outputbranchsel=None,  # Ensures all new branches are included
    modules=[LeptonAnalysis(histOutputFile)],
    noOut=False,
    justcount=False
)

p.run()



