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
            "trigger_pass": 0,
            "min_leptons_pass": 0,
            "z_candidate_found": 0,
            "final_events": 0
        }


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("invariant_mass", "F")  # Define new branch
        self.out.branch("A_Zmass", "F") #Masa invariante para canal A=eeen
        self.out.branch("B_Zmass", "F") #Masa invariante para canal A=eemn
        self.out.branch("C_Zmass", "F") #Masa invariante para canal A=mmen
        self.out.branch("D_Zmass", "F") #Masa invariante para canal A=mmmn
        self.out.branch("A_Dr_Z", "F") #Delta R para candidatos a Z canal A
        self.out.branch("B_Dr_Z", "F") #Delta R para candidatos a Z canal B
        self.out.branch("C_Dr_Z", "F") #Delta R para candidatos a Z canal C
        self.out.branch("D_Dr_Z", "F") #Delta R para candidatos a Z canal D
        self.out.branch("A_Sum_pt", "F") #Suma pt de leptones canal A
        self.out.branch("B_Sum_pt", "F") #Suma pt de leptones canal B
        self.out.branch("C_Sum_pt", "F") #Suma pt de leptones canal C
        self.out.branch("D_Sum_pt", "F") #Suma pt de leptones canal D
        self.out.branch("A_Sum_mass", "F") #Suma de masa de leptones canal A
        self.out.branch("B_Sum_mass", "F") #Suma de masa de leptones canal B
        self.out.branch("C_Sum_mass", "F") #Suma de masa de leptones canal C
        self.out.branch("D_Sum_mass", "F") #Suma de masa de leptones canal D
        
        
        inputTree.SetBranchStatus("Electron_pdgId",1)
        inputTree.SetBranchStatus("Muon_pdgId",1)
        inputTree.SetBranchStatus("MET_phi",1)
        

    def analyze(self, event):
        
        self.cutflow["total_events"] += 1
        
        
        #-------------------HLT triggers
        #2016:
        #commonDefinitions = [
        #Define("passTrigger","(HLT_TkMu50 || HLT_Mu50 || HLT_Ele27_XPTight_Gsf || HLT_Photon175)"),
        #]
        
        #2017:
        #commonDefinitions = [
        #Define("passTrigger","(HLT_Mu50 || HLT_OldMu100 || HLT_TkMu100 || HLT_Ele35_WPTight_Gsf || HLT_Photon200)"),
        #]
        
        #2018:
        HLT_triggers = [
        "HLT_Mu50",
        "HLT_OldMu100",
        "HLT_TkMu100",
        "HLT_Ele32_WPTight_Gsf",
        "HLT_Photon200"
        ]
        
        #Se activan las ramas de los triggers antes de acceder a ellas
        for trigger in HLT_triggers:
            event._tree.SetBranchStatus(trigger, 1) #habilita la lectura de las ramas
            
        pass_trigger = False
        
        
        #Iterar sobre los triggers de la lista
        for trigger in HLT_triggers:
            #acceder a la rama del trigger
            trigger_value = getattr(event, trigger, 0) #Si el trigger no existe se toma valor como 0
            
            if trigger_value == 1: #Si el valor es 1, trigger activado
                pass_trigger = True #El evento pasa por que almenos uno esta activado
                break #sale del ciclo ya que se encuentra que el evento se acepta
        
        if not pass_trigger:
            return False #si ningun trigger pasa, se rechaza el evento
        self.cutflow["trigger_pass"] += 1
        
        
        #2022:
        #commonDefinitions = [
        #Define("passTrigger","(HLT_Mu50 || HLT_Ele35_WPTight_Gsf || HLT_Photon200)"),
        #]
        
        
        
        
        
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leptons = list(muons) + list(electrons)
        
        
        
        met_pt = event.MET_pt #este tipo de variables no es un vector y no puede manejarse
        #como los demas ya que solo son un dato por evento
        met_phi = event.MET_phi
        
        
        #-------------------------MET
        
        #MET_pt_threshold = 40
        #if met_pt < MET_pt_threshold:
        #    return False
        
        #-----------------
    
        
        #---------------ELECTRONES
        
        #good_electrons = sorted([ele for ele in electrons if ele.pt > 10 and abs(ele.eta) < 2.5], key=lambda x: x.pt, reverse=True)
        
        
        # Asegurarse de que el primer electron tenga pT > 50 y el segundo > 10
        #if len(good_electrons) >= 2:
        
            
            #if good_electrons[0].pt > 50 and good_electrons[1].pt > 10:
                # Solo usamos los dos electrones mas importantes
            #    e1, e2 = good_electrons[0], good_electrons[1]
                # Calcular la masa invariante de estos dos electrones
            #    invariant_mass = self.computeInvariantMass(e1, e2)
            #    self.histograms["invariant_mass"].Fill(invariant_mass)
                
                # Llenar los histogramas de pT
            #self.histograms["electron_pt"].Fill(e1.pt)
            #self.histograms["electron_pt"].Fill(e2.pt)
            
         #   return True
                
             #print(good_electrons) 
                
        #-----------------MUONES
        
        
        #good_muons = sorted([mu for mu in muons if mu.pt > 20 and abs(mu.eta) < 2.4], key=lambda x: x.pt, reverse=True)
        
        
        # Asegurarse de que el primer muon tenga pT > 70 y el segundo > 20
        #if len(good_muons) >= 2:
            
            #Verificar si las cargas de los muones son opuestas entre si
         #   for i in range(len(good_muons)):
          #      for j in range(i+1, len(good_muons)):
                
           #         if good_muons[i].charge != good_muons[j].charge:
        
        
        
        ###HACER UN PRINT 
            #if good_muons[0].pt > 70 and good_muons[1].pt > 20:
                # Solo usamos los dos muones mas importantes
            #    m1, m2 = good_muons[0], good_muons[1]
                # Calcular la masa invariante de estos dos electrones
                #invariant_mass = self.computeInvariantMass(e1, e2)
                #self.histograms["invariant_mass"].Fill(invariant_mass)  
                
                # Llenar los histogramas de pT
                #self.histograms["muon_pt"].Fill(m1.pt)
                #self.histograms["muon_pt"].Fill(m2.pt)
                
          #             return True
                
        #return False
        
        #return len(good_electrons) + len(good_muons) >= self.minLeptons
        
        if len(leptons) < self.minLeptons:
            return False #No pasa el corte de leptones
            
            
            
        self.cutflow["min_leptons_pass"] += 1 #Si pasa el corte de leptones,
        #se cuenta este paso
        
        #Se busca un candidato a Z
        best_pair, best_mass = self.findBestZCandidate(leptons)
        if not best_pair:
            
            return False #No se encontro un candidato Z
        
        
        
        self.cutflow["z_candidate_found"] += 1 #Se encuentra un candidato a Z
        
        passed_in_channel = False #Marca si el evento ha pasado por algun canal
        
        #if len(leptons) >= 3:
           #------------------CANAL A
        if len(electrons) >= 3:
            passed_in_channel = True
            print("Eventos con 3 electrones - Canal A")  # Imprimir que es el canal A
            print(f"Numero de electrones en Canal A: {len(electrons)}")  # Imprimir numero de electrones
            lepton3 = electrons[2] #este corresponde al tercer electron del canal A
            dr = self.dr_l1l2_Z(best_pair)
            total_pt = self.Pt_threshold(leptons)
            invariant_mass = self.computeInvariantMass(best_pair[0], best_pair[1])    
            W_mass = self.WMass(lepton3, met_pt, met_phi)
            total_mass = invariant_mass + W_mass
              
            self.out.fillBranch("A_Zmass", best_mass)
            self.out.fillBranch("A_Sum_mass", total_mass)
            self.out.fillBranch("A_Dr_Z", dr)
            self.out.fillBranch("A_Sum_pt", total_pt)
              
            self.cutflow["final_events"] += 1 #Se cuenta solo una vez si pasa por este canal
            #passed_in_channel = True
            print(f"Evento Canal A: Z_mass={best_mass}, Sum_mass={total_mass}, Dr_Z={dr}, Sum_pt={total_pt}")
    
              
           #------------------CANAL B
        if len(electrons) >= 2 and len(muons) >= 1 and not passed_in_channel:
            passed_in_channel = True
            print("Evento con 2 electrones y 1 muon - Canal B")  # Imprimir que es el canal B
            print(f"Numero de electrones en Canal B: {len(electrons)}")  # Imprimir numero de electrones
            print(f"Numero de muones en Canal B: {len(muons)}")  # Imprimir numero de muones
    
            lepton3 = muons[0]
            dr = self.dr_l1l2_Z(best_pair)
            total_pt = self.Pt_threshold(leptons)
            invariant_mass = self.computeInvariantMass(best_pair[0], best_pair[1])
            W_mass = self.WMass(lepton3, met_pt, met_phi)
            total_mass = invariant_mass + W_mass   
              
              
            self.out.fillBranch("B_Zmass", best_mass)
            self.out.fillBranch("B_Sum_mass", total_mass)
            self.out.fillBranch("B_Dr_Z", dr) 
            self.out.fillBranch("B_Sum_pt", total_pt)
              
            self.cutflow["final_events"] += 1 #El evento pasa lor cortes del Canal B
            
            print(f"Evento Canal B: Z_mass={best_mass}, Sum_mass={total_mass}, Dr_Z={dr}, Sum_pt={total_pt}")

        
           #------------------CANAL C
        if len(muons) >= 2 and len(electrons) >= 1 and not passed_in_channel:
            passed_in_channel = True
            print("Evento con 2 muones y 1 electron - Canal C")  # Imprimir que es el canal C
            print(f"Numero de electrones en Canal C: {len(electrons)}")  # Imprimir numero de electrones
            print(f"Numero de muones en Canal C: {len(muons)}")  # Imprimir numero de muones
    
            lepton3 = electrons[0]
            dr = self.dr_l1l2_Z(best_pair)
            total_pt = self.Pt_threshold(leptons)
            invariant_mass = self.computeInvariantMass(best_pair[0], best_pair[1])
            W_mass = self.WMass(lepton3, met_pt, met_phi)
            total_mass = invariant_mass + W_mass
                     
            self.out.fillBranch("C_Zmass", best_mass)
            self.out.fillBranch("C_Sum_mass", total_mass)
            self.out.fillBranch("C_Dr_Z", dr)
            self.out.fillBranch("C_Sum_pt", total_pt)
              
            self.cutflow["final_events"] += 1 #El evento pasa lor cortes del Canal C
            print(f"Evento Canal C: Z_mass={best_mass}, Sum_mass={total_mass}, Dr_Z={dr}, Sum_pt={total_pt}")

                    
        
           #------------------CANAL D
        if len(muons) >= 3 and not passed_in_channel:
            passed_in_channel = True
            print("Evento con 3 muones - Canal D")  # Imprimir que es el canal D
            print(f"Numero de muones en Canal D: {len(muons)}")  # Imprimir numero de muones
    
            lepton3 = muons[2]
            dr = self.dr_l1l2_Z(best_pair)
            total_pt = self.Pt_threshold(leptons)
            invariant_mass = self.computeInvariantMass(best_pair[0], best_pair[1])
            W_mass = self.WMass(lepton3, met_pt, met_phi)
            total_mass = invariant_mass + W_mass 
              
            self.out.fillBranch("D_Zmass", best_mass)
            self.out.fillBranch("D_Sum_mass", total_mass)
            self.out.fillBranch("D_Dr_Z", dr)
            self.out.fillBranch("D_Sum_pt", total_pt)
              
            self.cutflow["final_events"] += 1 #El evento pasa lor cortes del Canal D
            print(f"Evento Canal D: Z_mass={best_mass}, Sum_mass={total_mass}, Dr_Z={dr}, Sum_pt={total_pt}")
              
        
        if not passed_in_channel:
            return False
        
        self.cutflow["final_events"] += 1
        
        
        return True
              
        
           
        
        
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
        
    def findBestZCandidate(self, leptons):
        Z_MASS = 91.2  # Z boson mass in GeV
        #Finds the lepton pair with invariant mass closest to the Z boson mass"""
        best_pair = None
        best_mass = float("inf")
        min_diff = float("inf")

        for lepton1, lepton2 in itertools.combinations(leptons, 2):
            if lepton1.charge + lepton2.charge != 0:  # Require opposite charge
                continue
            if abs(lepton1.pdgId) != abs(lepton2.pdgId): 
                continue
            
            mass = self.computeInvariantMass(lepton1, lepton2)
            diff = abs(mass - Z_MASS)

            if diff < min_diff:
                min_diff = diff
                best_mass = mass
                best_pair = (lepton1, lepton2)
                
        if best_pair:        
           return best_pair, best_mass  
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
           
    def Pt_threshold(self, leptons):
        total_pt = 0
        for lepton in leptons:
            total_pt += lepton.pt #Suma el pt de cada lepton
        return total_pt
        
               
    def endJob(self):
        self.outputFile.cd()
        for hist in self.histograms.values():
            hist.Write()
            
        print("Cutflow:")
        for cut, count in self.cutflow.items():
            print(f"{cut}: {count}")
            
        self.outputFile.Close()
        
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



