def cortes_por_canales (electrons, muons):
    #---------------ELECTRONES
        
        good_electrons = sorted([ele for ele in electrons if ele.pt > 10 and abs(ele.eta) < 2.5], key=lambda x: x.pt, reverse=True)
        
        
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
        
        
        good_muons = sorted([mu for mu in muons if mu.pt > 20 and abs(mu.eta) < 2.4], key=lambda x: x.pt, reverse=True)
        
        
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