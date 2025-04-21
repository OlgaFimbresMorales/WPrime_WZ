import ROOT
import glob
import os
import sys
import subprocess
import cmsstyle as CMS

# Enable multithreading for performance
ROOT.EnableImplicitMT()

ROOT.gROOT.SetBatch(True)

# Recibir el argumento
#if len(sys.argv) < 2:
#    print("Falta el argumento del directorio con archivos skim.root")
#    sys.exit(1)

#input_dir = sys.argv[1]
# Enable multithreading for performance
ROOT.EnableImplicitMT()

# File and Tree Name (Update if needed)
file_path = [
            "filteredNanoAOD/WprimeToWZToWlepZlep_narrow_M600_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1_NANOAODSIM/0/BB27DD3F-B227-4B49-8AA5-B595818B0E78_Skim.root",
            "filteredNanoAOD/WprimeToWZToWlepZlep_narrow_M600_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1_NANOAODSIM/1/27F93223-A2A0-494B-BF3C-F77D4B7235BE_Skim.root",            "filteredNanoAOD/WprimeToWZToWlepZlep_narrow_M600_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1_NANOAODSIM/2/F2488ADF-1FBF-904C-89DF-091EFF3D3946_Skim.root",
            "filteredNanoAOD/WprimeToWZToWlepZlep_narrow_M600_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1_NANOAODSIM/3/626B5066-7071-E24E-ACF8-4967F03AF4AB_Skim.root",
            "filteredNanoAOD/WprimeToWZToWlepZlep_narrow_M600_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1_NANOAODSIM/4/28432F1A-1665-DF4E-B2C7-72BEF517856F_Skim.root",
            "filteredNanoAOD/WprimeToWZToWlepZlep_narrow_M600_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1_NANOAODSIM/5/3FA01DED-F59E-134C-9A71-55E759763C50_Skim.root"            
            ]

# Base directory for the ROOT files
#base_dir = "filteredNanoAOD/WZTo3LNu_mllmin01_NNPDF31_TuneCP5_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1_NANOAODSIM"

# Use glob to find all .root files in subdirectories 0 to 65
#file_path = []
#for root, dirs, files in os.walk(input_dir):
#    for file in files:
#        if file.endswith("Skim.root"):
#            file_paths.append(os.path.join(root, file))
#            
#if not file_paths:
#    sys.exit(1)


tree_name = "Events"

# Load the ROOT file into an RDataFrame
df = ROOT.RDataFrame(tree_name, file_path)

# Define the cuts
NLEPTON_CUT = "nLeptons >= 3"
CUT_A = "nev_passA == 1"
CUT_B = "nev_passB == 1"
CUT_C = "nev_passC == 1"
CUT_D = "nev_passD == 1"
CUT_Z_A = "A_Zmass>70 && A_Zmass<110"
CUT_Z_B = "B_Zmass>70 && B_Zmass<110"
CUT_Z_C = "C_Zmass>70 && C_Zmass<110"
CUT_Z_D = "D_Zmass>70 && D_Zmass<110"
CUT_DR_A = "A_Dr_Z < 1.5"
CUT_DR_B = "B_Dr_Z < 1.5"
CUT_DR_C = "C_Dr_Z < 1.5"
CUT_DR_D = "D_Dr_Z < 1.5"
CUT_Mass_A = "A_Sum_mass > 120"
CUT_Mass_B = "B_Sum_mass > 120"
CUT_Mass_C = "C_Sum_mass > 120"
CUT_Mass_D = "D_Sum_mass > 120"
CUT_PT_A = "A_Sum_pt > 110"
CUT_PT_B = "B_Sum_pt > 110"
CUT_PT_C = "C_Sum_pt > 110"
CUT_PT_D = "D_Sum_pt > 110"


def apply_cuts_and_generate_histograms(df, cut_name, hist_var, hist_name, x_title, n_bins, x_min, x_max):
# df: RDataFrame con los datos.
# cut_name: Corte que se va a aplicar.
# hist_var: Nombre de la variable que se va a usar para el histograma.
# hist_name: Nombre del histograma.
# x_title: Titulo del eje X.
# n_bins, x_min, x_max: Parametros del histograma.

    # Aplicar el corte al RDataFrame
    df_filtered = df.Filter(cut_name)
    
    # Generar el histograma
    hist = df_filtered.Histo1D((hist_name, f"{hist_name};{x_title};Events", n_bins, x_min, x_max), hist_var)
    hist_obj = hist.GetValue()
    
    # Estilizado CMS ROOT-style
    CMS.SetExtraText("")
    CMS.SetCmsText("Private work (CMS simulation)")
    CMS.SetCmsTextFont(52)
    CMS.SetCmsTextSize(0.75 * 0.76)
    CMS.SetLumi("")
    
    # Canvas CMS style
    canv = CMS.cmsCanvas('', 0, 1, 0, 1, '', '', square=CMS.kSquare, extraSpace=0.01, iPos=0) 
    hist_obj.Draw("HIST")
    canv.Draw()
    
    # Guardar imagen como PNG
    canv.SaveAs(f"/eos/user/o/olfimbre/NanoAOD_Filtered/histogramas/{hist_name}.png")
    canv.Close()
    
    return hist

# Initial number of events
total_events = df.Count().GetValue()

# Apply nLepton cut
df_nlepton = df.Filter(NLEPTON_CUT)
events_nlepton = df_nlepton.Count().GetValue()


channels = ["A","B","C","D"]
# Apply channel cuts
df_A = df_nlepton.Filter(CUT_A)
df_B = df_nlepton.Filter(CUT_B)
df_C = df_nlepton.Filter(CUT_C)
df_D = df_nlepton.Filter(CUT_D)

events_A = df_A.Count().GetValue()
events_B = df_B.Count().GetValue()
events_C = df_C.Count().GetValue()
events_D = df_D.Count().GetValue()

# Apply Z window cut
events_A_Z = df_A.Filter(CUT_Z_A).Count().GetValue()
events_B_Z = df_B.Filter(CUT_Z_B).Count().GetValue()
events_C_Z = df_C.Filter(CUT_Z_C).Count().GetValue()
events_D_Z = df_D.Filter(CUT_Z_D).Count().GetValue()


#Apply dr cut
events_A_Z_DR = df_A.Filter(CUT_Z_A + " && " + CUT_DR_A).Count().GetValue();
events_B_Z_DR = df_B.Filter(CUT_Z_B + " && " + CUT_DR_B).Count().GetValue();
events_C_Z_DR = df_C.Filter(CUT_Z_C + " && " + CUT_DR_C).Count().GetValue();
events_D_Z_DR = df_D.Filter(CUT_Z_D + " && " + CUT_DR_D).Count().GetValue();

#Apply mass threshold cut

events_A_Mass_Threshold = df_A.Filter(CUT_Z_A + " && " + CUT_DR_A + "&&" + CUT_Mass_A).Count().GetValue();
events_B_Mass_Threshold = df_B.Filter(CUT_Z_B + " && " + CUT_DR_B + "&&" + CUT_Mass_B).Count().GetValue();
events_C_Mass_Threshold = df_C.Filter(CUT_Z_C + " && " + CUT_DR_C + "&&" + CUT_Mass_C).Count().GetValue();
events_D_Mass_Threshold = df_D.Filter(CUT_Z_D + " && " + CUT_DR_D + "&&" + CUT_Mass_D).Count().GetValue();

#Apply pt threshold cut

events_A_Pt_Threshold = df_A.Filter(CUT_Z_A + " && " + CUT_DR_A + "&&" + CUT_Mass_A + "&&" + CUT_PT_A ).Count().GetValue();
events_B_Pt_Threshold = df_B.Filter(CUT_Z_B + " && " + CUT_DR_B + "&&" + CUT_Mass_B + "&&" + CUT_PT_B ).Count().GetValue();
events_C_Pt_Threshold = df_C.Filter(CUT_Z_C + " && " + CUT_DR_C + "&&" + CUT_Mass_C + "&&" + CUT_PT_C ).Count().GetValue();
events_D_Pt_Threshold = df_D.Filter(CUT_Z_D + " && " + CUT_DR_D + "&&" + CUT_Mass_D + "&&" + CUT_PT_D ).Count().GetValue();


# Apply cuts and generate histograms

# Apply Z window cut and generate histograms
hist_A_Zmass = apply_cuts_and_generate_histograms(df, CUT_Z_A, "A_Zmass", "hist_A_Zmass", "A Z mass [GeV]", 40, 70, 110)
hist_B_Zmass = apply_cuts_and_generate_histograms(df, CUT_Z_B, "B_Zmass", "hist_B_Zmass", "B Z mass [GeV]", 40, 70, 110)
hist_C_Zmass = apply_cuts_and_generate_histograms(df, CUT_Z_C, "C_Zmass", "hist_C_Zmass", "C Z mass [GeV]", 40, 70, 110)
hist_D_Zmass = apply_cuts_and_generate_histograms(df, CUT_Z_D, "D_Zmass", "hist_D_Zmass", "D Z mass [GeV]", 40, 70, 110)

# Apply DR cut and generate histograms
hist_A_DR = apply_cuts_and_generate_histograms(df, CUT_Z_A + " && " + CUT_DR_A, "A_Dr_Z", "hist_A_DR", "A ?R [GeV]", 50, 0, 2)
hist_B_DR = apply_cuts_and_generate_histograms(df, CUT_Z_B + " && " + CUT_DR_B, "B_Dr_Z", "hist_B_DR", "B ?R [GeV]", 50, 0, 2)
hist_C_DR = apply_cuts_and_generate_histograms(df, CUT_Z_C + " && " + CUT_DR_C, "C_Dr_Z", "hist_C_DR", "C ?R [GeV]", 50, 0, 2)
hist_D_DR = apply_cuts_and_generate_histograms(df, CUT_Z_D + " && " + CUT_DR_D, "D_Dr_Z", "hist_D_DR", "D ?R [GeV]", 50, 0, 2)

# Apply Mass cut and generate histograms
hist_A_Mass = apply_cuts_and_generate_histograms(df, CUT_Z_A + " && " + CUT_DR_A + " && " + CUT_Mass_A, "A_Sum_mass", "hist_A_Mass", "A Sum Mass [GeV]", 50, 0, 500)
hist_B_Mass = apply_cuts_and_generate_histograms(df, CUT_Z_B + " && " + CUT_DR_B + " && " + CUT_Mass_B, "B_Sum_mass", "hist_B_Mass", "B Sum Mass [GeV]", 50, 0, 500)
hist_C_Mass = apply_cuts_and_generate_histograms(df, CUT_Z_C + " && " + CUT_DR_C + " && " + CUT_Mass_C, "C_Sum_mass", "hist_C_Mass", "C Sum Mass [GeV]", 50, 0, 500)
hist_D_Mass = apply_cuts_and_generate_histograms(df, CUT_Z_D + " && " + CUT_DR_D + " && " + CUT_Mass_D, "D_Sum_mass", "hist_D_Mass", "D Sum Mass [GeV]", 50, 0, 500)

# Apply Pt cut and generate histograms
hist_A_Pt = apply_cuts_and_generate_histograms(df, CUT_Z_A + " && " + CUT_DR_A + " && " + CUT_Mass_A + " && " + CUT_PT_A, "A_Sum_pt", "hist_A_Pt", "A Sum Pt [GeV]", 50, 0, 500)
hist_B_Pt = apply_cuts_and_generate_histograms(df, CUT_Z_B + " && " + CUT_DR_B + " && " + CUT_Mass_B + " && " + CUT_PT_B, "B_Sum_pt", "hist_B_Pt", "B Sum Pt [GeV]", 50, 0, 500)
hist_C_Pt = apply_cuts_and_generate_histograms(df, CUT_Z_C + " && " + CUT_DR_C + " && " + CUT_Mass_C + " && " + CUT_PT_C, "C_Sum_pt", "hist_C_Pt", "C Sum Pt [GeV]", 50, 0, 500)
hist_D_Pt = apply_cuts_and_generate_histograms(df, CUT_Z_D + " && " + CUT_DR_D + " && " + CUT_Mass_D + " && " + CUT_PT_D, "D_Sum_pt", "hist_D_Pt", "D Sum Pt [GeV]", 50, 0, 500)

# Save all histograms to a ROOT file
output_file = ROOT.TFile("/eos/user/o/olfimbre/NanoAOD_Filtered/histogramas/histograms_after_cuts.root", "RECREATE")



# Zmass histograms
hist_A_Zmass.Write()
hist_B_Zmass.Write()
hist_C_Zmass.Write()
hist_D_Zmass.Write()

# DR histograms
hist_A_DR.Write()
hist_B_DR.Write()
hist_C_DR.Write()
hist_D_DR.Write()

# Mass histograms
hist_A_Mass.Write()
hist_B_Mass.Write()
hist_C_Mass.Write()
hist_D_Mass.Write()

# Pt histograms
hist_A_Pt.Write()
hist_B_Pt.Write()
hist_C_Pt.Write()
hist_D_Pt.Write()

output_file.Close()


# Print the Cutflow Table
print("\nCUTFLOW TABLE:")
print(f"{'Selection':<20}{'Events':>10}")
print("-" * 40)
print(f"{'Total Events':<20}{total_events:>10}")
print(f"{'nLepton = 3':<20}{events_nlepton:>10}")
print(f"{'Channel':<20}{'A':>10}{'B':>10}{'C':>10}{'D':>10}")
print(f"{'Events':<20}{events_A:>10}{events_B:>10}{events_C:>10}{events_D:>10}")
print(f"{'Z window':<20}{events_A_Z:>10}{events_B_Z:>10}{events_C_Z:>10}{events_D_Z:>10}")
print(f"{'dr cut':<20}{events_A_Z_DR:>10}{events_B_Z_DR:>10}{events_C_Z_DR:>10}{events_D_Z_DR:>10}")
print(f"{'Mass cut':<20}{events_A_Mass_Threshold:>10}{events_B_Mass_Threshold:>10}{events_C_Mass_Threshold:>10}{events_D_Mass_Threshold:>10}")
print(f"{'Pt cut':<20}{events_A_Pt_Threshold:>10}{events_B_Pt_Threshold:>10}{events_C_Pt_Threshold:>10}{events_D_Pt_Threshold:>10}")
