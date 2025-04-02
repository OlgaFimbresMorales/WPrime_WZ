import ROOT

# Enable multithreading for performance
ROOT.EnableImplicitMT()

# File and Tree Name (Update if needed)
file_path = "filteredNanoAOD/1/1/BB27DD3F-B227-4B49-8AA5-B595818B0E78_Skim.root"
tree_name = "Events"

# Load the ROOT file into an RDataFrame
df = ROOT.RDataFrame(tree_name, file_path)

# Define the cuts
MET_CUT = "MET_pt > 40"
NLEPTON_CUT = "nLeptons >= 3"

CUT_A = "nev_passA == 1"
CUT_A_MZ = "A_Zmass>70 && A_Zmass<110"

# Initial number of events
total_events = df.Count().GetValue()

# Apply each cut sequentially

df_nlepton = df.Filter(NLEPTON_CUT)
events_nlepton = df_nlepton.Count().GetValue()

df_nlepton_passA = df_nlepton.Filter(CUT_A)
events_passA = df_nlepton_passA.Count().GetValue()

df_A_Zmass = df_nlepton_passA.Filter(CUT_A_MZ)
events_passA_MZ = df_A_Zmass.Count().GetValue()

# Cortes por canal

#df_mass = df_nlepton.Filter(MASS_CUT)
#events_mass = df_mass.Count().GetValue()

# Print the Cutflow Table
print("\nCUTFLOW TABLE:")
print(f"{'Cut':<25}{'Events':>10}")
print("-" * 35)
print(f"{'Total Events':<25}{total_events:>10}")
print(f"{'nLepton ≥ 3':<25}{events_nlepton:>10}")
print(f"{'Events in A':<25}{events_passA:>10}")
print(f"{'Z window':<25}{events_passA_MZ:>10}")


# # Optional: Save cutflow to a histogram
# hist = ROOT.TH1F("cutflow", "Cutflow;Cut Stage;Events", 4, 0, 4)
# hist.GetXaxis().SetBinLabel(1, "Total Events")
# hist.GetXaxis().SetBinLabel(2, "MET > 40")
# hist.GetXaxis().SetBinLabel(3, "nLepton ≥ 3")
# #hist.GetXaxis().SetBinLabel(4, "60 ≤ mass ≤ 120")

# hist.SetBinContent(1, total_events)
# hist.SetBinContent(2, events_met)
# hist.SetBinContent(3, events_nlepton)
# #hist.SetBinContent(4, events_mass)

# # Draw and save histogram
# canvas = ROOT.TCanvas("canvas", "Cutflow Histogram", 800, 600)
# hist.Draw()
# canvas.SaveAs("cutflow.png")  # Saves histogram as an image

# # Keep the script interactive (if running in a Python script)
# input("Press Enter to exit...")
