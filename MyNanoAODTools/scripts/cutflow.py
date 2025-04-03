import ROOT

# Enable multithreading for performance
ROOT.EnableImplicitMT()

# File and Tree Name (Update if needed)
file_path = "filteredNanoAOD/1/1/BB27DD3F-B227-4B49-8AA5-B595818B0E78_Skim.root"
tree_name = "Events"

# Load the ROOT file into an RDataFrame
df = ROOT.RDataFrame(tree_name, file_path)

# Define the cuts
NLEPTON_CUT = "nLeptons >= 3"
CUT_A = "nev_passA == 1"
CUT_B = "nev_passB == 1"
CUT_C = "nev_passC == 1"
CUT_D = "nev_passD == 1"
CUT_Z = "A_Zmass>70 && A_Zmass<110"

# Initial number of events
total_events = df.Count().GetValue()

# Apply nLepton cut
df_nlepton = df.Filter(NLEPTON_CUT)
events_nlepton = df_nlepton.Count().GetValue()

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
events_A_Z = df_A.Filter(CUT_Z).Count().GetValue()
events_B_Z = df_B.Filter(CUT_Z).Count().GetValue()
events_C_Z = df_C.Filter(CUT_Z).Count().GetValue()
events_D_Z = df_D.Filter(CUT_Z).Count().GetValue()

# Print the Cutflow Table
print("\nCUTFLOW TABLE:")
print(f"{'Selection':<20}{'Events':>10}")
print("-" * 40)
print(f"{'Total Events':<20}{total_events:>10}")
print(f"{'nLepton ≥ 3':<20}{events_nlepton:>10}")
print(f"{'Channel':<20}{'A':>10}{'B':>10}{'C':>10}{'D':>10}")
print(f"{'Events':<20}{events_A:>10}{events_B:>10}{events_C:>10}{events_D:>10}")
print(f"{'Z window':<20}{events_A_Z:>10}{events_B_Z:>10}{events_C_Z:>10}{events_D_Z:>10}")

