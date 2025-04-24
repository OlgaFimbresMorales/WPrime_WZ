import ROOT as r
from ROOT import RDataFrame
import cmsstyle as CMS
import glob

# Activar multithreading
r.ROOT.EnableImplicitMT()
r.gROOT.SetBatch(True)

def is_non_empty_root_file(file_path):
    # Intenta abrir el archivo y verificar si tiene el arbol "Events" con entradas
    try:
        f = r.TFile(file_path)
        tree = f.Get("Events")
        if tree and tree.GetEntries() > 0:
            return True
    except Exception as e:
        print(f"Error al abrir el archivo {file_path}: {e}")
    return False

# Estilo CMS
CMS.SetExtraText("Simulation")
CMS.SetLumi("")
CMS.SetEnergy("13")
CMS.ResetAdditionalInfo()

# Cargar archivos ROOT
signal_files = glob.glob("/eos/user/o/olfimbre/NanoAOD_Filtered/WprimeToWZToWlepZlep_narrow_M600*/*Skim.root", recursive=True)
background_files = glob.glob("/eos/user/o/olfimbre/NanoAOD_Filtered/WZTo3LNu_mllmin01*/*Skim.root", recursive=True)

signal_files = [f for f in signal_files if is_non_empty_root_file(f)]
background_files = [f for f in background_files if is_non_empty_root_file(f)]

# Cargar arboles
df_signal = RDataFrame("Events", signal_files)
df_background = RDataFrame("Events", background_files)

# Definicion de cortes y variables
CUTS = {
    "A": {
        "base": "nev_passA == 1 && nLeptons >= 3",
        "Z": "A_Zmass > 70 && A_Zmass < 110",
        "DR": "A_Dr_Z < 1.5",
        "Mass": "A_Sum_mass > 120",
        "Pt": "A_Sum_pt > 110",
        "vars": {
            "Zmass": ("A_Zmass", "A Z Mass [GeV]", 20, 50, 150),
            "DR": ("A_Dr_Z", "A dR", 25, 0, 2),
            "Sum_mass": ("A_Sum_mass", "A Sum Mass [GeV]", 25, 0, 700),
            "Sum_pt": ("A_Sum_pt", "A Sum pT [GeV]", 25, 0, 800)
        }
    },
    "B": {
        "base": "nev_passB == 1 && nLeptons >= 3",
        "Z": "B_Zmass > 70 && B_Zmass < 110",
        "DR": "B_Dr_Z < 1.5",
        "Mass": "B_Sum_mass > 120",
        "Pt": "B_Sum_pt > 110",
        "vars": {
            "Zmass": ("B_Zmass", "B Z Mass [GeV]", 20, 50, 150),
            "DR": ("B_Dr_Z", "B dR", 25, 0, 2),
            "Sum_mass": ("B_Sum_mass", "B Sum Mass [GeV]", 25, 0, 700),
            "Sum_pt": ("B_Sum_pt", "B Sum pT [GeV]", 25, 0, 800)
        }
    },
    "C": {
        "base": "nev_passC == 1 && nLeptons >= 3",
        "Z": "C_Zmass > 70 && C_Zmass < 110",
        "DR": "C_Dr_Z < 1.5",
        "Mass": "C_Sum_mass > 120",
        "Pt": "C_Sum_pt > 110",
        "vars": {
            "Zmass": ("C_Zmass", "C Z Mass [GeV]", 20, 50, 150),
            "DR": ("C_Dr_Z", "C dR", 25, 0, 2),
            "Sum_mass": ("C_Sum_mass", "C Sum Mass [GeV]", 25, 0, 700),
            "Sum_pt": ("C_Sum_pt", "C Sum pT [GeV]", 25, 0, 800)
        }
    },
    "D": {
        "base": "nev_passD == 1 && nLeptons >= 3",
        "Z": "D_Zmass > 70 && D_Zmass < 110",
        "DR": "D_Dr_Z < 1.5",
        "Mass": "D_Sum_mass > 120",
        "Pt": "D_Sum_pt > 110",
        "vars": {
            "Zmass": ("D_Zmass", "D Z Mass [GeV]", 20, 50, 150),
            "DR": ("D_Dr_Z", "D dR", 25, 0, 2),
            "Sum_mass": ("D_Sum_mass", "D Sum Mass [GeV]", 25, 0, 700),
            "Sum_pt": ("D_Sum_pt", "D Sum pT [GeV]", 25, 0, 800)
        }
    },
    # Puedes duplicar este bloque para canales B, C, y D si lo deseas
}


def plot_signal_vs_background(df_sig, df_bkg, cut, var, title, name, bins, xmin, xmax):
    df_sig_cut = df_sig.Filter(cut)
    df_bkg_cut = df_bkg.Filter(cut)

    model_sig = r.ROOT.RDF.TH1DModel(f"h_{name}_sig", "", bins, xmin, xmax)
    model_bkg = r.ROOT.RDF.TH1DModel(f"h_{name}_bkg", "", bins, xmin, xmax)

    h_sig = df_sig_cut.Histo1D(model_sig, var).GetValue()
    h_bkg = df_bkg_cut.Histo1D(model_bkg, var).GetValue()

    if h_sig.Integral() > 0:
        h_sig.Scale(1.0 / h_sig.Integral())
    if h_bkg.Integral() > 0:
        h_bkg.Scale(1.0 / h_bkg.Integral())

    # Estilo
    h_sig.SetFillColor(r.kRed + 1)
    h_sig.SetFillStyle(3004)  # Rayado diagonal
    h_sig.SetLineColor(r.kBlack)
    h_sig.SetMarkerStyle(20)

    h_bkg.SetFillColor(r.kBlue + 1)
    h_bkg.SetLineColor(r.kBlack)
    h_bkg.SetMarkerStyle(20)
    
    sig_max = h_sig.GetMaximum()
    bkg_max = h_bkg.GetMaximum()
    ymax = max(bkg_max, sig_max) * 1.2

    # Canvas
    stack = r.THStack(f"stack_{name}", "")
    stack.Add(h_bkg)
    stack.SetMaximum(bkg_max*1.2)

    canv = CMS.cmsCanvas(name, xmin, xmax, 0.0, ymax, title, "Events", square=CMS.kSquare, extraSpace=0.05, iPos=0)
    leg = CMS.cmsLeg(0.7, 0.78, 0.95, 0.88, textSize=0.04)

    
    # Ajustar eje Y
    #h_bkg.GetYaxis().SetTitleOffset(1.6)
    
    #stack.SetMaximum(max_y * 1.2)
    
    CMS.cmsDrawStack(stack, leg, {"Background": h_bkg})
    h_sig.Draw("HIST same")
    h_sig.Draw("P same")
    CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleOffset(1.6)


    leg.AddEntry(h_sig, "Signal", "lep")
    leg.Draw()
    CMS.CMS_lumi(canv, 0, 0)

    canv.SaveAs(f"/eos/user/o/olfimbre/NanoAOD_Filtered/histogramas/{name}.pdf")
    canv.Close()

cutflow_txt_path = "/eos/user/o/olfimbre/NanoAOD_Filtered/histogramas/cutflow.txt"
cutflow_latex_path = "/eos/user/o/olfimbre/NanoAOD_Filtered/histogramas/cutflow_latex.txt"


with open(cutflow_txt_path, "w") as f_txt, open(cutflow_latex_path, "w") as f_latex:
    #Encabezado txt normal
    f_txt.write("\nCUTFLOW TABLE:\n")
    f_txt.write(f"{'Channel':<8}{'Step':<15}{'Events Signal':>15}{'Events Background':>20}\n")
    f_txt.write("-" * 60 + "\n")
    
    #Encabezado latex
    f_latex.write("\\begin{table}[ht!]\n\\centering\n")
    f_latex.write("\\begin{tabular}{|c|c|c|c|}\n")
    f_latex.write("\\hline\n")
    f_latex.write("Channel & Step & Events Signal & Events Background \\\\\n")
    f_latex.write("\\hline\n")


    # Loop para generar todos los histogramas
    for ch in CUTS:
        c = CUTS[ch]
        base = c["base"]
    
        cut_Z = f"{base} && {c['Z']}"
        cut_DR = f"{cut_Z} && {c['DR']}"
        cut_mass = f"{cut_DR} && {c['Mass']}"
        cut_pt = f"{cut_mass} && {c['Pt']}"

        cuts_map = {
            "Z1": cut_Z,
            "DR2": cut_DR,
            "Mass3": cut_mass,
            "Pt4": cut_pt
        }

        
        
        
        for stage, full_cut in cuts_map.items():
            sig_events = df_signal.Filter(full_cut).Count().GetValue()
            bkg_events = df_background.Filter(full_cut).Count().GetValue()
            
            #escribir en txt normal
            f_txt.write(f"{ch:<8}{stage:<15}{sig_events:>15}{bkg_events:>20}\n")
            
            #escribir en latex
            f_latex.write(f"{ch} & {stage} & {sig_events} & {bkg_events} \\\\\n")
            f_latex.write("\\hline\n")
        
            for var_key, (var, title, bins, xmin, xmax) in c["vars"].items():
                name = f"{ch}_{stage}_{var_key}"
                print(f"Plotting {name}...")
                plot_signal_vs_background(df_signal, df_background, full_cut, var, title, name, bins, xmin, xmax)
    
    f_latex.write("\\end{tabular}\n")
    f_latex.write("\\caption{Cutflow para cada canal y paso de seleccion.}\n")
    f_latex.write("\\label{tab:cutflow}\n")
    f_latex.write("\\end{table}\n")