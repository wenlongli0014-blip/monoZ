import argparse
import ROOT
import ctypes

def main():
    rdf = ROOT.RDataFrame("Vars", f"../batch_singleelectron_{year}/merged/DYJetsToLL_Data-driven_GJets_nvtx_eta_pt_reweighted.root")

    rdf_new = (rdf.Filter('lepton_cat != 2', 'ossf')
                .Filter('jet_cat == 2.', '2j  ')
                .Filter('dijet_mass > 400.', 'dijet_mass_gt400')
                .Filter('std::abs(jet_eta[0] - jet_eta[1]) > 2.5', 'dijet_deta_gt2p5')
                .Filter('ptmiss > 120.', 'ptmiss_gt120')
                )

    list_binEdges_ptmiss = [120., 160., 200., 1200.]

    binEdges_ptmiss = (ctypes.c_double * len(list_binEdges_ptmiss)) (*list_binEdges_ptmiss)
    nbins_ptmiss = len(list_binEdges_ptmiss) - 1

    list_binEdges_llpt = [60., 82.5,  150., 300., 800.]
    binEdges_llpt = (ctypes.c_double * len(list_binEdges_llpt)) (*list_binEdges_llpt)
    nbins_llpt = len(list_binEdges_llpt) - 1

    model = ROOT.RDF.TH2DModel("histo2D", f"G+Jets ({year}) in SR", nbins_ptmiss, binEdges_ptmiss, nbins_llpt, binEdges_llpt)
    hist = rdf_new.Histo2D(model, "ptmiss", "ll_pt", "weight")
    bin_160_200_GJets = hist.GetBinContent(2,1)
    bin_200_1200_GJets = hist.GetBinContent(3,1)

    input_file = ROOT.TFile(f"histo2d_{year}_SR.root")
    hist_DD = input_file.Get("histo2D")

    bin_120_160_DD = hist_DD.GetBinContent(1,1)

    list_binEdges_ptmiss = [120., 200., 1200.]
    binEdges_ptmiss = (ctypes.c_double * len(list_binEdges_ptmiss)) (*list_binEdges_ptmiss)
    nbins_ptmiss = len(list_binEdges_ptmiss) - 1

    list_binEdges_llpt = [60.,  150., 300., 800.]
    binEdges_llpt = (ctypes.c_double * len(list_binEdges_llpt)) (*list_binEdges_llpt)
    nbins_llpt = len(list_binEdges_llpt) - 1

    model = ROOT.RDF.TH2DModel("histo2D", f"G+Jets ({year}) in SR", nbins_ptmiss, binEdges_ptmiss, nbins_llpt, binEdges_llpt)
    hist = rdf_new.Histo2D(model, "ptmiss", "ll_pt", "weight")
    bin_160_200_GJets = hist.GetBinContent(2,1)
    bin_200_1200_GJets = hist.GetBinContent(3,1)

    input_file = ROOT.TFile(f"histo2d_{year}_SR_v2.root")
    hist_DD = input_file.Get("histo2D")

    bin_120_160_DD = hist_DD.GetBinContent(1,1)

    list_binEdges_ptmiss = [120., 200., 1200.]
    binEdges_ptmiss = (ctypes.c_double * len(list_binEdges_ptmiss)) (*list_binEdges_ptmiss)
    nbins_ptmiss = len(list_binEdges_ptmiss) - 1

    list_binEdges_llpt = [60.,  150., 300., 800.]
    binEdges_llpt = (ctypes.c_double * len(list_binEdges_llpt)) (*list_binEdges_llpt)
    nbins_llpt = len(list_binEdges_llpt) - 1

    final_hist = ROOT.TH2D("histo2D", f"Data-driven DY {year} estimate in SR", nbins_ptmiss, binEdges_ptmiss, nbins_llpt, binEdges_llpt)

    final_hist.SetBinContent(1,1,bin_120_160_DD+bin_160_200_GJets+hist_DD.GetBinContent(1,2)+hist_DD.GetBinContent(2,2))
    final_hist.SetBinContent(2,1,bin_200_1200_GJets+hist_DD.GetBinContent(3,2))
    final_hist.SetBinContent(1,2,hist_DD.GetBinContent(1,3)+hist_DD.GetBinContent(2,3))
    final_hist.SetBinContent(2,2,hist_DD.GetBinContent(3,3))
    final_hist.SetBinContent(1,3,hist_DD.GetBinContent(1,4)+hist_DD.GetBinContent(2,4))
    final_hist.SetBinContent(2,3,hist_DD.GetBinContent(3,4))

    output_file = ROOT.TFile(f"histo2d_{year}_SR_final.root", "RECREATE")

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    final_hist.GetYaxis().SetTitle("ll pt [GeV]")
    final_hist.GetXaxis().SetTitle("ptmiss [GeV]")
    final_hist.Draw("colz,error,text")

    ROOT.gPad.SaveAs(f"{year}_DataDriven_SR_final.pdf")
    final_hist.Write()
    output_file.Close()
    print("Updated file saved")



if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-y', '--year', choices={'2016HIPM', '2016noHIPM', '2016', '2016_combined_before_reweighting', '2017', '2018', '2017_and_2018'},
                            default='2018',
                            help='YEAR')
    args = argparser.parse_args()

    year = args.year

    main()

