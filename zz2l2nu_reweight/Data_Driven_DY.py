import argparse
import ROOT
import ctypes



def main():
    names = [
            'DYJetsToLL_Data-driven_Data_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_QCD_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_TGJets_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_TTGJets_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WGToLNuG_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_Data_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_DYJetsToLL_PtZ_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_GJets_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_QCD_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_ST_s-channel_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_ST_t-channel_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_ST_tW_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_TTTo2L2Nu_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_TTToSemiLeptonic_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_WWTo1L1Nu2Q_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_WWTo2L2Nu_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_WZTo1L1Nu2Q_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_WZTo1L3Nu_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_WZTo2Q2L_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_WZTo3LNu_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_ZZTo2L2Nu_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_ZZTo2Q2L_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_WJets_Data-driven_ZZTo4L_nvtx_eta_pt_reweighted.root',
            'DYJetsToLL_Data-driven_ZNuNuGJets_Semi-data-driven_nvtx_eta_pt_reweighted.root',
        ]

    paths = [f'../batch_singlephoton_{year}/merged/{name}' for name in names]

    rdf = ROOT.RDataFrame("Vars", paths)

    rdf_new = (rdf.Filter('lepton_cat != 2', 'ossf')
                .Filter('jet_cat == 2.', '2j  ')
                .Filter('dijet_mass > 400.', 'dijet_mass_gt400')
                .Filter('std::abs(jet_eta[0] - jet_eta[1]) > 2.5', 'dijet_deta_gt2p5')
                .Filter('ptmiss > 120.', 'ptmiss_gt120')
                # .Filter('ptmiss <= 120.', 'ptmiss_lt120')
                )
    rdf_new_DY = (rdf.Filter('lepton_cat != 2', 'ossf')
                .Filter('jet_cat == 2.', '2j  ')
                .Filter('ptmiss < 120.', 'ptmiss_gt120')
                )
    # rdf.Report().Print()

    # hist = rdf_new.Histo1D("ptmiss", "weight")
    # hist.Draw()
    # input()

    list_binEdges_ptmiss = [120.,160., 200., 1200.]
    list_binEdges_ptmiss_DY = [0.,20.,40.,60.,80.,100.,120.]

    binEdges_ptmiss = (ctypes.c_double * len(list_binEdges_ptmiss)) (*list_binEdges_ptmiss)
    nbins_ptmiss = len(list_binEdges_ptmiss) - 1

    binEdges_ptmiss_DY = (ctypes.c_double * len(list_binEdges_ptmiss_DY)) (*list_binEdges_ptmiss_DY)
    nbins_ptmiss_DY = len(list_binEdges_ptmiss_DY) - 1

    # list_binEdges_llpt = [60., 150., 300., 800.]
    list_binEdges_llpt = [60., 82.5, 150., 300., 800.]
    binEdges_llpt = (ctypes.c_double * len(list_binEdges_llpt)) (*list_binEdges_llpt)
    nbins_llpt = len(list_binEdges_llpt) - 1


    model = ROOT.RDF.TH2DModel("histo2D", f"Data-driven DY estimate ({year}) in SR", nbins_ptmiss, binEdges_ptmiss, nbins_llpt, binEdges_llpt)
    hist = rdf_new.Histo2D(model, "ptmiss", "ll_pt", "weight")

    model_DY = ROOT.RDF.TH2DModel("histo2D", f"Data-driven DY estimate ({year}) in lowMET", nbins_ptmiss_DY, binEdges_ptmiss_DY, nbins_llpt, binEdges_llpt)
    hist_DY = rdf_new_DY.Histo2D(model_DY, "ptmiss", "ll_pt", "weight")
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    canvas1 = ROOT.TCanvas("canvas1", "canvas1", 800, 800)
    canvas2 = ROOT.TCanvas("canvas2", "canvas2", 800, 800)

    canvas1.cd()
    hist.GetYaxis().SetTitle("ll pt [GeV]")
    hist.GetXaxis().SetTitle("ptmiss [GeV]")
    hist.Draw("colz,error,text")
    canvas1.SaveAs(f"{year}_ptllbinnedat80_SR.pdf")

    canvas2.cd()
    hist_DY.GetYaxis().SetTitle("ll pt [GeV]")
    hist_DY.GetXaxis().SetTitle("ptmiss [GeV]")
    hist_DY.Draw("colz,error,text")
    canvas2.SaveAs(f"{year}_ptllbinnedat80_lowMET.pdf")

    output_file = ROOT.TFile(f'histo2d_{year}_SR.root', 'recreate')
    output_file.cd()
    hist.Write()
    output_file.Close()

    output_file_DY = ROOT.TFile(f'histo2d_{year}_lowMET.root', 'recreate')
    output_file_DY.cd()
    hist_DY.Write()

    output_file_DY.Close()

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-y', '--year', choices={'2016HIPM', '2016noHIPM', '2016', '2016_combined_before_reweighting', '2017', '2018', '2017_and_2018'},
                            default='2018',
                            help='YEAR')
    args = argparser.parse_args()

    year = args.year

    main()