import math
import ROOT
from ROOT import TFile, TTree, TBranch
from array import array
import argparse
def main():
    print(year)

    for jet_cat in [0, 1, 2]:
        print(f"jet_cat:  {jet_cat}")

        input_filename_data = f"../batch_zgamma_{year}/merged/Data.root"
        rdf_data = ROOT.RDataFrame("Vars", input_filename_data)
        rdf_data_new = rdf_data.Define("pass_cut", f"photon_pt > 60 && jet_cat == {jet_cat}")
        count_data = rdf_data_new.Filter("pass_cut").Count().GetValue()
        unc_data = math.sqrt(count_data)
        # print(f"{input_filename_data} Number of Data events with photon_pt > 60 && jet_cat == {jet_cat}: {count_data}")
        # print(f"uncertainty: {unc_data}")

        input_filename_DY = f"../batch_zgamma_{year}/merged/DYJetsToLL_PtZ.root"
        rdf_DY = ROOT.RDataFrame("Vars", input_filename_DY)
        rdf_DY_new = rdf_DY.Define("pass_cut", f"photon_pt > 60 && jet_cat == {jet_cat}")
        weighted_count_DY = rdf_DY_new.Filter("pass_cut").Sum("weight").GetValue()
        unc_DY = math.sqrt(rdf_DY_new.Filter("pass_cut").Define("w2","weight*weight").Sum("w2").GetValue())
        # print(f"{input_filename_data} Number of DYJetsToLL_PtZ events with photon_pt > 60 && jet_cat == {jet_cat}: {weighted_count_DY}")
        # print(f"uncertainty: {unc_DY}")

        input_filename_TTG = f"../batch_zgamma_{year}/merged/TTGJets.root"
        rdf_TTG = ROOT.RDataFrame("Vars", input_filename_TTG)
        rdf_TTG_new = rdf_TTG.Define("pass_cut", f"photon_pt > 60 && jet_cat == {jet_cat}")
        weighted_count_TTG = rdf_TTG_new.Filter("pass_cut").Sum("weight").GetValue()
        unc_TTG = math.sqrt(rdf_TTG_new.Filter("pass_cut").Define("w2","weight*weight").Sum("w2").GetValue())
        # print(f"{input_filename_data} Number of TTGJets events with photon_pt > 60 && jet_cat == {jet_cat}: {weighted_count_TTG}")
        # print(f"uncertainty: {unc_TTG}")
        try:
            input_filename_ZG = f"../batch_zgamma_{year}/merged/ZLLGJets.root"
        except:
            input_filename_ZG = f"../batch_zgamma_{year}/merged/ZGTo2LG.root"
        rdf_ZG = ROOT.RDataFrame("Vars", input_filename_ZG)
        rdf_ZG_new = rdf_ZG.Define("pass_cut", f"photon_pt > 60 && jet_cat == {jet_cat}")
        weighted_count_ZG = rdf_ZG_new.Filter("pass_cut").Sum("weight").GetValue()
        unc_ZG = math.sqrt(rdf_ZG_new.Filter("pass_cut").Define("w2","weight*weight").Sum("w2").GetValue())
        # print(f"{input_filename_data} Number of ZGTo2LG events with photon_pt > 60 && jet_cat == {jet_cat}: {weighted_count_ZG}")
        # print(f"uncertainty: {unc_ZG}")
        ratio = count_data / weighted_count_ZG
        # ratio = (count_data - weighted_count_DY - weighted_count_TTG) / weighted_count_ZG
        print("ratio: ", ratio)
        print("stat unc.: ", ratio * math.sqrt( (unc_data ** 2 + unc_DY ** 2 + unc_TTG ** 2) / ((count_data - weighted_count_DY - weighted_count_TTG) ** 2) + (unc_ZG / weighted_count_ZG) ** 2) )
        if jet_cat == 2:
            file_in = TFile.Open(f'../batch_singlephoton_{year}/merged/ZNuNuGJets.root', 'UPDATE')
            tree_in = file_in.Get('Vars')

            tree_in.SetBranchStatus("weight", 0)
            file_out = TFile.Open(f'../batch_singlephoton_{year}/merged/ZNuNuGJets_Semi-data-driven.root', 'RECREATE')

            tree_out = tree_in.CloneTree(0)
            tree_in.SetBranchStatus("weight", 1)
            new_weight = array('f', [0.0])
            tree_out.Branch('weight', new_weight, 'weight/F')

            for ieve in range(tree_in.GetEntries()):
                tree_in.GetEntry(ieve)
                new_weight[0] = tree_in.weight * ratio 
                tree_out.Fill()

            tree_out.Write()
            file_out.Close()
            file_in.Close()

            print("New tree with updated 'weight' branch has been created and saved.")


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-y', '--year', choices={'2016', '2016HIPM', '2016noHIPM', '2017', '2018'},
                            default='2018',
                            help='YEAR')
    args = argparser.parse_args()

    year = args.year

    main()

