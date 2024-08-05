import ROOT
import os
import yaml
import argparse

def main():
    ROOT.gROOT.LoadMacro(f'WJets_{year}.C')

    with open(f'config_{year}.yaml', 'r') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)

    input_dir = yaml_data["input_dir"]
    output_dir = yaml_data["output_dir"]
    samples = yaml_data["samples"]

    original_weight_name = 'original_weight_in_eCR'
    transfer_factor_name = 'transfer_factor_eCR_to_photonCR'

    for sample in samples:
        print(sample)

        name = sample["name"]
        multiplier = sample["multiplier"]

        input_filepath = os.path.join(input_dir, f'{name}.root')
        output_filepath = os.path.join(output_dir, f'WJets_Data-driven_{name}.root')

        df = ROOT.RDataFrame("Vars", input_filepath)

        transfer_factor = f"Float_t ({multiplier}) * GetFakeRatio(electron_pt, electron_eta) * GetTrigEffRatio(electron_pt, electron_eta)"

        df_new = ((df
                   if df.HasColumn("weight")
                   else df.Define("weight", "1")
                   )
                  .Filter("abs(electron_eta) < 2.5", "photon_eta_cut")
                  .Filter("abs(electron_eta_sc) <= 1.4442 or jet_cat == 2",
                          "reject_endcap_photon_01j")
                  .Define(original_weight_name, "weight")
                  .Define(transfer_factor_name, transfer_factor)
                  .Redefine("weight", f"{transfer_factor_name} * weight")
                  .Define("photon_pt", "electron_pt")
                  .Define("photon_eta", "electron_eta")
                  .Define("photon_phi", "electron_phi")
                  )

        df.Report().Print()
        df_new.Snapshot("Vars", output_filepath)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-y', '--year', choices={'2016', '2016HIPM', '2016noHIPM', '2017', '2018'},
                            default='2018',
                            help='YEAR')
    args = argparser.parse_args()

    year = args.year

    main()