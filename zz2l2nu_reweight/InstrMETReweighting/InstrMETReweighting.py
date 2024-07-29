import argparse
import ROOT
import os
import yaml


def main(step):

    print(step)

    ROOT.gROOT.LoadMacro('InstrMETReweighting.C')

    with open('config.yaml', 'r') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)

    input_dir = yaml_data["input_dir"]
    output_dir = yaml_data["output_dir"]
    samples = yaml_data["samples"]

    # original_weight_name = 'original_weight_in_photonCR'
    # transfer_factor_name = 'transfer_factor_photonCR_to_DY'
    # transfer_factor_name = 'nvtx_weight'

    for sample in samples:
        print(sample)

        name = sample["name"]
        # multiplier = sample["multiplier"]

        input_filepath = os.path.join(input_dir, f'{name}.root')
        output_filepath = os.path.join(output_dir, f'{name}.root')

        df = ROOT.RDataFrame("Vars", input_filepath)

        df_new = ((df
                    .Define("photon_nvtx_reweighting", f"GetNvtxWeight(num_pv_good, jet_cat)")
                    ) if step == 'nvtx'
                  else (df
                        .Define("photon_nvtx_reweighting", f"GetNvtxWeight(num_pv_good, jet_cat)")
                        .Define("photon_eta_reweighting", f"GetEtaWeight(photon_eta, jet_cat)")
                        ) if step == 'eta'
                  else (df
                       .Define("photon_nvtx_reweighting", f"GetNvtxWeight(num_pv_good, jet_cat)")
                       .Define("photon_eta_reweighting", f"GetEtaWeight(photon_eta, jet_cat)")
                       .Define("photon_pt_reweighting", f"GetPtWeight(photon_pt, jet_cat)")
                      ) if step == 'pt'
                  else None)

        df_new = df_new.Define("blah", "1")

        df.Report().Print()
        df_new.Snapshot("Vars", output_filepath)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(__doc__)
    arg_parser.add_argument('-s', '--step', choices={'nvtx', 'eta', 'pt'},
                            # default='nvtx',
                            # action='append',
                            help='Reweighting step to use. "nvtx" must be run'
                            'before "eta" and "pt" (it is necessary to re-run'
                            'the photon trees with other weights already'
                            'applied).')
    args = arg_parser.parse_args()

    main(args.step)
