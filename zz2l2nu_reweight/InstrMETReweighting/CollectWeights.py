import ROOT
import os
import yaml


def main():
    ROOT.gROOT.LoadMacro('InstrMETReweighting.C')

    with open('config.yaml', 'r') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)

    input_dir = yaml_data["input_dir"]
    output_dir = yaml_data["output_dir"]
    samples = yaml_data["samples"]

    original_weight_name = 'original_weight_in_photonCR'
    transfer_factor_name = 'transfer_factor_photonCR_to_DY'

    for sample in samples:
        print(sample)

        name = sample["name"]
        multiplier = sample["multiplier"]

        input_filepath = os.path.join(input_dir, f'{name}.root')
        # output_filepath = os.path.join(output_dir, f'{name}_nvtx_reweighted.root')
        # output_filepath = os.path.join(output_dir, f'{name}_nvtx_eta_reweighted.root')
        # output_filepath = os.path.join(output_dir, f'{name}_nvtx_eta_pt_reweighted.root')
        output_filepath = os.path.join(output_dir, f'DYJetsToLL_Data-driven_{name}_nvtx_eta_pt_reweighted.root')

        df = ROOT.RDataFrame("Vars", input_filepath)

        transfer_factor = f"Float_t ({multiplier}) * photon_nvtx_reweighting * photon_eta_reweighting * photon_pt_reweighting"

        df_new = ((df
                   if df.HasColumn("weight")
                   else df.Define("weight", "trigger_weight * beam_halo_weight")
                   )
                  .Define("photon_nvtx_reweighting", f"GetNvtxWeight(num_pv_good, jet_cat)")
                  .Define("photon_eta_reweighting", f"GetEtaWeight(photon_eta, jet_cat)")
                  .Define("photon_pt_reweighting", f"GetPtWeight(photon_pt, jet_cat)")
                  .Define(original_weight_name, "weight")
                  .Define(transfer_factor_name, transfer_factor)
                  .Redefine("weight", f"{transfer_factor_name} * weight")
                #   .Filter("ptmiss > 120.", "ptmiss_geq_120")
                  .Define("ll_pt", "photon_pt")
                  .Define("ll_eta", "photon_eta")
                  .Define("ll_phi", "photon_phi")
                  .Define("lepton_cat", "-1")
                  )

        df.Report().Print()
        df_new.Snapshot("Vars", output_filepath)


if __name__ == "__main__":
    main()
