
import ROOT
from ROOT import TFile, TH1F, TGraphErrors
import numpy as np

def main():
    print(year)
    # 打开输入 ROOT 文件
    file_in = TFile.Open(f'../batch_egammafrommisid_{year}/merged/Data.root', 'READ')
    tree_in = file_in.Get('Vars')

    # 创建一个新的 ROOT 文件用于保存输出直方图
    output_file = TFile(f'ratiosWJets/fakeRatio/fake_ratio_{year}_hist.root', 'RECREATE')

    # 定义自定义分bin边界
    pt_bins = np.array([50, 75, 100, 135, 220, 450], dtype='float64')

    # 创建直方图用于存储 event_cat == 1 和 event_cat == 0 的计数
    hist_cat1_endcap = TH1F('hist_cat1_endcap', 'Event Cat 1 Count for Endcap', len(pt_bins)-1, pt_bins)
    hist_cat0_endcap = TH1F('hist_cat0_endcap', 'Event Cat 0 Count for Endcap', len(pt_bins)-1, pt_bins)
    hist_ratio_endcap = TH1F('ratio_endcap', 'Fake Ratio for Endcap', len(pt_bins)-1, pt_bins)

    hist_cat1_barrel = TH1F('hist_cat1_barrel', 'Event Cat 1 Count for Barrel', len(pt_bins)-1, pt_bins)
    hist_cat0_barrel = TH1F('hist_cat0_barrel', 'Event Cat 0 Count for Barrel', len(pt_bins)-1, pt_bins)
    hist_ratio_barrel = TH1F('ratio_barrel', 'Fake Ratio for Barrel', len(pt_bins)-1, pt_bins)

    # 遍历树并填充直方图
    for i in range(tree_in.GetEntries()):
        tree_in.GetEntry(i)
        if abs(tree_in.probe_eta) >= 1.5660 and 81.1876 <= tree_in.tot_mass <= 101.1876: # endcap
            if tree_in.event_cat == 1:
                hist_cat1_endcap.Fill(tree_in.probe_pt)
            elif tree_in.event_cat == 0:
                hist_cat0_endcap.Fill(tree_in.probe_pt)
        elif abs(tree_in.probe_eta) <= 1.4442 and 81.1876 <= tree_in.tot_mass <= 101.1876: # barrel
            if tree_in.event_cat == 1:
                hist_cat1_barrel.Fill(tree_in.probe_pt)
            elif tree_in.event_cat == 0:
                hist_cat0_barrel.Fill(tree_in.probe_pt)
    for bin_idx in range(1, hist_ratio_endcap.GetNbinsX() + 1):
        count_cat1 = hist_cat1_endcap.GetBinContent(bin_idx)
        count_cat0 = hist_cat0_endcap.GetBinContent(bin_idx)
        if count_cat0 != 0:
            ratio = count_cat1 / count_cat0
        else:
            ratio = 0
        hist_ratio_endcap.SetBinContent(bin_idx, ratio)

    for bin_idx in range(1, hist_ratio_barrel.GetNbinsX() + 1):
        count_cat1 = hist_cat1_barrel.GetBinContent(bin_idx)
        count_cat0 = hist_cat0_barrel.GetBinContent(bin_idx)
        if count_cat0 != 0:
            ratio = count_cat1 / count_cat0
        else:
            ratio = 0
        hist_ratio_barrel.SetBinContent(bin_idx, ratio)

    # 写入直方图到输出文件
    output_file.cd()
    hist_ratio_endcap.Write()
    hist_ratio_barrel.Write()

    # 关闭文件
    output_file.Close()
    file_in.Close()

    print(f"{year} ROOT file 'fake_ratio_file.root' with histograms 'ratio_endcap' and 'ratio_barrel' created.")

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-y', '--year', choices={'2016', '2016HIPM', '2016noHIPM', '2017', '2018'},
                            default='2018',
                            help='YEAR')
    args = argparser.parse_args()

    year = args.year

    main()