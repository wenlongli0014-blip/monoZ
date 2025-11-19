#!/usr/bin/env python

"""Plots distributions of data and simulation."""

import argparse
from array import array
import collections.abc
import copy
import itertools
import math
import numbers
import os

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

import ROOT

import yaml

from hzz import Hist1D, mpl_style

ROOT.PyConfig.IgnoreCommandLineOptions = True


class Selection:
    """Event selection as described in configuration.

    Properties:
        formula:  String with formula for event selection.
        weight:   String with formula for event weight in simulation.
        weight_data:  String with formula for event weight in data.
        tag:      String uniquely identifying this selection.
        label:    LaTeX label for this selection to be shown in plots.
        blind:    Boolean showing whether this region is blinded.
    """

    def __init__(self, config):
        """Initialize from configuration fragment."""
        self.formula = config['formula']
        self.weight = config['weight']
        self.weight_data = config.get('weight_data', '1')
        self.tag = config['tag']
        self.label = config.get('label', self.tag)
        self.blind = config.get('blind', False)

        if self.blind and self.weight_data != '1':
            raise RuntimeError(
                'Cannot use weighted data with a blind selection.')


class Variable:
    """Variable as described in configuration.

    Properties:
        formula:  String with formula to compute this variable.
        tag:      String uniquely identifying this variable.
        label:    LaTeX label for this variable to be used in plots.
        unit:     String representing unit of measurement.
        y_scale:  Scale for y axis of the panel with distributions.
    """

    def __init__(self, config):
        """Initialize from configuration fragment."""

        self.formula = config['formula']
        self.tag = config['tag']
        self.label = config.get('label', self.tag)
        self.unit = config.get('unit', '')
        self.y_scale = config.get('y_scale', 'log')
        self._default_x_scale = config.get('x_scale', 'linear')
        self._binning = config['binning']

    def binning(self, selection_tag):
        """Return binning for given selection tag.

        The binning is represented by a NumPy array.
        """

        binning_info = self._binning_info(selection_tag)
        if 'range' in binning_info:
            r = binning_info['range']
            n = binning_info['bins']
            if not isinstance(n, numbers.Integral):
                raise RuntimeError(
                    'Incorrect specification of binning: When "range" is '
                    'given, "bins" must be an integral number.')
            if binning_info.get('scale', self._default_x_scale) == 'log':
                return np.geomspace(r[0], r[1], n + 1)
            else:
                return np.linspace(r[0], r[1], n + 1)
        else:
            binning = binning_info['bins']
            if not isinstance(binning, collections.abc.Sequence):
                raise RuntimeError(
                    'Incorrect specification of binning: When "range" is '
                    'omitted, "bins" must be a sequence.')
            return np.asarray(binning)

    def x_scale(self, selection_tag):
        """Return x axis scale for given selection tag.

        Possible values are "linear" and "log".
        """

        binning_info = self._binning_info(selection_tag)
        return binning_info.get('scale', self._default_x_scale)

    def _binning_info(self, selection_tag):
        """Find binning configuration for given selection tag."""

        if selection_tag in self._binning:
            binning_info = self._binning[selection_tag]
        elif 'default' in self._binning:
            binning_info = self._binning['default']
        else:
            raise RuntimeError(
                f'Variable "{self.tag}" has no binning definition for '
                f'selection "{selection_tag}" nor a default one.')
        return binning_info


class Sample:
    """Sample as described in configuration.

    Properties:
        tree_name:  Name of tree to read from input ROOT files.
        files:  Paths to ROOT files included in this sample.
    """

    def __init__(self, config, path_prefix='', tree_name='Vars'):
        """Initialize from configuration fragment.

        To construct paths to files included in this sample, for each
        filename given in the configuration fragment prepend the
        provided prefix and each of possible endings.
        """

        self.tree_name = tree_name
        self.files = []
        for name in config['files']:
            file_found = False
            for ending in ['.root', '_weights.root']:
                try_path = f'{path_prefix}{name}{ending}'
                if os.path.isfile(try_path):
                    self.files.append(try_path)
                    file_found = True
                    break
            if not file_found:
                raise RuntimeError(
                    f'Failed to find file corresponding to name "{name}".')


class DecoratedSample(Sample):
    """Sample with decoration, as described in configuration.

    Properties (in addition to those in the base class):
        label:  LaTeX label for this sample to be used in plots.
        color:  Colour for this sample.
    """
    def __init__(self, config, path_prefix='', tree_name='Vars'):
        """Initialize from configuration fragment."""

        super().__init__(config, path_prefix, tree_name)
        self.label = config['label']
        self.color = config['color']
        self.tag = config['tag'] #only fixed in di-lepton tree configs


class Configuration:
    """Represents parsed top-level configuration.

    Properties:
        selections:   List of all defined selections.
        variables:    List of all defined variables.
        data_sample:  Sample defining real data.
        sim_samples:  List of DecoratedSample's defining simulated
                      processes.
        entries_label:  Label for the type of entries in the source
                        trees.
    """

    def __init__(self, path, sample_path_prefix=''):
        with open(path) as f:
            config = yaml.safe_load(f)

        self.selections = [Selection(cfg) for cfg in config['selections']]
        self.variables = [Variable(cfg) for cfg in config['variables']]
        if not sample_path_prefix:
            sample_path_prefix = config['samples'].get('path_prefix', '')
        tree_name = config['samples'].get('tree_name', 'Vars')
        self.data_sample = Sample(
            config['samples']['data'], sample_path_prefix, tree_name)
        self.sim_samples = [
            DecoratedSample(cfg, sample_path_prefix, tree_name)
            for cfg in config['samples']['simulation']
        ]
        self.entries_label = config['samples'].get('entries_label', 'Events')


class HistogramBuilder:
    """Fills histograms and provides access to them."""

    def __init__(self, config):
        """Initialize from an instance of Configuration."""

        self._config = config
        self._data_hists = None
        self._sim_hists = None

    def build(self):
        """Fill histograms.

        Process all combinations of selections and variables.
        """

        self._data_hists = self._process_sample(
            self._config.data_sample, is_sim=False)
        sim_hists = [
            self._process_sample(sample)
            for sample in self._config.sim_samples
        ]

        # Convert sim_hists from a list of dicts to a dict of lists
        self._sim_hists = {}
        for selection, variable in itertools.product(
            self._config.selections, self._config.variables
        ):
            self._sim_hists[selection.tag, variable.tag] = [
                h[selection.tag, variable.tag]
                for h in sim_hists
            ]

    def data_hist(self, selection_tag, variable_tag):
        """Return data histogram for given selection and variable.

        Represented with class Hist1D.
        """

        return self._data_hists[selection_tag, variable_tag]

    def sim_hists(self, selection_tag, variable_tag):
        """Return simulation histograms.

        Return a list of histograms for given selection and variable.
        The list is sorted in the same way as simulated samples in the
        configuration.  The histograms are represented with class
        Hist1D.
        """

        return self._sim_hists[selection_tag, variable_tag]

    def _process_sample(self, sample, is_sim=True):
        """Fill histograms for given sample.

        Arguments:
            sample:  An instance of Sample.
            is_sim:  Flag indicating whether this sample is data or
                     simulation.

        Return value:
            Dictionary that maps (selection.tag, variable.tag) to
            Hist1D or None.  The latter is used for data if the
            selection is blind.
        """

        histograms = {}
        chain = ROOT.TChain(sample.tree_name)
        for path in sample.files:
            chain.AddFile(path)
        data_frame = ROOT.RDataFrame(chain)

        proxies = {}
        for selection in self._config.selections:
            if not is_sim and selection.blind:
                for variable in self._config.variables:
                    histograms[selection.tag, variable.tag] = None
                continue

            df_filtered = data_frame.Filter(selection.formula).Define(
                '_weight',
                selection.weight if is_sim else selection.weight_data)
            for variable in self._config.variables:
                binning = array('d', variable.binning(selection.tag))
                proxies[selection.tag, variable.tag] = df_filtered.Define(
                    '_var', variable.formula
                ).Histo1D(
                    ROOT.RDF.TH1DModel('', '', len(binning) - 1, binning),
                    '_var', '_weight')

        for key, proxy in proxies.items():
            hist = Hist1D(proxy.GetValue())
            self._tidy_hist(hist)
            histograms[key[0], key[1]] = hist
        return histograms

    @staticmethod
    def _tidy_hist(hist):
        """Adjust the histogram in place.

        Include under- and overflow bins into neighbouring bins.  Clip
        negative bin contents.
        """

        hist.contents[1] += hist.contents[0]
        hist.contents[-2] += hist.contents[-1]
        hist.contents[0] = hist.contents[-1] = 0.

        hist.errors[1] = math.hypot(hist.errors[0], hist.errors[1])
        hist.errors[-2] = math.hypot(hist.errors[-2], hist.errors[-1])
        hist.errors[0] = hist.errors[-1] = 0.

        np.clip(hist.contents, 0., None, out=hist.contents)


def plot_data_sim(variable, data_hist, sim_hists_infos, selection,
                  save_path, formats=['pdf'],
                  entries_label='Events', info_label='', root_path=''):
    """Plot and compare distributions in data and simulation.

    Arguments:
        variable:    Description of the variable to be plotted.
        data_hist:   Hist1D representing distribution of data.  Ignored
                     if the selection is blind.
        sim_hists_infos:  List of pairs (Hist1D, DecoratedSample)
                          representing distributions of different
                          processes in simulation.
        selection:   Destription of the event selection used.
        save_path:   Path under which to save the figure, without the
                     file extension.
        formats:     Formats in which to save the figure.
        entries_label:  Label denoting type of entries in histograms.
        info_label:  LaTeX description to be added in the plot.

    Return value:
        None.
    """

    fig = plt.figure(figsize=(7.2, 4.8))
    fig.patch.set_alpha(0.)
    gs = mpl.gridspec.GridSpec(
        3, 2, hspace=0., wspace=0.,
        height_ratios=[2, 1, 1], width_ratios=[5.5, 1]
    )
    axes_distributions = fig.add_subplot(gs[0, 0])
    axes_composition = fig.add_subplot(gs[1, 0])
    axes_residuals = fig.add_subplot(gs[2, 0])

    sim_hist_total = Hist1D()
    for hist, _ in sim_hists_infos:
        sim_hist_total += hist

    if selection.blind:
        data_hist = copy.deepcopy(sim_hist_total)
        data_hist.errors = np.sqrt(data_hist.contents)

    binning = data_hist.binning
    widths = binning[1:] - binning[:-1]
    bin_centres = (binning[:-1] + binning[1:]) / 2 

    if root_path is not None:

        rf = ROOT.TFile(root_path, "RECREATE")
        if data_hist is not None:
            th1 = ROOT.TH1D("data", "data",
                            len(data_hist.binning)-1,
                            array('d', data_hist.binning))
            for i in range(len(data_hist.contents)):
                th1.SetBinContent(i, data_hist.contents[i])
                th1.SetBinError(i, data_hist.errors[i])
            th1.Write()

        for hist, sample in sim_hists_infos:

            th1 = ROOT.TH1D(sample.tag, sample.tag,
                            len(hist.binning)-1,
                            array('d', hist.binning))

            for i in range(len(hist.contents)):
                
                val = hist.contents[i]
                err = hist.errors[i]
                
                th1.SetBinContent(i, val)
                th1.SetBinError(i, err)

            th1.Write()

        rf.Close()

    # Distributions of data and total simulation
    handle_sim = axes_distributions.hist(
        binning[:-1], weights=sim_hist_total.contents[1:-1] / widths,
        bins=binning, histtype='stepfilled',
        ec='dodgerblue', fc=mpl.colors.to_rgba('dodgerblue', 0.2)
    )[2][0]
    handle_data = axes_distributions.errorbar(
        bin_centres, data_hist.contents[1:-1] / widths,
        yerr=data_hist.errors[1:-1] / widths,
        marker='o', color='black', ls='none'
    )
    axes_distributions.set_yscale(variable.y_scale)

    # The total expectation in not guaranteed to be positive in each
    # bin.  In the following will only plot bins where it is.
    good_bins = np.arange(len(binning) - 1)[sim_hist_total.contents[1:-1] > 0.]
    sim_total_filtered = sim_hist_total.contents[good_bins + 1]

    # Composition of simulation
    axes_composition.hist(
        [binning[good_bins]] * len(sim_hists_infos), bins=binning,
        weights=[
            hist.contents[good_bins + 1] / sim_total_filtered
            for hist, _ in reversed(sim_hists_infos)
        ],
        color=[info.color for _, info in reversed(sim_hists_infos)],
        stacked=True
    )
    axes_composition.set_ylim(0., 1.)
    axes_composition.yaxis.set_major_locator(mpl.ticker.NullLocator())
    axes_composition.yaxis.set_minor_locator(mpl.ticker.NullLocator())

    # Residuals
    sim_errors = np.zeros(len(binning))
    sim_errors[good_bins] = sim_hist_total.errors[good_bins + 1] \
        / sim_total_filtered
    sim_errors[-1] = sim_errors[-2]  # Needed for fill_between
    axes_residuals.fill_between(
        binning, 1. + sim_errors, 1. - sim_errors,
        step='post', lw=0., color='0.75'
    )
    axes_residuals.errorbar(
        bin_centres[good_bins],
        data_hist.contents[good_bins + 1] / sim_total_filtered,
        yerr=data_hist.errors[good_bins + 1] / sim_total_filtered,
        marker='o', color='black', ls='none'
    )
    axes_residuals.axhline(1., ls='dashed', lw=0.8, color='black')
    axes_residuals.set_ylim(0., 2.)

    for axes in [axes_distributions, axes_composition, axes_residuals]:
        axes.set_xlim(binning[0], binning[-1])

    if variable.x_scale(selection.tag) == 'log':
        for axes in [axes_distributions, axes_composition, axes_residuals]:
            axes.set_xscale('log')

        # Provide a formatter for minor ticks so that they get labelled.
        # Also change the formatter for major ticks in order to get a
        # consistent formatting (1000 intead of 10^3).
        axes_residuals.xaxis.set_major_formatter(mpl.ticker.LogFormatter())
        axes_residuals.xaxis.set_minor_formatter(
            mpl.ticker.LogFormatter(minor_thresholds=(2, 0.4)))

    for axes in [axes_distributions, axes_composition]:
        axes.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
        axes.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    if variable.unit:
        xlabel = f'{variable.label} [{variable.unit}]'
    else:
        xlabel = variable.label
    ylabel_distributions = r'$\langle\mathrm{{{}}}\//\/' \
        r'\mathrm{{{}}}\rangle$'.format(
            entries_label, variable.unit if variable.unit else 'unit')

    axes_residuals.set_xlabel(xlabel)
    axes_distributions.set_ylabel(ylabel_distributions)
    axes_composition.set_ylabel('Compos.')
    axes_residuals.set_ylabel(r'$\mathrm{Data}\//\/\mathrm{Exp.}$')
    fig.align_ylabels()

    if info_label:
        axes_distributions.text(
            1., 1.002, info_label,
            ha='right', va='bottom', transform=axes_distributions.transAxes
        )

    # Manually construct a common legend for all panels
    handles = [handle_data, handle_sim]
    labels = ['Pseudodata' if selection.blind else 'Data', 'Exp.']
    for _, info in sim_hists_infos:
        handles.append(mpl.patches.Patch(color=info.color))
        labels.append(info.label)
    axes_distributions.legend(
        handles, labels,
        bbox_to_anchor=(1., 1.), loc='upper left', frameon=False
    )

    for fmt in formats:
        fig.savefig(f'{save_path}.{fmt}')
    plt.close(fig)


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument('config', help='Path to configuration file.')
    arg_parser.add_argument(
        '-p', '--prefix', default='', help='Prefix for ROOT files.')
    arg_parser.add_argument(
        '-o', '--output', default='fig',
        help='Directory for produced figures.')
    arg_parser.add_argument('-f', '--formats', nargs='+', default=['pdf'],
                            help='Formats for plots.')
    arg_parser.add_argument(
        '-y', '--year', default='2016',
        help='Year of data taking.')
    args = arg_parser.parse_args()

    config = Configuration(args.config, args.prefix)
    plt.style.use(mpl_style)
    for selection in config.selections:
        try:
            os.makedirs(os.path.join(args.output, selection.tag))
        except FileExistsError:
            pass

    histograms = HistogramBuilder(config)
    histograms.build()

    for variable, selection in itertools.product(
        config.variables, config.selections
    ):
        plot_data_sim(
            variable,
            histograms.data_hist(selection.tag, variable.tag),
            list(zip(histograms.sim_hists(selection.tag, variable.tag),
                     config.sim_samples)),
            selection,
            os.path.join(args.output, selection.tag, variable.tag),
            formats=args.formats,
            entries_label=config.entries_label,
            info_label=', '.join(
                token for token in [selection.label, args.year] if token
            ),
            root_path=os.path.join(args.output, selection.tag, variable.tag)+".root",
        )
