import os
import re

import yaml


class Dataset:
    """Describes a dataset.

    An object of this class contains the following fields:
        path:    Path to the dataset definition file.
        name:    String identifying the dataset.
        files:   Paths to input ROOT files included in the dataset.
        is_sim:  Boolean indicating whether the dataset is real data or
            simulation.
        parameters:  Mapping with all parameters extracted from the
            definition file.  Integer and floating-point values are
            converted into the corresponding native representations.

    This class offers functionality similar although not identical to
    class DatasetInfo in the C++ part.  It parses YAML dataset
    definition files [1].

    [1] https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/wikis/dataset-definitions
    """

    def __init__(self, path, stems={}):
        """Initialize from a dataset definition file.

        Arguments:
            path:   Path to a dataset definition file.  Fragments are
                supported.
            stems:  Stems to be used to construct full dataset
                definitions from fragments.  Must be represented by a
                mapping with the names of the stems as the keys.
        """

        self.path = path
        self.stems = stems
        self.name = ''
        self.parameters = {}
        self.is_sim = None
        self.files = []

        self._from_yaml(path)

    def save(self, path):
        """Save to a file.

        Write the full definition file even if the dataset has been
        constructed from a definition fragment.
        """

        self._save_yaml(path)

    def _from_yaml(self, path):
        """Initialize from YAML format."""

        with open(path, 'r') as f:
            loaded_dict = yaml.safe_load(f)

        if 'stem' in loaded_dict:
            self._splice_yaml(loaded_dict)

        for parameter_name in ['name', 'is_sim', 'files']:
            if parameter_name not in loaded_dict:
                raise RuntimeError(
                    'Mandatory parameter "{}" is missing in dataset '
                    'definition file "{}".'.format(parameter_name, path)
                )
            else:
                setattr(self, parameter_name, loaded_dict.pop(parameter_name))

        self.parameters = loaded_dict

    def _save_yaml(self, path):
        """Save in YAML format."""

        write_dict = {}
        write_dict.update(self.parameters)

        if 'name' not in write_dict:
            write_dict['name'] = self.name

        if 'is_sim' not in write_dict:
            write_dict['is_sim'] = self.is_sim

        write_dict['files'] = self.files

        with open(path, 'w') as f:
            yaml.dump(write_dict, f, default_flow_style=False)

    def _splice_yaml(self, config):
        """Incorporate the stem into the loaded configuration."""

        if not self.stems:
            raise RuntimeError(
                'Dataset definition file "{}" is a fragment. '
                'Cannot construct the full definition because no '
                'stems have been provided.'.format(self.path)
            )

        try:
            stem = self.stems[config['stem']]
        except KeyError:
            raise RuntimeError(
                'Stem "{}" required by datsaet definition fragment "{}" '
                'is not found.'.format(config['stem'], self.path)
            )

        del config['stem']
        config.update(stem)


def parse_datasets_file(path, config_path=''):
    """Parse file with a list of dataset definition files.

    Paths to dataset definition files (DDF) are given in a text file,
    one per line.  A common directory with respect to which these paths
    are resolved can be specified optionally.  Empty lines and lines
    that only contain comments are skipped.

    Arguments:
        path:  Path to a file with a list of dataset definition files.
        config_path:  Path to the master configuration file.  Resolved
            with respect to $HZZ2L2NU_BASE/config.

    Return value:
        List of constructed datasets.
    """

    # Extract paths to dataset definition files
    ddfs = []

    blank_regex = re.compile(r'^\s*$')
    comment_regex = re.compile(r'^\s*#')
    directory_regex = re.compile(r'^\s*catalogPath\s*=\s*(.+)\s*$')

    datasets_file = open(path, 'r')
    directory = ''

    for line in datasets_file:
        if blank_regex.match(line) or comment_regex.search(line):
            continue

        match = directory_regex.match(line)

        if match:
            directory = match.group(1)
        else:
            ddfs.append(os.path.join(directory, line.strip()))

    datasets_file.close()

    # Read stem dataset definitions if available
    stems = {}

    if config_path:
        config_dir = os.path.join(os.environ['HZZ2L2NU_BASE'], 'config/')

        with open(config_dir + config_path) as f:
            config = yaml.safe_load(f)

        if 'dataset_stems' in config:
            stems = {}
            for stem_filename in config['dataset_stems']:
                with open(config_dir + stem_filename) as f:
                    for stem in yaml.safe_load(f):
                        stems[stem['name']] = stem

    datasets = [Dataset(ddf, stems) for ddf in ddfs]
    return datasets
