#!/usr/bin/env python

"""Prepares HTCondor jobs for the analysis."""


import argparse
import math
import os

from hzz import SystDatasetSelector, parse_datasets_file


class JobBuilder:
    def __init__(
        self, task_dir, config_path, prog, prog_args=[], output_prefix='',
        events_per_job=500000
    ):
        """Initialize the builder.

        Arguments:
            task_dir:  Directory for the task.  Created if needed.
            config_path:  Path to master configuration for the analysis.
                Forwarded to the program to be run.
            prog:  Program to be run in a job.
            prog_args:  List with arguments to be forwared to the
                program to be run.
            output_prefix:  Prefix to be added to names of output files.
        """

        # Make sure task_dir is an absolute path because it is used to
        # define the output path for jobs
        self.task_dir = os.path.abspath(task_dir)

        self.config_path = config_path
        self.prog = prog
        self.prog_args = prog_args
        self.output_prefix = output_prefix
        self.install_path = os.environ['HZZ2L2NU_BASE']

        self._create_directories()
        # self.submit_commands = []
        self.job_names = []
        self.events_per_job = events_per_job

    def prepare_jobs(self, dataset, syst):
        """Construct jobs for given dataset.

        Decide on the job splitting and delegate the construction of
        individual jobs to the dedicated function.

        Arguments:
            dataset:  Dataset for which to construct jobs.
            syst:     Label of requested systematic variation.

        Return value:
            None.
        """

        print(
            'Preparing scripts for dataset '
            '\033[1;33m{}\033[0;m'.format(dataset.name)
        )

        # For some systematic variations all files within a dataset have to
        # be processed within a single job
        if syst and ('pdf' in syst or 'QCDscale' in syst):
            self._prepare_job_script(dataset, syst)
        else:
            # Choose the number of files per job based on the average number
            # of events per file, if relevant information is available
            num_events = dataset.parameters.get('num_selected_events', None)

            if num_events is not None:
                events_per_job = self.events_per_job

                if num_events > events_per_job:
                    job_splitting = int(math.ceil(
                        len(dataset.files) / (num_events / events_per_job)
                    ))
                else:
                    job_splitting = len(dataset.files)
            else:
                job_splitting = 25

            num_jobs = int(math.ceil(len(dataset.files) / job_splitting))

            for job_id in range(num_jobs):
                self._prepare_job_script(
                    dataset, syst, job_id,
                    skip_files=job_id * job_splitting,
                    max_files=job_splitting
                )

    def write_submit_script(self, script_path):
        """Write script to submit jobs.

        All jobs created so far by method prepare_jobs are included.
        """

        job_names_file_path = '{}/job_names.dat'.format(
            self.task_dir
        )

        with open(job_names_file_path, 'w') as f:
            for line in self.job_names:
                f.write(line)
                f.write('\n')

        submit_file_path = '{}/{}jobs.sub'.format(
            self.task_dir, self.output_prefix
        )

        submit_file_content = [
            'executable              = {}/jobs/scripts/runOnBatch_{}$(Job_name).sh'.format(self.task_dir, self.output_prefix),
            'log                     = {}/jobs/logs/runOnBatch_{}$(Job_name).$(Cluster).log'.format(self.task_dir, self.output_prefix),
            'output                  = {}/jobs/logs/runOnBatch_{}$(Job_name).outfile.$(Cluster).txt'.format(self.task_dir, self.output_prefix),
            'error                   = {}/jobs/logs/runOnBatch_{}$(Job_name).errors.$(Cluster).txt'.format(self.task_dir, self.output_prefix),
            'should_transfer_files   = Yes',
            'when_to_transfer_output = ON_EXIT',
            'queue Job_name from {}'.format(job_names_file_path)
        ]

        with open(submit_file_path, 'w') as f:
            for line in submit_file_content:
                f.write(line)
                f.write('\n')

        with open(script_path, 'w') as f:
            f.write('condor_submit {}'.format(submit_file_path))
            f.write('\n')

    def _create_directories(self):
        """Create directory structure for the task if needed."""

        if not os.path.exists(self.task_dir):
            print('\033[1;31mWill create task directory '
                  '"{}"\033[0;m'.format(self.task_dir))
            os.makedirs(self.task_dir)

        sub_dirs = ['output', 'jobs', 'jobs/scripts', 'jobs/logs']

        for sub_dir in sub_dirs:
            try:
                os.makedirs(os.path.join(self.task_dir, sub_dir))
            except OSError:
                pass

    def _prepare_job_script(
        self, dataset, syst, job_id=0, skip_files=0, max_files=-1
    ):
        """Create script defining the job.

        The script performs set-up and executes specified program.  It
        is saved in subdirectory jobs/scripts.

        Arguments:
            dataset:  Dataset to be processed.
            syst:     Label of requested systematic variation.
            job_id:   Number for the current job among all jobs in the given
                dataset.
            skip_files:  Number of input files from the dataset to skip.
            max_files:   Maximal number of input files to process in this
                job.  A value of -1 means all remaining files.

        Return value:
            None.
        """

        if syst:
            job_name = '{}_{}_{}'.format(dataset.name, syst, job_id)
        else:
            job_name = '{}_{}'.format(dataset.name, job_id)

        script_commands = [
          'export INITDIR={}'.format(self.install_path),
          'cd $INITDIR',
          '. ./env.sh',
          'cd -',
          'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi',
          'hostname',
          'date'
        ]

        # Construct options for the program
        options = [
            '--config={}'.format(self.config_path),
            '--ddf={}'.format(dataset.path),
            '--output={}{}.root'.format(self.output_prefix, job_name),
            '--skip-files={}'.format(skip_files),
            '--max-files={}'.format(max_files), '--max-events=-1'
        ]

        options += self.prog_args

        if syst:
            options.append('--syst={}'.format(syst))

        # Use debug-level verbosity when running on the batch system
        # since one is expected to consult logs only in case of a
        # problem
        options.extend(['-v', '2'])

        run_application_command = ' '.join([self.prog] + options)
        script_commands.append('echo ' + run_application_command)
        script_commands.append(run_application_command + ' || exit $?')

        script_commands.append('cp {}{}.root {}/output'.format(self.output_prefix, job_name,"/eos/user/h/hgao/batch_dilepton_2018")
        # script_commands.append('cp {}{}.root {}/output'.format(
        #     self.output_prefix, job_name, self.task_dir
        # ))

        script_path = '{}/jobs/scripts/runOnBatch_{}{}.sh'.format(
            self.task_dir, self.output_prefix, job_name
        )

        with open(script_path, 'w') as f:
            for command in script_commands:
                f.write(command)
                f.write('\n')

        self.job_names.append(
            job_name
        )


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument(
        'datasets', help='File with a list of datasets to process.'
    )
    arg_parser.add_argument(
        'prog_args', nargs='*',
        help='Arguments for program to run. Some standard arguments are added '
        'automatically.'
    )
    arg_parser.add_argument(
        '-d', '--task-dir', default='task',
        help='Directory for scripts and results of this task.'
    )
    arg_parser.add_argument(
        '--config', default='2016.yaml',
        help='Master configuration for the analysis.'
    )
    arg_parser.add_argument(
        '--prog', default='runHZZanalysis',
        help='Program to run. Full path must be given if it is not placed '
        'in a standard location.'
    )
    arg_parser.add_argument(
        '--syst', default='',
        help='Requested systematic variation or a group of them.'
    )
    arg_parser.add_argument(
        '--prefix', default='',
        help='Prefix for names of output files.'
    )
    arg_parser.add_argument(
        '--split-weights', action='store_true',
        help='Requests that weight-based systematic variations are processed '
        'in separate jobs.'
    )
    arg_parser.add_argument(
        '--events-perjob', default=500000, type=int,
        help='Maximum number of events per job.'
    )
    args = arg_parser.parse_args()

    if args.syst == 'no':
        args.syst = ''

    job_builder = JobBuilder(
        args.task_dir, args.config, args.prog, prog_args=args.prog_args,
        output_prefix=args.prefix, events_per_job=args.events_perjob
    )
    datasets = parse_datasets_file(args.datasets, args.config)

    # Central configuration for systematic variations
    if not args.syst or args.syst in {'all', 'weights'}:
        combine_weights = args.syst and not args.split_weights
        for dataset in datasets:
            job_builder.prepare_jobs(
                dataset,
                'weights' if dataset.is_sim and combine_weights else ''
            )

    if args.syst and args.syst != 'weights':
        # Systematic variations
        dataset_selector = SystDatasetSelector(
            os.path.join(job_builder.install_path, 'config/syst.yaml')
        )

        for variation, dataset in dataset_selector(
            datasets, args.syst, skip_nominal=True,
            combine_weights=not args.split_weights
        ):
            job_builder.prepare_jobs(dataset, variation)

    job_builder.write_submit_script(
        os.path.join(args.task_dir, 'send_jobs.sh')
    )
