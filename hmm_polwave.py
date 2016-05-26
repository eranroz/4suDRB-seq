"""
A script to define the polymerase wave frontier based on 2 state HMM, for 4SU-DRB-seq.

"""
import os
import requests

import numpy as np
from scipy.io import loadmat
# reuse the code for hmm and sequence publishing
import pandas as pd
from hmm_kit import HMMModel
import hmm_kit.bwiter as bwiter

from config import *


class DiscreteTransformer(object):
    """
    A class to transform sequences to discrete alphabet
    @param percentiles: percentiles according to which assign letter (0-100).
                        Values between 0 to value of percentiles[0] assign letter 0. etc

    @note For consistence transform the class use the values of percentiles of the first sequence,
          for latter called sequences
    """

    def __init__(self, percentiles=[60, 75, 90], percentiles_values=None, log_transform=True):
        self.percentiles = percentiles
        self.percentiles_values = percentiles_values
        self.log_transform = log_transform

    def __call__(self, *args, **kwargs):
        """
        Transforms the sequence to discrete alphabet
        """
        if self.log_transform:
            lg_data_n_zero = np.log1p(np.asfarray(args[0]))
        else:
            lg_data_n_zero = np.asfarray(args[0])

        if self.percentiles_values is None:
            self.percentiles_values = np.percentile(lg_data_n_zero, q=self.percentiles)

        mapped2 = np.zeros(len(lg_data_n_zero), dtype=int)  # * len(self.percentiles_values)
        for i, bounds in enumerate(self.percentiles_values):
            mapped2[lg_data_n_zero > bounds] = i + 1
        return mapped2


def fit_hmm_multi(cond=None, exp_time=None, rep_samples=None, suffix='', genome='hg19'):
    """
    Finds Pol II for multiple experiments. Invokes hmm_fit for each of them
    :param cond: condition such as control or knockdown
    :param exp_time: time of the experiments (usually in minutes) between DRB wash and harvesting
    :param rep_samples: replicates
    :param suffix: suffix to add to the output file
    :param genome: genome as appear in the UCSC
    :return:
    """
    from config import ALL_REPLICATES, KNOWN_GENES

    if rep_samples is None or rep_samples[0] is None:
        rep_samples = ALL_REPLICATES

    print('Loading genes')
    download_clustersWithNames(genome)  # validate we have the metadata from UCSC
    genes = pd.read_csv(KNOWN_GENES, sep='\t',
                        usecols=['#hg19.knownGene.name', 'hg19.knownGene.chrom', 'hg19.knownGene.strand',
                                 'hg19.knownGene.txStart', 'hg19.knownGene.txEnd', 'hg19.knownGene.exonStarts',
                                 'hg19.knownGene.exonEnds', 'hg19.kgXref.geneSymbol', 'hg19.knownIsoforms.clusterId'])

    genes.columns = ['name', 'chrom', 'strand', 'TSS', 'TES', 'exst', 'exen', 'symbol', 'clusterId']

    relevant_genes = np.abs(genes['TES'] - genes['TSS']) > MIN_GENE_LENGTH  # reduce to large enough genes
    print('%i transcripts will be analyzed, %i will be ignored (too short)' % (relevant_genes.sum(),
                                                                               (~relevant_genes).sum()))
    genes = genes[relevant_genes]

    diff_data = dict()

    print('Loading data for {}{}'.format(cond, exp_time))
    diff_data_time = dict()
    for replicate in rep_samples:
        try:
            startdata_path = TRANSCRIPTION_DATA_DIR + '{}_{}_{}.mat'.format(replicate, cond, '0')
            enddata_path = TRANSCRIPTION_DATA_DIR + '{}_{}_{}.mat'.format(replicate, cond, exp_time)
            print('Opening data for start time: %s' % startdata_path)
            start_data = loadmat(TRANSCRIPTION_DATA_DIR + '{}_{}_{}.mat'.format(replicate, cond, '0'))['chip']
            print('Opening data for end time: %s' % enddata_path)
            end_data = loadmat(enddata_path)['chip']
        except Exception as e:
            print('Couldn\'t load file in %s. Either missing or with wrong format' % TRANSCRIPTION_DATA_DIR)
            raise e
        for chrom in start_data.dtype.names:
            chrom_data = np.nan_to_num(end_data[chrom][0, 0][0]) - np.nan_to_num(start_data[chrom][0, 0][0])
            if chrom not in diff_data_time:
                diff_data_time[chrom] = []
            diff_data_time[chrom].append(
                np.maximum(chrom_data, 0))  # data in start>data in end => ineffective area

    key_name = '{}_{}'.format(exp_time, rep_samples[0])
    diff_data[key_name] = dict()
    for chrom, chrom_data in diff_data_time.items():
        diff_data[key_name][chrom] = np.mean(chrom_data, 0)

    for exp_name, cond_diff_data in diff_data.items():
        hmm_fit(cond_diff_data, genes, exp_name, suffix)


def constaint_training(new_mode):
    new_mode.state_transition[0, 1:] = [0.9, 0.1]


def hmm_fit(diff_data, genes, exp_name, suffix, one_model=True):
    """
    Discretization and hmm fit/predict for each transcript  in the genome
    :param diff_data: array of read coverage (data) after harvesting minus the read coverage in the begining
    :param genes: genes metadata from UCSC
    :param exp_name: name of the experiment
    :param suffix: suffix to append to the output file
    :param one_model: whether to train model for each transcript or use one model for all
    """
    percentiles_for_discretization = [5, 50, 75, 90]  # [50, 75, 90]
    discretor = DiscreteTransformer(percentiles_for_discretization)
    train_chrom = 'chr2'
    discrete_data_train = diff_data[train_chrom]

    one_model_train_data = []  # used for training
    # zero exons of chr2
    train_genes = genes[genes.loc[:, 'chrom'] == train_chrom]
    for clusterId, transcripts in train_genes.groupby('clusterId'):
        cluster_chrom = transcripts['chrom'].iloc[0]
        cluster_tss = transcripts['TSS'].min()
        cluster_tes = transcripts['TES'].max()
        cluster_strand = transcripts['strand'].iloc[0]
        indics = np.arange(np.ceil(cluster_tss / jump), np.floor(cluster_tes / jump), dtype=int)

        Icds = np.zeros(indics.shape, dtype=bool)

        cluster_exst = transcripts['exst'].apply(lambda x: x.split(','))
        cluster_exen = transcripts['exen'].apply(lambda x: x.split(','))
        for exst, exen in zip(cluster_exst, cluster_exen):
            exst = np.ceil(np.array(exst[:-1], dtype=int) / jump) - indics[0]
            exen = np.floor(np.array(exen[:-1], dtype=int) / jump) - indics[0]
            for ex, en in zip(exst, exen):
                Icds[ex:en] = True
        # zero exons data
        discrete_data_train[indics[Icds]] = 0

        gene_data = diff_data[cluster_chrom][indics]

        for tx_key, tx in transcripts.iterrows():
            tx_slice = slice(np.ceil(tx.TSS / jump) - np.ceil(cluster_tss / jump),
                             gene_data.shape[0] - (np.floor(cluster_tes / jump) - np.floor(tx.TES / jump)))

            tx_indics_rel = np.arange(gene_data.shape[0])[tx_slice]
            tx_indics_rel = tx_indics_rel[~(Icds[tx_indics_rel])]
            gene_data_no_icds = gene_data[tx_indics_rel]  # remove exons

            if gene_data_no_icds.shape[0] < min_sequence_length:
                print('Skipping - to small sequence')
                continue
            if cluster_strand == '-':
                gene_data_no_icds = gene_data_no_icds[::-1]

            if one_model:
                one_model_train_data.append(gene_data_no_icds)

    discretor(discrete_data_train[discrete_data_train > 0])  # fit the discretion to values higher than 0
    print('Discretization values')
    print(np.exp(discretor.percentiles_values) - 1)  # note: fixed to -1
    # train hmm for all genes

    if one_model:
        # discertize all of them!
        train_data = [discretor(seq) for seq in one_model_train_data]
        # ignore non-transcribed (at least 4*500=2000 bases)
        train_data = [seq for seq in train_data if np.sum(seq > 0) > 3]
        print('Distribution of discrete alphabet in training data (chr2, transcribed genes, no exons)')
        print(np.sum([np.bincount(seq, minlength=len(percentiles_for_discretization) + 1) for seq in train_data], 0))
        # multi train
        state_transition = np.array([
            [0.0, 0.9, 0.1],  # begin state moves to transcribed state
            [0.1, 0.999, 0.001],  # transcribed state - small P to get to non transcribed
            [0.9, 0, 1]  # 3' end
        ])

        emission = np.array([
            [0.1, 0.1, 0.2, 0.2, 0.2],  # begin state (dummy)
            [0.1, 0.1, 0.3, 0.4, 0.1],  # transcription
            [0.8, 0.1, 0.09, 0.009, 0.001]  # non transcription
        ])
        init_model = HMMModel.DiscreteHMM(state_transition, emission)
        print('Training HMM model using training set of %i transcripts. This may take few minutes' % len(train_data))
        model, prob = bwiter.bw_iter_multisequence(train_data, init_model, bwiter.DiffCondition(1),
                                                   constraint_func=constaint_training)
        print(model)

    # for each gene we train and fit seperatly, similar to Danko. et al.
    gene_hmm_data = []
    # group by cluster to avoid exons of other overlapping transcripts
    for clusterId, transcripts in genes.groupby('clusterId'):
        cluster_chrom = transcripts['chrom'].iloc[0]
        cluster_tss = transcripts['TSS'].min()
        cluster_tes = transcripts['TES'].max()
        cluster_strand = transcripts['strand'].iloc[0]

        indics = np.arange(np.ceil(cluster_tss / jump), np.floor(cluster_tes / jump), dtype=int)

        Icds = np.zeros(indics.shape, dtype=bool)

        cluster_exst = transcripts['exst'].apply(lambda x: x.split(','))
        cluster_exen = transcripts['exen'].apply(lambda x: x.split(','))
        for exst, exen in zip(cluster_exst, cluster_exen):
            exst = np.ceil(np.array(exst[:-1], dtype=int) / jump) - indics[0]
            exen = np.floor(np.array(exen[:-1], dtype=int) / jump) - indics[0]
            for ex, en in zip(exst, exen):
                Icds[ex:en] = True
        # should end with 31949
        gene_data = diff_data[cluster_chrom][indics]

        for tx_key, tx in transcripts.iterrows():
            tx_slice = slice(np.ceil(tx.TSS / jump) - np.ceil(cluster_tss / jump),
                             gene_data.shape[0] - (np.floor(cluster_tes / jump) - np.floor(tx.TES / jump)))
            tx_indics_rel = np.arange(gene_data.shape[0])[tx_slice]
            tx_indics_rel = tx_indics_rel[~(Icds[tx_indics_rel])]
            gene_data_no_icds = gene_data[tx_indics_rel]  # remove exons
            tx_indics = indics[tx_indics_rel]
            if gene_data_no_icds.shape[0] < min_sequence_length:
                print('Skipping - to small sequence')
                continue
            if cluster_strand == '-':
                gene_data_no_icds = gene_data_no_icds[::-1]
                tx_indics = tx_indics[::-1]

            tx_data = discretor(gene_data_no_icds)
            if np.sum(tx_data > 0) < 3:
                continue  # not expressed

            if not one_model:
                state_transition = np.array([
                    [0.0, 0.9, 0.1],  # begin state moves to transcribed state
                    [0.1, 0.8, 0.2],  # transcribed state - small P to get to non transcribed
                    [0.9, 0, 1]  # 3' end
                ])

                emission = np.array([
                    [0.25, 0.25, 0.25, 0.25],  # begin state
                    [0.1, 0.3, 0.4, 0.2],  # transcription
                    [0.9, 0.098, 0.001, 0.001]  # non transcription
                ])
                init_model = HMMModel.DiscreteHMM(state_transition, emission)
                model, p = bwiter.bw_iter(tx_data, init_model, 3, constraint_func=constaint_training)
                if p == np.nan:
                    raise Exception('bug')
            else:
                p = 1
            decoded_seq = model.viterbi(tx_data)
            if np.any(decoded_seq == 0):
                pol_wave = tx_indics[np.argmax(np.where(decoded_seq == 0)[0])]
            else:
                pol_wave = tx_indics[0]

            pol_wave *= jump
            # validity check
            assert (pol_wave >= tx.TSS and pol_wave <= tx.TES)
            # save relevant data for later processing
            g_hmm = [clusterId, tx.symbol, tx['name'], tx.TSS, tx.TES, pol_wave, cluster_chrom, p, cluster_strand]

            gene_hmm_data.append(g_hmm)
        print('---------------')

    col_names = ['clusterId', 'geneSymbol', 'name', 'TSS', 'TES', 'PolIIWave', 'chrom', 'p', 'strand']

    out_file = HMM_RESULT_DIR + 'gene_hmm_data{}{}.csv'.format(exp_name, suffix)
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    print('Saving output to %s' % out_file)
    pd.DataFrame(gene_hmm_data, columns=col_names).to_csv(out_file, index=False)


def download_clustersWithNames(genome='hg19'):
    from io import BytesIO
    print('Checking for metadata file in %s' % KNOWN_GENES)
    if os.path.exists(KNOWN_GENES):
        return
    else:
        print('downloading metadata from UCSC')
        print('Downloading knownGene table')
        known_gene = pd.read_csv(BytesIO(
            requests.get('http://hgdownload.soe.ucsc.edu/goldenPath/%s/database/knownGene.txt.gz' % genome).content),
            compression='gzip', sep='\t', names=[col % genome for col in [
                '#%s.knownGene.name', '%s.knownGene.chrom', '%s.knownGene.strand',
                '%s.knownGene.txStart', '%s.knownGene.txEnd',
                '%s.knownGene.exonStarts', '%s.knownGene.exonEnds']]
            , usecols=[0, 1, 2, 3, 4, 8, 9])
        print('Downloading kgxref table')
        kgxref = pd.read_csv(BytesIO(
            requests.get('http://hgdownload.soe.ucsc.edu/goldenPath/%s/database/kgXref.txt.gz' % genome).content),
            compression='gzip', sep='\t', usecols=[0, 4],
            names=['#%s.knownGene.name' % genome, '%s.kgXref.geneSymbol' % genome])
        print('Downloading knownIsoforms table')
        isoforms = pd.read_csv(
            BytesIO(requests.get(
                'http://hgdownload.soe.ucsc.edu/goldenPath/%s/database/knownIsoforms.txt.gz' % genome).content),
            compression='gzip', sep='\t', names=['%s.knownIsoforms.clusterId' % genome, '#%s.knownGene.name' % genome])
        known_gene = pd.merge(known_gene, kgxref, on='#%s.knownGene.name' % genome)
        known_gene = pd.merge(known_gene, isoforms, on='#%s.knownGene.name' % genome)
        os.makedirs(os.path.dirname(os.path.dirname(KNOWN_GENES)), exist_ok=True)
        known_gene.to_csv(KNOWN_GENES, index=False, sep='\t')
    print('Done')


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='HMM to create segmentation to transcribed/non-transcribed for each transcript')
    parser.add_argument("cond", help='L or R', default='L')
    parser.add_argument("exp_time", help='Time of elongation (DRB wash until harvesting)', default='8')
    parser.add_argument("--rep_sample", help='Replicate name or None (for averaging)', default=None)
    parser.add_argument("--suffix", help='unique suffix to the output files', default="")
    parser.add_argument("--genome", help='Name of the genome in UCSC', default="hg19")
    args = parser.parse_args()
    fit_hmm_multi(args.cond, args.exp_time, [args.rep_sample], args.suffix, genome=args.genome)


if __name__ == '__main__':
    main()
