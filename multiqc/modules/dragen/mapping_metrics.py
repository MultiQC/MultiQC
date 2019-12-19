#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, beeswarm

# Initialise the logger
import logging
log = logging.getLogger(__name__)


class DragenMappingMetics(BaseMultiqcModule):
    def parse_mapping_metrics(self):
        all_data_by_sample = dict()
        general_info_by_sample_by_id = defaultdict(dict)

        for f in self.find_log_files('dragen/mapping_metrics'):
            data_by_sample, general_info = parse_mapping_metrics_file(f)

            if data_by_sample:
                for sn, data in data_by_sample.items():
                    if sn in all_data_by_sample:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                    self.add_data_source(f, section='stats')

                    all_data_by_sample[sn] = data

                    for k, v in general_info.items():
                        general_info_by_sample_by_id[k][sn] = v

        # Filter to strip out ignored sample names:
        all_data_by_sample = self.ignore_samples(all_data_by_sample)

        if not all_data_by_sample:
            return
        log.info("Found mapping metrics for {} samples".format(len(all_data_by_sample)))

        all_metric_names = set()
        for sn, sdata in all_data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        headers, beeswarm_keys = make_mapping_stats_headers(all_metric_names)

        all_data_by_sample = maybe_add_general_metrics(general_info_by_sample_by_id, all_data_by_sample)

        self.general_stats_addcols(all_data_by_sample, headers, 'Mapping metrics')
        # Make bargraph plots of mapped, dupped and paired reads
        self.__map_dup_read_chart(all_data_by_sample)
        self.__map_pair_read_chart(all_data_by_sample)

        self.add_section(
            name='Mapping metrics',
            anchor='dragen-mapping-metrics',
            description="This module parses the output from <code>Dragen mapping metrics</code>. All numbers in millions.",
            plot=beeswarm.plot(all_data_by_sample, beeswarm_keys, {'id': 'dragen-mapping-metrics-dp'})
        )

    def __map_dup_read_chart(self, data_by_sample):
        chart_data = dict()
        for sample_id, data in data_by_sample.items():
            if data['Number of unique & mapped reads (excl. duplicate marked reads)']\
                    + data['Number of duplicate marked reads']\
                    + data['Unmapped reads'] != data['Total reads in RG']:
                log.warning("sum of unique/duplicate/unmapped reads not matching total, "
                            "skipping mapping/duplicates percentages plot for: {}".format(sample_id))
            else:
                chart_data[sample_id] = data
        self.add_section(
            name='Mapped and duplicates percentages',
            anchor='dragen-mapping-dup-percentage',
            description='Mapping and duplicate reads from Dragen: uniqly mapped vs. duplicate vs. unmapped reads.',
            plot=bargraph.plot(chart_data, {
                'Number of unique & mapped reads (excl. duplicate marked reads)': {'color': '#437bb1', 'name': 'Mapped'},
                'Number of duplicate marked reads':                               {'color': '#f5a742', 'name': 'Duplicated'},
                'Unmapped reads':                                                 {'color': '#b1084c', 'name': 'Unmapped'},
            }, {
                'id': 'mapping_dup_percentage_plot',
                'title': 'Dragen mapping metrics: duplicate reads',
                'ylab': '# Reads',
                'cpswitch_counts_label': 'Number of reads'
            })
        )

    def __map_pair_read_chart(self, data_by_sample):
        chart_data = dict()
        for sample_id, data in data_by_sample.items():
            if data['Not properly paired reads (discordant)'] + data['Properly paired reads']\
                    + data['Singleton reads (itself mapped; mate unmapped)']\
                    + data['Unmapped reads'] != data['Total reads in RG']:
                log.warning("sum of unpaired/discordant/proppaired/unmapped reads not matching total, "
                            "skipping mapping/paired percentages plot for: {}".format(sample_id))
            else:
                chart_data[sample_id] = data
        self.add_section(
            name='Mapped and paired percentages',
            anchor='dragen-mapping-paired-percentage',
            description="Mapping and pairing read metrics from Dragen: properly paired vs. discordant vs. unpaired vs. unmapped reads.",
            plot=bargraph.plot(chart_data, {
                'Properly paired reads':                             {'color': '#00cc00', 'name': 'Paired, properly'},
                'Not properly paired reads (discordant)':            {'color': '#ff9900', 'name': 'Paired, discordant'},
                'Singleton reads (itself mapped; mate unmapped)':    {'color': '#ff33cc', 'name': 'Singleton'},
                'Unmapped reads':                                    {'color': '#b1084c', 'name': 'Unmapped'},
            }, {
                'id': 'mapping_paired_percentage_plot',
                'title': 'Dragen mapping metrics: paired reads',
                'ylab': '# Reads',
                'cpswitch_counts_label': 'Number of reads'
            })
        )


def maybe_add_general_metrics(general_info_by_sample_by_id, all_data_by_sample):
    if general_info_by_sample_by_id:
        for id_in_data, title, fmt, modify in general_metrics:
            v_by_sn = general_info_by_sample_by_id[id_in_data]
            if len(list(set(v_by_sn.values()))) > 1:
                # if general stats are different for different samples - adding them into the general stats table
                for sn, v in v_by_sn.items():
                    if v != 'NA':
                        all_data_by_sample[sn][id_in_data] = v
            else:
                # else - adding them into the report header:
                if config.report_header_info is None:
                    config.report_header_info = list()

                v = list(set(v_by_sn.values()))[0]
                if v != 'NA':
                    if modify:
                        v = modify(v)
                    config.report_header_info.append({title: fmt.format(v)})

    return all_data_by_sample


def parse_mapping_metrics_file(f):
    """
    T_SRR7890936_50pc.mapping_metrics.csv

    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Total reads in RG,1100000000,100.00
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Number of duplicate marked reads,NA
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Number of duplicate marked and mate reads removed,NA
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Number of unique reads (excl. duplicate marked reads),NA
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Reads with mate sequenced,1100000000,100.00
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Reads without mate sequenced,0,0.00
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,QC-failed reads,0,0.00
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mapped reads,1066826874,96.98
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mapped reads R1,533963768,97.08
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mapped reads R2,532863106,96.88
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Number of unique & mapped reads (excl. duplicate marked reads),1066826874,96.98
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Unmapped reads,33173126,3.02
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Singleton reads (itself mapped; mate unmapped),1585054,0.14
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Paired reads (itself & mate mapped),1065241820,96.84
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Properly paired reads,1055464776,95.95
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Not properly paired reads (discordant),9777044,0.89
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Paired reads mapped to different chromosomes,6948576,0.65
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Paired reads mapped to different chromosomes (MAPQ>=10),3629346,0.34
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Reads with indel R1,13131067,2.46
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Reads with indel R2,12308296,2.31
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Total bases,165000000000
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Total bases R1,82500000000
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Total bases R2,82500000000
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mapped bases R1,80094565200
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mapped bases R2,79929465900
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Soft-clipped bases R1,700070179,0.87
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Soft-clipped bases R2,1419400306,1.77
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mismatched bases R1,294215249,0.37
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mismatched bases R2,610527210,0.76
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mismatched bases R1 (excl. indels),252464476,0.32
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Mismatched bases R2 (excl. indels),567877112,0.71
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Q30 bases,147576651866,89.44
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Q30 bases R1,77087948549,93.44
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Q30 bases R2,70488703317,85.44
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Q30 bases (excl. dups & clipped bases),150322268169
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Total alignments,1094089541
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Secondary alignments,0
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Supplementary (chimeric) alignments,27262667
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Estimated read length,150.00
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Average sequenced coverage over genome,51.41
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Insert length: mean,383.14
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Insert length: median,377.00
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Insert length: standard deviation,82.67
    NORMAL MAPPING/ALIGNING SUMMARY,,Bases in reference genome,3209286105
    NORMAL MAPPING/ALIGNING SUMMARY,,Bases in target bed [% of genome],NA
    NORMAL MAPPING/ALIGNING SUMMARY,,Provided sex chromosome ploidy,NA
    NORMAL MAPPING/ALIGNING SUMMARY,,DRAGEN mapping rate [mil. reads/second],0.53
    """

    data_by_sample = defaultdict(dict)
    general_info = dict()
    for line in f['f'].splitlines():
        fields = line.split(',')
        analysis = fields[0].split('/')[1]  # ALIGNING SUMMARY, ALIGNING PER RG
        # sample-unspecific metrics are reported only in ALIGNING SUMMARY sections
        if analysis == 'ALIGNING SUMMARY':
            metric = fields[2]
            if metric in [id_in_data for (id_in_data, title, fmt, modify) in general_metrics]:
                value = fields[3]
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                general_info[metric] = value

        # for sample-specific metrics, using ALIGNING PER RG because it has the sample name in the 2nd col
        if analysis == 'ALIGNING PER RG':
            sample = fields[1]
            metric = fields[2]
            value = fields[3]
            try:
                value = int(value)
            except ValueError:
                pass
            data_by_sample[sample][metric] = value

            if len(fields) > 4:  # percentage
                percentage = fields[4]
                try:
                    percentage = float(percentage)
                except ValueError:
                    pass
                data_by_sample[sample][metric + ' pct'] = percentage

    # adding some missing values that we wanna report for consistency
    for sname, data in data_by_sample.items():
        # fixing when deduplication wasn't performed
        if data['Number of duplicate marked reads'] == 'NA':
            data['Number of duplicate marked reads'] = 0
        if data['Number of duplicate marked and mate reads removed'] == 'NA':
            data['Number of duplicate marked and mate reads removed'] = 0
        if data['Number of unique reads (excl. duplicate marked reads)'] == 'NA':
            data['Number of unique reads (excl. duplicate marked reads)'] = data['Mapped reads']
        # adding alignment percentages
        if data['Total alignments'] > 0:
            data['Secondary alignments pct'] = data['Secondary alignments'] / data['Total alignments'] * 100.0
        # adding some missing bases percentages
        if data['Total bases'] > 0:
            data['Q30 bases (excl. dups & clipped bases) pct'] = data['Q30 bases (excl. dups & clipped bases)'] / data['Total bases'] * 100.0
            data['Mapped bases R1 pct'] = data['Mapped bases R1'] / data['Total bases'] * 100.0
            data['Mapped bases R2 pct'] = data['Mapped bases R2'] / data['Total bases'] * 100.0

    # removing NA values
    for sname, data in data_by_sample.items():
        for k, v in data.items():
            if v == 'NA':
                del data[k]

    return data_by_sample, general_info


read_format = '{:,.1f}'
if config.read_count_multiplier == 1:
    read_format = '{:,.0f}'
read_format += '&nbsp;' + config.read_count_prefix

base_format = '{:,.1f}&nbsp;'
if config.base_count_multiplier == 1:
    base_format = '{:,.0f}'
elif config.base_count_multiplier == 0.000000001:
    base_format = '{:,.2f}'
base_format += '&nbsp;' + config.base_count_prefix

mapping_metrics = [
    # id_in_data                                               # title (display name)        # show  # unit  # beeswarm  # description
    # Read stats:
    ('Total reads in RG'                                      , 'Reads'                      , '#',  'reads', True,  'Total number of reads in this read group (sample), {}'),
    ('Reads with mate sequenced'                              , 'Reads with mate'            , None, 'reads', True,  'Number of reads with a mate sequenced, {}'),
    ('Reads without mate sequenced'                           , 'Reads w/o mate'             , None, 'reads', False, 'Number of reads without a mate sequenced, {}'),
    ('QC-failed reads'                                        , 'QC-fail'                    , None, 'reads', True,  'Number of reads not passing platform/vendor quality checks (SAM flag 0x200), {}'),
    ('Mapped reads'                                           , 'Map'                        , '%',  'reads', True,  'Number of mapped reads, {}'),
    ('Mapped reads R1'                                        , 'Map R1'                     , None, 'reads', False, 'Number of mapped reads R1, {}'),
    ('Mapped reads R2'                                        , 'Map R2'                     , None, 'reads', False, 'Number of mapped reads R2, {}'),
    ('Reads with MAPQ [40:inf)'                               , 'MQ>=40'                     , None, 'reads', True,  'Number of reads with MAPQ [40:inf), {}'),
    ('Number of duplicate marked reads'                       , 'Dup'                        , '%',  'reads', False, 'Number of duplicate marked reads as a result of the --enable-duplicatemarking option being used, {}'),
    ('Number of duplicate marked and mate reads removed'      , 'Dup with mates removed'     , None, 'reads', False, 'Number of reads marked as duplicates, along with any mate reads, that are removed '
                                                                                                                        'when the --remove-duplicates option is used, {}'),
    ('Number of unique reads (excl. duplicate marked reads)'  , 'Uniq'                       , None, 'reads', True,  'Number of unique reads (all reads minus duplicate marked reads), {}'),
    ('Number of unique & mapped reads (excl. duplicate marked reads)', 'Uniq map'            , None, 'reads', True,  'Number of unique & mapped reads (mapped reads minus duplicate marked reads), {}'),
    ('Unmapped reads'                                         , 'Unmap'                      , None, 'reads', False, 'Number of unmapped reads, {}'),
    ('Paired reads (itself & mate mapped)'                    , 'Self & mate mapped'         , None, 'reads', True,  'Number of reads mapped in pairs (itself & mate mapped), {}'),
    ('Properly paired reads'                                  , 'Prop pair'                  , '%',  'reads', True,  'Number of properly paired reads, {} (both reads in pair are mapped and '
                                                                                                                        'fall within an acceptable range from each other based on the estimated insert length distribution)'),
    ('Not properly paired reads (discordant)'                 , 'Discord'                    , None, 'reads', True,  'Number of discordant reads: paired reads minus properly paired reads , {}'),
    ('Singleton reads (itself mapped; mate unmapped)'         , 'Singleton'                  , None, 'reads', True,  'Number of singleton reads: itself mapped; mate unmapped, {}'),
    ('Paired reads mapped to different chromosomes'           , 'Mate map to diff chrom'     , None, 'reads', True,  'Number of paired reads with a mate mapped to a different chromosome, {}'),
    ('Paired reads mapped to different chromosomes (MAPQ>=10)', 'Mate diff chrom, MQ>=10'    , None, 'reads', True,  'Number of paired reads, mapped with MAPQ>=10 and with a mate mapped to a different chromosome, {}'),
    ('Reads with indel R1'                                    , 'Reads with indel R1'        , None, 'reads', False, 'Number of R1 reads containing at least 1 indel, {}'),
    ('Reads with indel R2'                                    , 'Reads with indel R2'        , None, 'reads', False, 'Number of R2 reads containing at least 1 indel, {}'),
    # Alignments stats:
    ('Total alignments'                                       , 'Alignments'                 , None, 'reads', True,  'Total number of alignments with MQ > 0, {}'),
    ('Secondary alignments'                                   , 'Second\'ry'                 , '%',  'reads', True,  'Number of secondary alignments, {}. Secondary alignment occurs when '
                                                                                                                        'a given read could align reasonably well to more than one place. '
                                                                                                                        'One of the possible reported alignments is termed "primary" and '
                                                                                                                        'the others will be marked as "secondary".'),
    ('Supplementary (chimeric) alignments'                    , 'Suppl\'ry'                  , None, 'reads', True,  'Number of supplementary (chimeric) alignments, {}. A chimeric read is split '
                                                                                                                        'over multiple loci (possibly due to structural variants). One alignment is '
                                                                                                                        'referred to as the representative alignment, the other are supplementary'),
    # Read length stats:
    ('Estimated read length'                                  , 'Read len'                   , '#',  'len',   False, 'Estimated read length. Total number of input bases divided by the number of reads'),
    ('Insert length: mean'                                    , 'Avg IS'                     , None, 'len',   False, 'Insert length: mean'),
    ('Insert length: median'                                  , 'Med IS'                     , '#',  'len',   False, 'Insert length: median'),
    ('Insert length: standard deviation'                      , 'IS std'                     , None, 'len',   False, 'Insert length: standard deviation'),
    # Coverage:
    ('Average sequenced coverage over genome'                 , 'Cov'                        , '#',  'x',     True,  'Average sequenced coverage over genome'),
    # Bases stats:
    ('Total bases'                                            , 'Bases'                      , '#',  'bases', True,  'Total number of bases sequenced, {}'),
    ('Total bases R1'                                         , 'Bases R1'                   , None, 'bases', False, 'Total number of bases sequenced on R1 reads, {}'),
    ('Total bases R2'                                         , 'Bases R2'                   , None, 'bases', False, 'Total number of bases sequenced on R2 reads, {}'),
    ('Mapped bases R1'                                        , 'Mapped bases R1'            , None, 'bases', False, 'Number of mapped bases on R1 reads, {}'),
    ('Mapped bases R2'                                        , 'Mapped bases R2'            , None, 'bases', False, 'Number of mapped bases on R2 reads, {}'),
    ('Soft-clipped bases R1'                                  , 'Soft-clip bases R1'         , None, 'bases', False, 'Number of soft-clipped bases on R1 reads, {}'),
    ('Soft-clipped bases R2'                                  , 'Soft-clip bases R2'         , None, 'bases', False, 'Number of soft-clipped bases on R2 reads, {}'),
    ('Mismatched bases R1'                                    , 'MM bases R1'                , None, 'bases', False, 'Number of mismatched bases on R1, {}, which is the sum of SNP count and indel lengths. '
                                                                                                                        'It does not count anything within soft clipping, or RNA introns. '
                                                                                                                        'It also does not count a mismatch if either the reference base or read base is N'),
    ('Mismatched bases R2'                                    , 'MM bases R2'                , None, 'bases', False, 'Number of mismatched bases on R2, {}, which is the sum of SNP count and indel lengths. '
                                                                                                                        'It does not count anything within soft clipping, or RNA introns. '
                                                                                                                        'It also does not count a mismatch if either the reference base or read base is N'),
    ('Mismatched bases R1 (excl. indels)'                     , 'MM bases R1 excl indels'    , None, 'bases', False, 'Number of mismatched bases on R1, {}. The indels lengts are ignored. '
                                                                                                                        'It does not count anything within soft clipping, or RNA introns. '
                                                                                                                        'It also does not count a mismatch if either the reference base or read base is N'),
    ('Mismatched bases R2 (excl. indels)'                     , 'MM bases R2 excl indels'    , None, 'bases', False, 'Number of mismatched bases on R2, {}. The indels lengts are ignored. '
                                                                                                                        'It does not count anything within soft clipping, or RNA introns. '
                                                                                                                        'It also does not count a mismatch if either the reference base or read base is N'),
    ('Q30 bases'                                              , 'Q30'                        , '%',  'bases', True,  'Number of bases with BQ >= 30, {}'),
    ('Q30 bases R1'                                           , 'Q30 R1'                     , None, 'bases', False, 'Number of bases on R1 reads with BQ >= 30, {}'),
    ('Q30 bases R2'                                           , 'Q30 R2'                     , None, 'bases', False, 'Number of bases on R2 reads with BQ >= 30, {}'),
    ('Q30 bases (excl. dups & clipped bases)'                 , 'Q30 excl dup & clipped'     , None, 'bases', True,  'Number of non-clipped bases with BQ >= 30 on non-duplicate reads, {}'),
    # General metrics. Showing only when general metrics are different for different samples, otherwise showing in the header
    ('Bases in reference genome'                              , 'Bases in ref. genome'       , '#',  'bases', False, 'Bases in reference genome'              ),
    ('Bases in target bed [% of genome]'                      , 'Bases in target bed'        , '#',  '%'    , False, 'Bases in target bed [% of genome]'      ),
    ('Provided sex chromosome ploidy'                         , 'Provided sex chrom ploidy'  , None,  None  , False, 'Provided sex chromosome ploidy'         ),
    ('DRAGEN mapping rate [mil. reads/second]'                , 'DRAGEN map rate'            , None,  None  , False, 'DRAGEN mapping rate [mil. reads/second]'),
]

general_metrics = [
    # id_in_data                              # title                            # format                         # modify
    ('Bases in reference genome'              , 'Bases in ref. genome'           , base_format                    , lambda v: v * config.base_count_multiplier ),
    ('Bases in target bed [% of genome]'      , 'Bases in target bed'            , '{} % of genome'               , None                                       ),
    ('Provided sex chromosome ploidy'         , 'Provided sex chrom ploidy'      , '{}'                           , None                                       ),
    ('DRAGEN mapping rate [mil. reads/second]', 'DRAGEN mapping rate'            , '{:,.2f} [mil. reads/second]'  , None                                       ),
]

def make_mapping_stats_headers(metric_names):
    # Init general stats table
    stats_headers = OrderedDict()

    # Init beeswarm plot
    beeswarm_keys = OrderedDict()

    for id_in_data, title, showing, unit, show_in_beeswarm, descr in mapping_metrics:
        col = dict(
            title=title,
            description=descr,
            min=0,
        )
        if unit == 'reads':
            col['scale'] = 'RdYlGn'
        if unit == 'bases':
            col['scale'] = 'RdBu'
        if unit == 'len':
            col['scale'] = 'BrBG'
        if unit == 'x':
            col['scale'] = 'PiYG'

        if id_in_data + ' pct' in metric_names:
            # if % value is available, showing it instead of the number value; the number value will be hidden
            pct_col = dict(
                col,
                description=descr.replace(', {}', '').replace('Number of ', '% of '),
                max=100,
                suffix='%',
                hidden=showing != '%'
            )
            stats_headers[id_in_data + ' pct'] = pct_col

        col['hidden'] = showing != '#'
        if unit == 'reads':
            col['description'] = col['description'].format(config.read_count_desc)
            col['modify'] = lambda x: x * config.read_count_multiplier
            col['shared_key'] = 'read_count'
            col['format'] = read_format
        if unit == 'bases':
            col['description'] = col['description'].format(config.base_count_desc)
            col['modify'] = lambda x: x * config.base_count_multiplier
            col['shared_key'] = 'base_count'
            col['format'] = base_format
        if unit == 'len':
            col['suffix'] = ' bp'
            col['format'] = '{:,.0f}'
        if unit == 'x':
            col['suffix'] = ' x'
            col['format'] = '{:,.1f}'
        if unit == '%':
            col['suffix'] = ' %'
            col['format'] = '{:,.1f}'

        stats_headers[id_in_data] = col

        if show_in_beeswarm:
            suffix = ''
            if unit == 'reads':
                suffix = ' ' + config.read_count_prefix
            if unit == 'bases':
                suffix = ' ' + config.base_count_prefix
            beeswarm_keys[id_in_data] = dict(col,
                suffix=suffix
            )

    return stats_headers, beeswarm_keys






























