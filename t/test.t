#!/bin/env perl
use strict;
use warnings;
use autodie;
use v5.10;

use Test2::Bundle::Extended;
use Test::LongString;

use lib 'lib';

use GBS::Tassel::PipelineV2 qw(get_script_for read_config);

my %config = read_config('t/pipeline.json'); 

for my $step ( @{ $config{STEPS} } ) {
    my ($number, $lil_name) = split '_', $step;

    for my $minTag ( @{$config{MINTAGS}} ) {
        is( get_script_for(STEP   => $step,
                           MINTAG => $minTag,
                           %config), 
            expected_for($number),
            "Step $number script correct (Test2 comparison)"
        );
        is_string( get_script_for(STEP   => $step,
                           MINTAG => $minTag,
                           %config), 
            expected_for($number),
            "Step $number script correct (Test::LongString comparison)"
        );
    }
}

is("THIS IS A TEST\n", expected_for('minitest'), 'mini test worked');

done_testing();

sub expected_for {
    my $step = shift;

    my %expected_for = (

'minitest' => <<'END',
THIS IS A TEST
END

'100' => <<'END',
#!/bin/bash
#SBATCH -J 100_GBSToTag_1minTag_MyNIL
#SBATCH -o 100_GBSToTag_1minTag_MyNIL.o_%j
#SBATCH -e 100_GBSToTag_1minTag_MyNIL.e_%j
#SBATCH --partition=BioCompute,Lewis
#SBATCH --nodes=1
#SBATCH --mem=460G
#SBATCH --cpus-per-task=6
#SBATCH --time=1-0

module load Tassel/tassel-5.2.35

run_pipeline.pl -Xms200G -Xmx450G  \
    -fork1                         \
        -GBSSeqToTagDBPlugin       \
            -k  keyfile.txt        \
            -db GBS.db             \
            -c  1                  \
            -e  PstI-MspI          \
            -i  fastq              \
            -mxKmerNum 500000000   \
        -endPlugin                 \
    -runfork1

# Defaults
#            -c 10                                                                                                                \
#            -mnQS 0                                                                                                              \
#            -kmerLength 64                                                                                                       \
#            -minKmerL 20                                                                                                         \
#            -mxKmerNum 50000000                                                                                                  \
#            -batchSize 8                                                                                                         \
END

'200' => <<'END',
#!/bin/bash
#SBATCH -J 200_TagToFASTQ_1minTag_MyNIL
#SBATCH -o 200_TagToFASTQ_1minTag_MyNIL.o_%j
#SBATCH -e 200_TagToFASTQ_1minTag_MyNIL.e_%j
#SBATCH --partition=BioCompute,Lewis
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --time=01:00:00

module load Tassel/tassel-5.2.35

# Create link to FASTQ directory
ln -s ../fastq

# Copy database
cp ../GBS.db ./GBS.db

run_pipeline.pl -Xms5G -Xmx5G       \
    -fork1                          \
        -TagExportToFastqPlugin     \
            -c  1                   \
            -db GBS.db              \
            -o  1minTag.fastq       \
        -endPlugin                  \
    -runfork1                       \

# Defaults
#            -c 1

sbatch cmd_300_bowtie2_1minTag_MyNIL.sbatch
END

'300' => <<'END',
#!/bin/bash
#SBATCH -J 300_bowtie2_1minTag_MyNIL
#SBATCH -o 300_bowtie2_1minTag_MyNIL.o_%j
#SBATCH -e 300_bowtie2_1minTag_MyNIL.e_%j
#SBATCH --partition=BioCompute,Lewis
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=10:00:00

module load bowtie2/bowtie2-2.3.1

bowtie2 --no-unal --un-gz unaligned.fastq --threads 12 --reorder -x /group/ircf/dbase/genomes/maize/maizesequence.org/AGPv2/bowtie2.3.1-index/AGPv2 1minTag.fastq > 1minTag.sam

sbatch cmd_400_SAMToGBS_1minTag_MyNIL.sbatch
END

'400' => <<'END',
#!/bin/bash
#SBATCH -J 400_SAMToGBS_1minTag_MyNIL
#SBATCH -o 400_SAMToGBS_1minTag_MyNIL.o_%j
#SBATCH -e 400_SAMToGBS_1minTag_MyNIL.e_%j
#SBATCH --partition=BioCompute,Lewis
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --time=01:00:00

module load Tassel/tassel-5.2.35

run_pipeline.pl -Xms40G -Xmx40G  \
    -fork1                       \
        -SAMToGBSdbPlugin        \
            -db GBS.db           \
            -i  1minTag.sam      \
        -endPlugin               \
    -runfork1                    \

# Defaults
# -aProp 0.0
# -aLen 0

sbatch cmd_500_DiscoverySNP_1minTag_MyNIL.sbatch
END

'500' => <<'END',
#!/bin/bash
#SBATCH -J 500_DiscoverySNP_1minTag_MyNIL
#SBATCH -o 500_DiscoverySNP_1minTag_MyNIL.o_%j
#SBATCH -e 500_DiscoverySNP_1minTag_MyNIL.e_%j
#SBATCH --partition=BioCompute,Lewis
#SBATCH --nodes=1
#SBATCH --mem=160G
#SBATCH --nodes=1
#SBATCH --time=2-0

module load Tassel/tassel-5.2.35

run_pipeline.pl -Xms120G -Xmx159G    \
    -fork1                           \
        -DiscoverySNPCallerPluginV2  \
            -db GBS.db               \
        -endPlugin                   \
    -runfork1                        \

#The parameters to this plugin are:
# -callBiSNPsWGap <true | false> : Include sites where the third allele is a GAP (mutually exclusive with inclGaps) (Default: false) This option has not yet been implemented.
# -db <Input GBS Database> : Input Database file if using SQLite (REQUIRED)
# -gapAlignRatio <Gap Alignment Threshold> : Gap alignment threshold ratio of indel contrasts to non indel contrasts: IC/(IC + NC). Tags will be excluded from any loci that has a tag with a gap ratio that exceeds the threshold.
# -inclGaps <true | false > : Include sites where major or minor allele is a GAP (Default: false) This option has not yet been implemented.
# -inclRare <true | false> : Include the rare alleles at site (3 or 4th states) (Default: false)
# -maxTagsPerCutSite <Max tags per cut site position> : Maximum number of tags allowed per cut site when aligning tags . All cut site positions and their tags are stored in the database, but alignment of tags at a particular cut site position only occurs when the number of tags at this position is equal to or less than maxTagsPerCutSite. This guards against software degradation when a position has hundreds or thousands of associated tags. (Default: 64)
# -mnLCov <Min Locus Coverage> : Minimum locus coverage (proportion of Taxa with a genotype) (Default: 0.1)
# -mnMAF <Min Minor Allele Freq> : Minimum minor allele frequency (Default: 0.01)
# -ref <Reference Genome File> : Path to reference genome in fasta format. Ensures that a tag from the reference genome is always included when the tags at a locus are aligned against each other to call SNPs. (Default: Don't use reference genome)
# -sC <Start Chromosome> : Start Chromosome: If missing, processing starts with the first chromosome (lexicographically) in the database.
# -eC <End Chromosome> : End Chromosome : If missing, plugin processing ends with the last chromosome (lexicographically) in the database.
# -deleteOldData <true | false> : Whether to delete old SNP data from the data bases. If true, all data base tables previously populated from the DiscoverySNPCallerPluginV2 and later steps in the GBSv2 pipeline is deleted. This allows for calling new SNPs with different pipeline parameters. (Default: false)

sbatch cmd_600_SNPQuality_1minTag_MyNIL.sbatch
END

'600' => <<'END',
#!/bin/bash
#SBATCH -J 600_SNPQuality_1minTag_MyNIL
#SBATCH -o 600_SNPQuality_1minTag_MyNIL.o_%j
#SBATCH -e 600_SNPQuality_1minTag_MyNIL.e_%j
#SBATCH --partition=BioCompute,Lewis
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --time=10:00:00

module load Tassel/tassel-5.2.35

run_pipeline.pl -Xms60G -Xmx60G         \
    -fork1                              \
        -SNPQualityProfilerPlugin       \
            -db GBS.db                  \
            -statFile snpQualStats.txt  \
        -endPlugin                      \
    -runfork1                           \

# The parameters to this plugin are:
# -db <Output Database file> : Name of output file (e.g. GBSv2.db) (REQUIRED)
# -deleteOldData <true | false> : Delete existing SNP quality data from db tables. If true, SNP quality data is cleared from the snpQuality table in the database before adding new data. (Default: false)
# -taxa <Taxa sublist file> : Name of taxa list input file in taxa list format
# -tname <Taxa set name > : Name of taxa set for database.
# -statFile <Quality information output name > : Name of the output file containing the quality statistics in a tab delimited format.

sbatch cmd_700_ProductionSNP_1minTag_MyNIL.sbatch
END

'700' => <<'END',
#!/bin/bash
#SBATCH -J 700_ProductionSNP_1minTag_MyNIL
#SBATCH -o 700_ProductionSNP_1minTag_MyNIL.o_%j
#SBATCH -e 700_ProductionSNP_1minTag_MyNIL.e_%j
#SBATCH --partition=BioCompute,Lewis
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=6
#SBATCH --time=2-0

module load Tassel/tassel-5.2.35

run_pipeline.pl -Xms10G -Xmx100G      \
    -fork1                            \
        -ProductionSNPCallerPluginV2  \
            -k  ../keyfile.txt           \
            -db GBS.db                \
            -e  PstI-MspI             \
            -i  fastq                 \
            -o  1mintag.hdf5          \
        -endPlugin                    \
    -runfork1                         \

# The parameters to this plugin are:
# -batchSize <Batch Size> : Number of flow cells to process simultaneously. (Default: 8)
# -d <Max Divergence> : Maximum divergence (edit distance) between new read and previously mapped read (Default: 0 = perfect matches only) (Default: 0)
# -db <Input GBS Database> : Input Database file if using SQLite (REQUIRED)
# -e <Enzyme> : Enzyme used to create the GBS library (REGQUIRED)
# -eR <Ave Seq Error Rate> : Average sequencing error rate per base (used to decide between heterozygous and homozygous calls) (Default: 0.01)
# -i <Input Directory> : Input directory containing fastq AND/OR qseq files (REQUIRED)
# -k <Key File> : Key file listing barcodes distinguishing the sample (REQUIRED)
# -ko <true | false> : Keep HDF5 genotypes file open for future runs that add more taxa or more depth (Default: false)
# -kmerLength <Length of Kmer> : Lemgth of kmers to grab from fastQ sequences. This value should match the kmerLength parameter value used in the GBSSeqToTagDBPlugin step of the pipeline. Bad values may be stored in the HDF5 file if these values are inconsistent. (Default: 64)
# -minPosQS < Minimum Quality Score> : The minimum quality score a SNP position must have to be output to the HDF5 file. The quality scores are loaded into the database via the tab-delimited file read in from the UpdateSNPPositionQualityPlugin. A value of 0 indicates all positions should be processed. (Default: 0)
# -mnQS < Minimum Quality Score> : The minimum quality score within the barcode and read length for a position to be accepted. This filters the read values from the fastQ files. (Default: 0)
# -do <true | false> : depth output: write depths to the output HDF5 genotypes file (Default: true)
# -o <Output HDF5 Genotypes File> : Output (target) HDF5 genotypes file to which is added new genotypes (new file created if it doesn't exist) (REQUIRED)

# This was the last step, so finished
END
);
    return $expected_for{$step};
}

