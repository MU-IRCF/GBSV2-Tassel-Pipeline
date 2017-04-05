# NAME

GBSV2::Tassel::Pipeline  

# VERSION

version 0.0006

# DESCRIPTION

Run Tassel's GBSv2 pipeline using parameters from a json configuration file  

# SYNOPSIS

    GBSV2/Tassel/Pipeline.pm config.json  

Typical `config.json`:

    { 
         "MINTAG"             :           "1",
         "NAME"               :       "MyNIL",
         "NUM_OF_FASTQ_FILES" :           "6",
         "ENZYME_OR_ENZYMES"  :   "PstI-MspI",
         "KEYFILE"            : "keyfile.txt",
         "DATABASE"           :         "RAD",
         "FASTQ_DIR"          :       "fastq",
         "STEPS"              : ["100_GBSToTag",
                                 "200_TagToFASTQ",
                                 "300_bowtie2",
                                 "400_SAMToGBS",
                                 "500_DiscoverSNP",
                                 "600_SNPQuality",
                                 "700_ProductionSNP"]
    }

# DIAGNOSTICS

    None implemented  

# INCOMPATIBILITIES

    None known  

# BUGS AND LIMITATIONS

    There are no known bugs in this module.  
    Please report problems to the author.  
    Patches are welcome.  

# DEPENDENCIES

## Runtime requires  
    Perl 5.10.1 or later  
    SLURM scheduler  
    Linux  

## Tests require  
    Test2::Bundle::Extended  
    File::Slurper  
