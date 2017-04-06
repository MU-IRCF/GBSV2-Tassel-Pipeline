# NAME

GBSV2::Tassel::Pipeline  

# VERSION

version 0.0009

# DESCRIPTION

Run Tassel's GBSv2 pipeline using parameters from a json configuration file  

# SYNOPSIS

    GBSV2/Tassel/Pipeline.pm config.json  

Typical `config.json`:

    { 
         "MINTAGS"            : ["1","5","10"],
         "NAME"               :        "MyNIL",
         "ENZYME_OR_ENZYMES"  :    "PstI-MspI",
         "KEYFILE"            :  "keyfile.txt",
         "FASTQ_DIR"          :        "fastq"
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
