# RScript proteos_cleave_fasta_peptides_DIA_pipeline.R
# 2020-12-14
# (C) Signature Science, LLC 2020
# Myles Gardner
# Script to cleave proteins in a FASTA file into peptides
# using specified list of enzymes, filtering for unique peptides
# and peptides within a certain mass range.

# RScript --vanilla proteos_cleave_fasta_peptides_DIA_pipeline.R
# -f "L:/Proteos/fasta/human/Homo_sapiens_GRCh38_pep_all_proteomics_target.fasta"
# -c 2 -l 7 -L 30 -M "TRUE"
# -e "trypsin,trypsin-p"
# -E "U:/Projects/proteos_pipelines/config/enzyme_ref.csv"
# -o "Y:/fasta/human/"

`%notin%` <- function(x,y) !(x %in% y)

options(warn=-1)

check_libraries <- function() {
    # libraries to check for and install if not available
    r_packages <- c("Peptides", "stringr", "optparse", 'stringi',
                    "BiocManager", 'seqinr', "data.table", 'R.utils', 'yaml')
    bio_packages <- c("cleaver")
    new_packages <- r_packages[!(r_packages %in% installed.packages()[,"Package"])]
    if(length(new_packages)) install.packages(new_packages)
    library(BiocManager, warn.conflicts = FALSE, quietly = TRUE)
    new_packages <- bio_packages[!(bio_packages %in% installed.packages()[,"Package"])]
    if(length(new_packages)) {
        #source("https://bioconductor.org/biocLite.R")
        #biocLite(new_packages)
        BiocManager::install(new_packages)
    }
    # libraries to import    
    library(tools, warn.conflicts = FALSE, quietly = TRUE)
    library(R.utils, warn.conflicts = FALSE, quietly = TRUE)
    library(optparse, warn.conflicts = FALSE, quietly = TRUE)
    library(cleaver, warn.conflicts = FALSE, quietly = TRUE)
    library(Peptides, warn.conflicts = FALSE, quietly = TRUE)
    library(stringr, warn.conflicts = FALSE, quietly = TRUE)
    library(stringi, warn.conflicts = FALSE, quietly = TRUE)
    library(seqinr, warn.conflicts = FALSE, quietly = TRUE)
    library(yaml, warn.conflicts = FALSE, quietly = TRUE)
    # library(itertools, warn.conflicts = FALSE, quietly = TRUE)
    library(data.table, warn.conflicts = FALSE, quietly = TRUE)
    # library(foreach)
    # library(doSNOW)
}

get_os <- function() {
    if (.Platform$OS.type == "windows") { 
        rep_str <<- "\\\\"
        os <<- "windows"
    } else if (Sys.info()["sysname"] == "Darwin") {
        rep_str <<- "/"# "mac" 
        os <<- "mac"
    } else if (.Platform$OS.type == "unix") { 
        rep_str <<- "/"# "unix" 
        os <<- "linux"
    } else {
        rep_str <<- "/"
        os <<- 'unknown'
        # stop("Unknown OS")
    }
}

check_libraries()
get_os()
# threads <<- detectCores()

default_val <- list(
    missed_cleavages = 2,
    fasta_file = "Y:/fasta/n25_exomes/SA001.fa.gz",
    enzyme_file = "Y:/config/enzyme_ref.csv",
    length_min = 7,
    length_max = 50,
    Nterm_Met_cleavage = TRUE,
    enzymes = 'trypsin;trypsin-p',
    protein_donor_id_sep = '_',
    donor_id = NULL,
    output_dir = "Y:/fasta/HG38Proteomes"
)

get_args <- function(default_val, enzymes_allowed) {
    option_list = list(
        make_option(c("-f", "--fasta_file"), type="character", default=default_val$fasta_file,
                    help=paste("Path to fasta file to be cleaved.",
                               'Paths with spaces must be enclosed in double quotations ""'),
                    metavar="character"),
        make_option(c("-e", "--enzymes"), type="character", default=default_val$enzymes,
                    help=paste('Enzyme name(s). Allowable values include:',
                               paste(paste0('"', enzymes_allowed, '"'), collapse = "; "),
                               '\nEnzymes with spaces',
                               'must be enclosed in double quotations. Multiple enzymes\n',
                               'should be separated with commas or semicolons.')),
        make_option(c("-c", "--missed_cleavages"), type="integer", default=default_val$missed_cleavages,
                    help="number of missed cleavages"),
        make_option(c("-s", "--protein_donor_id_sep"), type="character", default=default_val$protein_donor_id_sep,
                    help="Character in FASTA protein identifier values that separates actual protein ID and donor ID."),
        make_option(c("-l", "--length_min"), type="integer", default=default_val$length_min,
                    help="minimum peptide length in amino acids"),
        make_option(c("-L", "--length_max"), type="integer", default=default_val$length_max,
                    help="maximum peptide length in amino acids"),
        make_option(c("-M", "--Nterm_Met_cleavage"), type="logical", default=default_val$Nterm_Met_cleavage,
                    help=paste("TRUE | FALSE value as to whether N-terminal methionine residues\n",
                               "should be cleaved prior to digestion in addition to not-cleaving\n",
                               "N-terminal methionines.")),
        make_option(c("-E", "--enzyme_file"), type="character", default=default_val$enzyme_file,
                    help=paste("Path to enzyme look-up csv file.",
                               'Paths with spaces must be enclosed in double quotations ""')),
        make_option(c("--donor_id"), type="character", default=default_val$donor_id,
                    help=paste("Accurate donor identifier value to use.", 
                               "If not supplied, will be inferred from protein sequences.",
                               "If supplied, will overwrite all protein sequences"),
                    metavar="character"),
        make_option(c("-o", "--output_dir"), type="character", default=default_val$output_dir,
                    help=paste("Path to output directory. Output file will be the same as the",
                               "input fasta file, appended with enzyme name(s) and min-max length.",
                               'Paths with spaces must be enclosed in double quotations ""'))
    ) # human_db
    opt_parser <- OptionParser(option_list=option_list)
    opt <- parse_args(opt_parser)
    opt$custom <- NULL
    return(opt)
}

check_args <- function(opt, enzymes_allowed, default_val) {
    # check for fasta file
    if (!file.exists(opt$fasta_file)) {
        stop(paste0("FASTA file not found: ", opt$fasta_file, "!"))
    }
    opt$fasta_file <- file_path_as_absolute(opt$fasta_file)
    # check for enzyme file
    if (!file.exists(opt$enzyme_file)) {
        stop(paste0("Enzyme reference lookup file not found: ", opt$enzyme_file, "!"))
    }
    opt$enzyme_file <- file_path_as_absolute(opt$enzyme_file)
    # check length_min
    if (opt$length_min > opt$length_max | opt$length_min < 5) {
        print(paste("Invalid length_min value:", opt$length_min, "! Setting to default value:",
                    default_val$length_min))
        opt$length_min <- default_val$length_min
    }
    # check length_max
    if (opt$length_max < opt$length_min | opt$length_max > 100) {
        print(paste("Invalid length_max value:", opt$length_max, "! Setting to default value:",
                    default_val$length_max))
        opt$length_max <- default_val$length_max
    }
    # check missed_cleavages
    if (opt$missed_cleavages < 0 | opt$missed_cleavages > 4) {
        print(paste("Invalid missed_cleavages value:", opt$missed_cleavages, "! Setting to default value:",
                    default_val$missed_cleavages))
        opt$missed_cleavages <- default_val$missed_cleavages
    }
    # check Nterm_Met_cleavage
    if (opt$Nterm_Met_cleavage %notin% c(TRUE, FALSE)) {
        print(paste("Invalid Nterm_Met_cleavage value:", opt$Nterm_Met_cleavage, "! Setting to default value:",
                    default_val$Nterm_Met_cleavage))
        opt$Nterm_Met_cleavage <- default_val$Nterm_Met_cleavage
    }
    # check enzymes_allowed
    enzyme_v <- sort(unique(unlist(str_split(opt$enzymes, pattern = c("[,;]")))))
    idx <- enzyme_v %notin% enzymes_allowed
    opt$enzymes <- enzyme_v[!idx]
    if (any(idx)) {
        print(paste("The following enzyme values are invalid:",
                    paste(enzyme_v[idx], collapse = ';'), "!"))
        if (all(idx)) {
            print(paste("Setting to default value:", paste(default_val$enzymes, collapse = ';'), "."))
            opt$enzymes <- default_val$enzymes
        } else {
            print(paste("Setting to:", paste(enzyme_v[!idx], collapse = ';'), "."))
        }
    }
    # check for output directory and set full path to output file
    if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = T)
    return(opt)
}

parse_fasta <- function(fasta_file) {
    start_char <- readLines(fasta_file, n = 1)
    str_sub(start_char, 1, 1) <- ""
    fasta_source <- ifelse(grepl("^(((tr)|(sp)|(sap))\\|)", start_char), 'uniprot',
                           ifelse(grepl("^(ens\\|)", start_char), 'ensembl_mod', 'ensembl'))
    print(paste0("Reading FASTA database ", fasta_file))
    raw_fasta <- read.fasta(fasta_file, seqtype = "AA", as.string = T, strip.desc = T)
    print(paste0("Parsing FASTA database ", fasta_file))
    annot <- as.character(sapply(getAnnot(raw_fasta), function(x){x[1:length(x)]}))
    raw_fasta <- as.data.frame(sapply(raw_fasta, function(x){x[1:length(x)]}), stringsAsFactors = F)
    nm <- row.names(raw_fasta)
    names(raw_fasta) <- "aa_sequence"
    if (fasta_source != 'ensembl') {
        # uniprot or uniprot-esque files with '>sp" or '>tr' headers
        # and '|' separators
        nmSplit <- str_split_fixed(nm, "\\|", 3)
        descSplit <- str_split_fixed(str_split_fixed(annot, "\\|", 3)[,3], " ", 2)
        rm(annot)
        raw_fasta <- setDT(raw_fasta)
        raw_fasta$fasta_header <- nmSplit[, 1]
        raw_fasta$protein_id <- nmSplit[, 2]
        if (fasta_source != 'ensembl_mod') {
            raw_fasta$gene_id <- nmSplit[, 3]
            raw_fasta$description <- descSplit[, 2]
            raw_fasta[grepl("OS\\=", description), `:=`(species_name = gsub("( OX\\=)|(OS\\=)", "", str_extract(description, "OS\\=.*OX\\=")),
                                                        species_number = as.numeric(gsub("OX\\=", "", str_extract(description, "OX\\=[0-9]+"))),
                                                        protein_name = gsub(" OS\\=", "", str_extract(description, ".* OS\\=")))]
            rm(nmSplit)
        } else {
            # MWG modified ensembl files to initially look like uniprot with
            # '|' separators and an '>ens' header, as opposed to no header and space separators
            # in actual ensembl fasta files
            rm(nmSplit)
            descSplit <- paste0(" ", descSplit[, 2])
            dt <- do.call(cbind, lapply(c("chromosome", "gene", "transcript", 
                                          "description", "gene_biotype", 
                                          "transcript_biotype", "gene_symbol"), 
                                        function (j) {
                                            # search for and add fasta description information
                                            regex_pattern <- ifelse(j != 'description',
                                                                    paste0(" ", j, "\\:[^ ]*"),
                                                                    paste0(" ", j, "\\:(.)*"))
                                            x <- gsub(paste0(" ", j, "\\:"), "", 
                                                      str_extract(descSplit, regex_pattern))
                                            x <- setDT(data.frame(x, stringsAsFactors = F))
                                            names(x) <- j
                                            x
                                            }))
            rm(descSplit)
            raw_fasta <- cbind(raw_fasta, dt)
            rm(dt)
        }
    } else {
        # ensembl fasta files with no header'>' and space separators between metadata
        raw_fasta$protein_id <- nm
        raw_fasta <- setDT(raw_fasta)
        dt <- do.call(cbind, lapply(c("chromosome", "gene", "transcript", 
                                      "description", "gene_biotype", 
                                      "transcript_biotype", "gene_symbol"), 
                                    function (j) {
                                        # search for and add fasta description information
                                        regex_pattern <- ifelse(j != 'description',
                                                                paste0(" ", j, "\\:[^ ]*"),
                                                                paste0(" ", j, "\\:(.)*"))
                                        x <- gsub(paste0(" ", j, "\\:"), "", 
                                                  str_extract(annot, regex_pattern))
                                        x <- setDT(data.frame(x, stringsAsFactors = F))
                                        names(x) <- j
                                        x
                                    }))
        raw_fasta <- cbind(raw_fasta, dt)
        rm(dt, annot)
    }
    rm(nm)
    print(paste0("Done parsing FASTA database ", fasta_file))
    return(setDT(raw_fasta))
}

get_enzyme_vals <- function(enzyme_val) {
    # split values in the cleaver column of the enzyme_df into two values to be
    # used by the cleaver::cleave function for digesting proteins
    if (nchar(enzyme_val) <= 3 | grepl(";|\\|", enzyme_val) | grepl("[[:upper:]]", str_sub(enzyme_val, 1, 1))) {
        if (grepl(";", enzyme_val)) {
            custom_cleave <- str_split_fixed(enzyme_val, ";", 2)[1,]
        } else {
            custom_cleave <- enzyme_val
        }
        enzyme <- "trypsin"
    } else {
        custom_cleave <- NULL
        enzyme <- enzyme_val
    }
    return(list(enzyme = enzyme, custom_cleave = custom_cleave))
}

cleave_single_protein <- function(protein_info, enzyme_info, missed_cleavages,
                                  length_min, length_max) {
    enzyme_vals <- get_enzyme_vals(enzyme_val = enzyme_info$cleaver)
    peptides <- data.frame(
        # cleave the protein to get all peptide sequences
        peptide_seq = unlist(cleave(protein_info$aa_sequence, enzym = enzyme_vals$enzyme,
                                    missedCleavages = c(0: missed_cleavages),
                                    unique = F, custom = enzyme_vals$custom_cleave)),
        # get the start/end positions of the peptide sequences produced by cleavage
        data.frame(cleavageRanges(protein_info$aa_sequence, enzym = enzyme_vals$enzyme,
                                  missedCleavages = c(0: missed_cleavages),
                                  custom = enzyme_vals$custom_cleave)[[1]],
                   stringsAsFactors = F),
        stringsAsFactors = F)
    # make the data.frame a data.table for faster processing
    peptides <- setDT(peptides)
    # filter peptides to retain only those between 7 and 50 a.a. in length
    peptides <- peptides[nchar(peptide_seq) >= length_min & 
                             nchar(peptide_seq) <= length_max]
    # add the pre and post amino acids (just useful to hold onto)
    peptides[, `:=`(protein_id = protein_info$protein_id,
                    pre_aa = str_sub(protein_info$aa_sequence, start - 1, start - 1),
                    post_aa = str_sub(protein_info$aa_sequence, end + 1, end + 1))]
    if (protein_info$protein_start == 1) {
        # adjust pre_aa for protein_start == 1 which means N-terminal Methionine
        peptides[pre_aa == "" & start == 1, pre_aa := "M"]
    }
    # adjust start/end positions to handle N-terminal Methionine cleavage
    peptides[, `:=`(start = start + protein_info$protein_start,
                    end = end + protein_info$protein_start)]
    peptides[, enzyme_name := enzyme_info$enzyme]
    # add missed cleavages
    peptides[, missed_cleavages := internal_cleavage(x = peptide_seq, 
                                                     enzyme_reg_exp = enzyme_info$mzid_regex)]
    return(peptides)
}

internal_cleavage <- function(x, enzyme_reg_exp = "(?<=[KR])(?!P)") {
    return(stri_count_regex(str_sub(x, 1, -2), enzyme_reg_exp))
}

cleave_proteins <- function(fasta_df, enzyme_df, length_min = 7, length_max = 50,
                            missed_cleavages = 2) {
    n_rows <- nrow(fasta_df)
    # cl <- makeCluster(threads - 1)
    # registerDoSNOW(cl)
    # setup a progress bar
    # pb <- txtProgressBar(max = length(n_rows), style = 3)
    # progress <- function(n) setTxtProgressBar(pb, n)
    # opts <- list(progress = progress)
    # n_enzyme <- nrow(enzyme_df)
    peptide_lst1 <- lapply(1:nrow(enzyme_df), function (i) { # loop through trypsin cleavages
        # get enzyme values appropriate for using cleaver::cleave
        print(paste0("Cleaving ", n_rows, " proteins in sequence using ", enzyme_df$enzyme[i],"...", Sys.time()))
        # enzyme_vals <- get_enzyme_vals(enzyme_val = enzyme_df$cleaver[i])
        peptide_lst <- lapply(1:n_rows, function(j) {
            x <- cleave_single_protein(protein_info = as.list(fasta_df[j]),
                                                                  enzyme_info = as.list(enzyme_df[i]),
                                                                  missed_cleavages = missed_cleavages,
                                                                  length_min = length_min,
                                                                  length_max = length_max)
            x
        })        
        print(paste0("Completed cleaving ", n_rows, " proteins using ", enzyme_df$enzyme[i],"...", Sys.time()))
        dt <- rbindlist(peptide_lst, use.names = T, fill = T)
        # add missed cleavages
        # dt[, missed_cleavages := internal_cleavage(x = peptide_seq, enzyme_reg_exp = enzyme_df$mzid_regex[i])]
        rm(peptide_lst)
        dt
    })
    peptide_dt <- rbindlist(peptide_lst1, use.names = T, fill = T)
    rm(peptide_lst1)
    # stop cluster
    # stopCluster(cl)
    # rm(cl)
    # gc()
    # remove duplicates (protein_id, peptide_seq, start, end, pre_aa, post_aa)
    peptide_dt <- peptide_dt[!duplicated(peptide_dt[, .(protein_id, peptide_seq, start, end,
                                                        pre_aa, post_aa)])]
    return(peptide_dt)
}

check_protein_id <- function(dt, sep = '_') {
    if (any(grepl(sep, dt$protein_id))) {
        new_prot <- str_split_fixed(dt$protein_id, sep, n = 2)
        new_prot <- setDT(as.data.frame(new_prot, stringsAsFactors = F))
        names(new_prot) <- c('protein_id', 'donor_id')
        setnames(dt, old = 'protein_id', new = 'protein_id_original')
        dt <- cbind(dt, new_prot)
    }
    return(dt)
}

main_cleave_fasta <- function(default_val) {
    # run all enzymes
    enzymes_allowed <- c("trypsin", "trypsin-p", "chymotrypsin", "pepsin", "V8-DE", 
                         "glutamyl_endopeptidase", "lys-c", "arg-c")
    # load option values
    opt <- get_args(default_val, enzymes_allowed)
    # check arguments
    opt <- check_args(opt, enzymes_allowed, default_val)
    # load protein sequences from fasta file
    fasta_df <- parse_fasta(fasta_file = opt$fasta_file)
    fasta_df[, protein_start := 0]
    if (opt$Nterm_Met_cleavage) {
        # get protein sequences that start with 'M' and 'cleave' that residue, adding protein to list
        fasta_df_M <- fasta_df[str_sub(aa_sequence, 1, 1) == 'M']
        fasta_df_M[, `:=`(aa_sequence = str_sub(aa_sequence, 2, -1),
                          protein_start = 1)]
        # bind rows of fasta data frames
        fasta_df <- rbindlist(list(fasta_df, fasta_df_M), use.names = T)
        rm(fasta_df_M)
    }
    # load enzymes
    enzyme_df <- fread(opt$enzyme_file, sep = ',', header = T, na.strings = c("NA", ""))
    enzyme_df <- enzyme_df[!is.na(cleaver)]
    # enzyme_df <- enzyme_df[enzyme %in% c('trypsin', 'trypsin-p')]
    enzyme_df <- enzyme_df[enzyme %in% opt$enzymes]
    # cleave proteins by allowed enzymes
    print(paste('Digesting FASTA into peptides...', opt$fasta_file, '...', Sys.time()))
    peptide_dt <- cleave_proteins(fasta_df = fasta_df, enzyme_df = enzyme_df, 
                                  length_min = opt$length_min, length_max = opt$length_max,
                                  missed_cleavages = opt$missed_cleavages)
    # separate protein and donor id values
    peptide_dt <- check_protein_id(dt = peptide_dt, sep = opt$protein_donor_id_sep)
    if (!is.null(opt$donor_id)) peptide_dt[, donor_id := opt$donor_id]
    # write to file
    output_file <- file.path(opt$output_dir, 
                             paste0(file_path_sans_ext(file_path_sans_ext(basename(opt$fasta_file))),
                                    "__", paste(opt$enzymes, collapse = '_'), "__",
                                    opt$length_min, '-', opt$length_max, 
                                    ifelse(opt$Nterm_Met_cleavage, '_NMet', ''),
                                    '.tsv.gz'))
    fwrite(peptide_dt, sep = '\t', na = "", file = output_file)
    print(paste('Completed digesting FASTA into peptides...', opt$fasta_file, '...', Sys.time()))
    cat('\n')
    cat(output_file)
}

main_cleave_fasta(default_val)

# RScript --vanilla proteos_cleave_fasta_peptides_DIA_pipeline.R
# -f "L:/Proteos/fasta/human/Homo_sapiens_GRCh38_pep_all_proteomics_target.fasta"
# -c 2 -l 7 -L 30 -M "TRUE"
# -e "trypsin,trypsin-p"
# -E "U:/Projects/proteos_pipelines/config/enzyme_ref.csv"
# -o "Y:/fasta/human/"

# DT <- data.table(fasta_file = list.files('Y:/fasta/HG38Proteomes',
#                                          pattern = '((Pr)|(PR)|(SA)).*\\.fa\\.gz',
#                                          full.names = T))
# DT[, donor_id := gsub('A','K', toupper(str_extract(basename(fasta_file), '(P[Rr])|(SA)')))]
# DT[, donor_id_tmp := as.integer(str_extract(basename(fasta_file), '[0-9]{1,3}'))]
# DT[, donor_id_num := character()]
# DT[donor_id == 'PR', donor_id_num := formatC(donor_id_tmp, width = 2,
#                                              format = 'd', flag = '0')]
# DT[donor_id == 'SK', donor_id_num := formatC(donor_id_tmp, width = 3,
#                                              format = 'd', flag = '0')]
# DT[, donor_id := paste0(donor_id, donor_id_num)]
# DT[, donor_id_tmp := NULL]
# DT[, donor_id_num := NULL]

# for (k in 1:nrow(DT)) {
#     default_val$missed_cleavages <- 2
#     default_val$length_min <- 7
#     default_val$length_max <- 50
#     default_val$Nterm_Met_cleavage <- TRUE
#     default_val$enzymes <- "trypsin;trypsin-p"
#     default_val$fasta_file <- DT$fasta_file[k]
#     default_val$donor_id <- DT$donor_id[k]
#     main_cleave_fasta(default_val)
# }

# default_val$fasta_file <- 'Y:/fasta/human/Homo_sapiens_GRCh38_pep_all_proteomics_target.fasta'
# default_val$output_dir <- "Y:/fasta/human/"
# default_val$donor_id <- 'GRCh38_ref'
# default_val$missed_cleavages <- 2
# default_val$length_min <- 7
# default_val$length_max <- 50
# default_val$Nterm_Met_cleavage <- TRUE
# default_val$enzymes <- "trypsin;trypsin-p"
# main_cleave_fasta(default_val)
