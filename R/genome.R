load_genome = function(genome_name = NULL, genome_path = NULL) {
  if (!is.null(genome_path) && nzchar(genome_path)) {
    if (!file.exists(genome_path)) {
      stop(sprintf("Genome FASTA does not exist: %s", genome_path))
    }
    cat ("Loading genome from FASTA file: ", genome_path, "\n")
    genome = Biostrings::readDNAStringSet(genome_path)
    names(genome) = sub(" .*", "", names(genome))
    return(genome)
  }

  cat ("Loading genome from BSgenome package: ", genome_name, "\n")
  genome_list = BSgenome::available.genomes()
  if (genome_name %in% genome_list) {
    genome_package = genome_name
  } else {
    if (startsWith(genome_name, "hg")) {
      organism_prefix = "BSgenome.Hsapiens.UCSC."
    } else if (startsWith(genome_name, "mm")) {
      organism_prefix = "BSgenome.Mmusculus.UCSC."
    } else if (startsWith(genome_name, "rn")) {
      organism_prefix = "BSgenome.Rnorvegicus.UCSC."
    } else if (startsWith(genome_name, "dm")) {
      organism_prefix = "BSgenome.Dmelanogaster.UCSC."
    } else if (startsWith(genome_name, "ce")) {
      organism_prefix = "BSgenome.Celegans.UCSC."
    } else if (startsWith(genome_name, "dr")) {
      organism_prefix = "BSgenome.Drerio.UCSC."
    } else {
      stop("Unknown genome name or unsupported organism.")
    }
    genome_package = paste0(organism_prefix, genome_name)
  }

  if (!requireNamespace(genome_package, quietly = TRUE)) {
    BiocManager::install(genome_package)
  }

  BSgenome::getBSgenome(genome_package)
}

get_genome_seq = function(genome, chr, start, end, strand = "+") {
  start = max(1L, as.integer(start))

  if (inherits(genome, "BSgenome")) {
    return(BSgenome::getSeq(genome, chr, start = start, end = end, strand = strand))
  }

  if (!chr %in% names(genome)) {
    stop(sprintf("Chromosome '%s' not found in supplied genome FASTA", chr))
  }

  seq = genome[[chr]]
  end = min(as.integer(end), length(seq))
  out = Biostrings::subseq(seq, start = start, end = end)
  if (strand == "-") {
    out = Biostrings::reverseComplement(out)
  }
  out
}
