#' Download genome from NCBI and returns preprocessed object
#'
#' @param accession refseq specific accession identifier
#' @param verbose TRUE/FALSE
#' @export
download_from_ncbi <- function(accession = "GCF_000007825.1", ...){
  require(biomartr)
  # process additional arguments
  dots = list(...)
  file_path <- biomartr::getGenome(db = "refseq",
                      organism = accession,
                      path = file.path("_ncbi_downloads","genomes"))
  # load fasta file
  genome <- biomartr::read_genome(file_path, format = "fasta")
  # get semi-redundant one-hot enconding
  def.vals = list(char = genome, maxlen = 30,
                  vocabulary = c("\n", "a", "c", "g", "t"))
  ind = unlist(lapply(dots[names(def.vals)], is.null))
  dots[names(def.vals)[ind]] = def.vals[ind]
  genome_preprocessed <- do.call(altum::preprocess, dots)
  return(genome_preprocessed)
}