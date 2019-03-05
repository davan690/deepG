#' Download genome from NCBI and returns preprocessed object
#'
#' @param accession RefSeq assembly accession e.g. 'GCA_000196515.1' forhttps://www.ncbi.nlm.nih.gov/assembly/GCA_000196515.1 . A collection of RefSeq accessions can be generated
#' using https://www.ncbi.nlm.nih.gov/genome/browse
#' @export
download_from_ncbi <- function(accession = "GCF_000007825.1", ...){
  require(biomartr)
  # process additional arguments
  dots = list(...)
  file_path <- biomartr::getGenome(db = "genbank",
                      organism = accession,
                      reference = FALSE,
                      path = file.path("_ncbi_downloads","genomes"))

  #file_path_gff <- getGFF(db = "genbank",
  #                        organism = accession,
  #                        reference = FALSE,
  #                        path = file.path("_ncbi_downloads", "annotation"))
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

