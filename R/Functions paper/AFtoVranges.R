#' Obtaining a Vranges object from the output of readAF
#'
#' @param AF Output of the function readAF.


AFtoVranges <- function(AF){
	require("VariantAnnotation")
  VRanges(seqnames = AF$chr,
          ranges = as.character(AF$pos),
          ref = AF$ref,
          alt = AF$alt,
          sampleID = AF$sample,
          treatment = AF$treatment,
          species = AF$species,
          af=AF$af,
          SBS=AF$SBS)
}

AFtoGranges <- function(AF){
  require("VariantAnnotation")
  GRanges(seqnames = AF$chr,
          ranges = as.character(AF$pos),
          ref = AF$ref,
          alt = AF$alt,
          sampleID = AF$sample,
          treatment = AF$treatment,
          species = AF$species,
          af=AF$af,
          SBS=AF$SBS)
}
