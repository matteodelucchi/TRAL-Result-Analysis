#' Get the Protein ID from the Number of TR
#'
#' @param tr_dataframe a data.frame with TR results imported through \code{\link{function}}
#' @param no_TR number of TRs per protein.
#'
#' @return the Protein ID toghether with the number of detected TRs
#' @export
protein_id_by_number_of_TR <- function(tr_dataframe, no_TR){
  conttab <- table(tr_dataframe$ID)
  # get ID of the one with == no_TR
  for(id in 1:nrow(conttab)){
    if (conttab[[id]]== no_TR){
      return(as.data.frame(conttab[id], col.names = c("protein id", "Number of TR")))
    }
  }
}
