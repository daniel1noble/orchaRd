#' @title submerge
#' @description Merge two model results tables (orchard objects).
#' @param object1  object of class \code{orchard}.
#' @param object2  object of class \code{orchard}.
#' @param ... Other arguments passed to submerge.
#' @param mix If \code{TRUE}, it will add the number to the moderator name.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame.
#' @export
#'
submerge <- function(object1, object2, ..., mix = FALSE){
  orchard_list <- list(object1, object2, ...)

  len <- length(orchard_list)
  # merging tables
  tables <- lapply(orchard_list, function(x) x$mod_table)
  tables <- do.call("rbind", tables)

  # merging data
  ## checking moderator names are the same or not
  datas <- lapply(orchard_list, function(x) x$data)
  datas <- do.call("rbind", datas)

  # renaming
  if(mix == TRUE){
  names <- lapply(orchard_list, function(x) x$data$moderator)
  names <- as.vector(unlist(mapply(function(x, y) paste0(x, y), x = names, y = 1:len)))
  datas$moderator <- factor(names)
  tables$name <- levels(factor(names))
  }

  model_results <- list(mod_table = tables, data = datas)

  class(model_results) <- "orchard"

  return(model_results)

}
