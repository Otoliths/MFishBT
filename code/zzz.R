library(biscale)

bi_class2 <- function (.data, x, y, style = "quantile", dim = 3, keep_factors = FALSE) 
{
  bi_x = bi_y = NULL
  if (missing(.data)) {
    stop("An object containing data must be specified for the '.data' argument.")
  }
  if (missing(x)) {
    stop("A variable must be given for the 'x' argument.")
  }
  if (missing(y)) {
    stop("A variable must be given for the 'y' argument.")
  }
  if (style %in% c("quantile", "equal", "fisher", "jenks") == 
      FALSE) {
    stop("The allowed styles are 'equal', 'fisher', 'jenks', or 'quantile'.")
  }
  # if (is.numeric(dim) == FALSE) {
  #   stop("The 'dim' argument only accepts the numeric values '2' or '3'.")
  # }
  # if (dim != 2 & dim != 3) {
  #   stop("The 'dim' argument only accepts the numeric values '2' or '3'.")
  # }
  if (is.logical(keep_factors) == FALSE) {
    stop("A logical scalar must be supplied for 'keep_factors'. Please provide either 'TRUE' or 'FALSE'.")
  }
  paramList <- as.list(match.call())
  if (!is.character(paramList$x)) {
    xQ <- rlang::enquo(x)
  }
  else if (is.character(paramList$x)) {
    xQ <- rlang::quo(!!rlang::sym(x))
  }
  xQN <- rlang::quo_name(rlang::enquo(x))
  if (!is.character(paramList$y)) {
    yQ <- rlang::enquo(y)
  }
  else if (is.character(paramList$y)) {
    yQ <- rlang::quo(!!rlang::sym(y))
  }
  yQN <- rlang::quo_name(rlang::enquo(y))
  if (xQN %in% names(.data) == FALSE) {
    stop(glue::glue("The given 'x' variable '{var}' is not found in the given data set.", 
                    var = xQN))
  }
  if (yQN %in% names(.data) == FALSE) {
    stop(glue::glue("The given 'y' variable '{var}' is not found in the given data set.", 
                    var = yQN))
  }
  bins_x <- dplyr::pull(.data, !!xQ)
  bins_y <- dplyr::pull(.data, !!yQ)
  if (style == "quantile") {
    bins_x <- classInt::classIntervals(bins_x, n = dim, style = "quantile")
    bins_y <- classInt::classIntervals(bins_y, n = dim, style = "quantile")
  }
  else if (style == "equal") {
    bins_x <- classInt::classIntervals(bins_x, n = dim, style = "equal")
    bins_y <- classInt::classIntervals(bins_y, n = dim, style = "equal")
  }
  else if (style == "fisher") {
    bins_x <- classInt::classIntervals(bins_x, n = dim, style = "fisher")
    bins_y <- classInt::classIntervals(bins_y, n = dim, style = "fisher")
  }
  else if (style == "jenks") {
    bins_x <- classInt::classIntervals(bins_x, n = dim, style = "jenks")
    bins_y <- classInt::classIntervals(bins_y, n = dim, style = "jenks")
  }
  bins_x <- bins_x$brks
  bins_y <- bins_y$brks
  out <- dplyr::mutate(.data, bi_x = cut(!!xQ, breaks = bins_x, 
                                         include.lowest = TRUE))
  out <- dplyr::mutate(out, bi_y = cut(!!yQ, breaks = bins_y, 
                                       include.lowest = TRUE))
  out <- dplyr::mutate(out, bi_class = paste0(as.numeric(bi_x), 
                                              "-", as.numeric(bi_y)))
  if (keep_factors == FALSE) {
    out <- dplyr::select(out, -c(bi_x, bi_y))
  }
  return(out)
}


bi_legend2 <- function (pal, dim = 3, xlab, ylab, size = 10) 
{
  bi_class = bi_fill = x = y = NULL
  x <- pal
  # if (missing(pal) == TRUE) {
  #   stop("A palette must be specified for the 'pal' argument.")
  # }
  # if ("bi_pal_custom" %in% class(pal) == TRUE) {
  #   if (dim == 2 & length(pal) != 4) {
  #     stop("There is a mismatch between the length of your custom palette object and the given dimensions.")
  #   }
  #   else if (dim == 3 & length(pal) != 9) {
  #     stop("There is a mismatch between the length of your custom palette object and the given dimensions.")
  #   }
  # }
  # else if ("bi_pal_custom" %in% class(pal) == FALSE) {
  #   if (pal %in% c("Brown", "DkBlue", "DkCyan", "DkViolet", 
  #                  "GrPink") == FALSE) {
  #     stop("The given palette is not one of the allowed options for bivariate mapping. Please choose one of: 'Brown', 'DkBlue', 'DkCyan', 'DkViolet', or 'GrPink'.")
  #   }
  # }
  # if (is.numeric(dim) == FALSE) {
  #   stop("The 'dim' argument only accepts the numeric values '2' or '3'.")
  # }
  # if (dim != 2 & dim != 3) {
  #   stop("The 'dim' argument only accepts the numeric values '2' or '3'.")
  # }
  if (missing(xlab) == TRUE) {
    xlab <- "x var "
  }
  if (is.character(xlab) == FALSE) {
    stop("The 'xlab' argument must be a character string.")
  }
  if (missing(ylab) == TRUE) {
    ylab <- "y var "
  }
  if (is.character(ylab) == FALSE) {
    stop("The 'ylab' argument must be a character string.")
  }
  if (is.numeric(size) == FALSE) {
    stop("The 'size' argument must be a numeric value.")
  }
  xQN <- rlang::quo_name(rlang::enquo(xlab))
  yQN <- rlang::quo_name(rlang::enquo(ylab))
  # if ("bi_pal_custom" %in% class(pal) == TRUE) {
  #   x <- pal
  # }
  # else if ("bi_pal_custom" %in% class(pal) == FALSE) {
  #   if (pal == "DkViolet") {
  #     x <- pal_dkviolet(n = dim)
  #   }
  #   else if (pal == "GrPink") {
  #     x <- pal_grpink(n = dim)
  #   }
  #   else if (pal == "DkBlue") {
  #     x <- pal_dkblue(n = dim)
  #   }
  #   else if (pal == "DkCyan") {
  #     x <- pal_dkcyan(n = dim)
  #   }
  #   else if (pal == "Brown") {
  #     x <- pal_brown(n = dim)
  #   }
  # }
  x <- dplyr::tibble(bi_class = names(x), bi_fill = x)
  leg <- tidyr::separate(x, bi_class, into = c("x", "y"), sep = "-")
  leg <- dplyr::mutate(leg, x = as.integer(x), y = as.integer(y))
  legend <- ggplot2::ggplot() + ggplot2::geom_tile(data = leg, 
                                                   mapping = ggplot2::aes(x = x, y = y, fill = bi_fill)) + 
    ggplot2::scale_fill_identity() + ggplot2::labs(x = substitute(paste(xQN, 
                                                                        "" %->% "")), y = substitute(paste(yQN, "" %->% ""))) + 
    #bi_theme() + 
    ggplot2::theme(axis.title = ggplot2::element_text(size = size),
                   plot.subtitle = element_text(face = "bold",size = size),
                   axis.text = element_text(size = size)) + 
    ggplot2::coord_fixed()
  return(legend)
}



