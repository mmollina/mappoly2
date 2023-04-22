#' Check if an object is of class mappoly2.data
#' @export
is.mappoly2.data <- function (x)
  inherits(x, "mappoly2.data")

#' Check if an object is of class mappoly2.map
#' @export
is.mappoly2.map <- function (x)
  inherits(x, "mappoly2.map")

#' Check if an object is of class mappoly2.sequence
#' @export
is.mappoly2.sequence <- function (x)
  inherits(x, "mappoly2.sequence")

#' Check if an object is of class mappoly2..twopt
#' @export
is.mappoly2.twopt <- function (x)
  inherits(x, "mappoly2.twopt")

#' Check if an object is of class mappoly2.pcmap
#' @export
is.mappoly2.pcmap <- function (x)
  inherits(x, "mappoly2.pcmap")

#' Check if an object is of class mappoly2.pcmap3d
#' @export
is.mappoly2.pcmap3d <- function (x)
  inherits(x, "mappoly2.pcmap3d")

#' Check if an object is of class mappoly2.group
#' @export
is.mappoly2.group <- function (x)
  inherits(x, "mappoly2.group")

#' Check if an object is of class mappoly2.chitest.seq
#' @export
is.mappoly2.chitest.seq <- function (x)
  inherits(x, "mappoly2.chitest.seq")

#' Check if an object is of class mappoly2.geno.ord
#' @export
is.mappoly2.geno.ord <- function (x)
  inherits(x, "mappoly2.geno.ord")
