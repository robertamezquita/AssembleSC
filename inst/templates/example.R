#' Combine values into a vector or list
#' 
#' @description
#' Only use the description tag if you have multiple paragraphs.
#' The map function transform the input, returning a vector the same length
#' as the input. 
#' 
#' * `map()` returns a list or a data frame
#' * `map_lgl()`, `map_int()`, `map_dbl()` and `map_chr()` return 
#'     vectors of the corresponding type (or die trying); 
#' * `map_dfr()` and `map_dfc()` return data frames created by row-binding 
#'    and column-binding respectively. They require dplyr to be installed.
#'
#' @details 
#' In the following, we present the bullets of the list:
#' * Four cats are few animals.
#' * forcats is a package.
#' 
#' @section Tidy data:
#' When applied to a data frame, row names are silently dropped. To preserve,
#' convert to an explicit variable with [tibble::rownames_to_column()].
#'
#' @section Scoped filtering:
#' The three [scoped] variants ([filter_all()], [filter_if()] and
#' [filter_at()]) make it easy to apply a filtering condition to a
#' selection of variables.
#' 
#' @param key The bare (unquoted) name of the column whose values will be used 
#'   as column headings.
#' @inheritParams function_to_inherit_from [use to avoid duplication]
#'
#' @seealso [fct_lump()] to automatically convert the rarest (or most common)
#'   levels to "other".
#' 
#' @seealso
#' * [tibble()] constructs from individual columns.
#' * [enframe()] converts a named vector into a two-column tibble (names and 
#'   values).
#' * [name-repair] documents the details of name repair.
#' 
#' @family single table verbs
#'
#' @return A full sentence.
#' 
#' @examples
#' 1 + 1
#' sin(pi)




#' Internal function example - drop last
#'
#' Drops the last element from a vector.
#'
#' @param x A vector object to be trimmed.
#'
#' @noRd



#' Example of combining documentation for related functions
#'
#' @param x,y numeric vectors.
#' @name arith
NULL
#> NULL

#' @rdname arith
add <- function(x, y) x + y

#' @rdname arith
times <- function(x, y) x * y
