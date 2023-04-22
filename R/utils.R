#'@export
get_seq_indices <- function(input.seq)
  match(input.seq$mrk.names, input.seq$data$mrk.names)
