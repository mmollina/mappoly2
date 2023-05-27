get_seq_indices <- function(input.seq)
  match(input.seq$mrk.names, input.seq$data$mrk.names)

embedded_to_numeric <- function(x) {
  as.integer(gsub("[^0-9]", "", x))
}
