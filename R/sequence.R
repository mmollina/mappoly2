
make_sequence <- function(input.data,
                          grouping.scope = c("lg.grouping",
                                             "chrom",
                                             "mrk.seq"),
                          genomic.info = NULL) {
  assert_that(has.mappoly2.rf(input.data))
  if (is.mappoly2.data(input.obj))
  {
    pattern <- "(ch|chr|CH|Chr|CHR|chrom|Chrom|Chromsome)"
    ## Sequence with all markers
    if (all(arg  ==  "all"))
    {
      mrk.names <- input.obj$mrk.names
      out.dat <- input.obj
    } ## If chromosome informed
    else if (all(is.character(arg)) &
             sum(grepl(pattern, arg, ignore.case = TRUE))  ==  length(arg) &
             all(!arg%in%rownames(input.obj$geno.dose)))
    {
      if (all(is.na(input.obj$chrom)))
        stop("There is no chromosome information.")
      ch.n.arg <- embedded_to_numeric(arg)
      ch.n.dat <- embedded_to_numeric(input.obj$chrom)
      ch.id <- ch.n.dat%in%ch.n.arg
      mrk.names <- input.obj$mrk.names[ch.id]
      chrom <- input.obj$chrom[mrk.names]
      if (any(!is.na(input.obj$genome.pos)))
        genome.pos <- input.obj$genome.pos[mrk.names]
      ch_geno <- data.frame(chrom, genome.pos)
      sorted_ch_geno <- ch_geno[with(ch_geno, order(chrom, genome.pos)),]
      mrk.names <- rownames(sorted_ch_geno)
      out.dat <- subset_data(input.obj, select.mrk = mrk.names)
    } ## sequence with specific markers
    else if (all(is.character(arg)) & (length(arg)  ==  length(arg %in% input.obj$mrk.names)))
    {
      mrk.names <- intersect(arg, input.obj$mrk.names)
      out.dat <- subset_data(input.obj, select.mrk = mrk.names)
    }
    else if (is.vector(arg) && all(is.numeric(arg)))
    {
      assert_that(max(arg) <= input.obj$n.mrk)
      mrk.names <- input.obj$mrk.names[arg]
      out.dat <- subset_data(input.obj, select.mrk = mrk.names)
    }
    else stop("Invalid argument to select markers")
  }
  else if (is.mappoly2.sequence(input.obj))
  {
    return(make_sequence(input.obj$data, arg, info.parent))
  }
  else if (is.mappoly2.group(input.obj))
  {
    lgs.idx <- names(input.obj$groups.snp[input.obj$groups.snp  %in%  arg])
    if(is.null(genomic.info)){
      return(make_sequence(input.obj = input.obj$input.seq,
                           arg = lgs.idx))
    } else {
      assert_that(is.numeric(genomic.info))
      chrom <- input.obj$input.seq$data$chrom[lgs.idx]
      chrom.table <- sort(table(chrom, useNA = "always"), decreasing = TRUE)
      seq.group <- names(chrom)[chrom %in% names(chrom.table[genomic.info])]
      return(make_sequence(input.obj$input.seq$data, seq.group))
    }
  }
  else if (is.mappoly2.geno.ord(input.obj))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the genome order instead.")
    return(make_sequence(input.obj$data, rownames(input.obj$ord)))
  }
  if (is.mappoly2.pcmap(input.obj) | is.mappoly2.pcmap3d(input.obj))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the MDS order instead.")
    return(input.obj$mds.seq)
  }
  d.p1 <- input.obj$dosage.p1[mrk.names]
  d.p2 <- input.obj$dosage.p2[mrk.names]
  if(info.parent == "p1"){
    mrk.names <- mrk.names[d.p2 == 0 | d.p2 == input.obj$ploidy.p2]
    out.dat <- subset_data(input.obj, select.mrk = mrk.names)
  }
  else if(info.parent == "p2"){
    mrk.names <- mrk.names[d.p1 == 0 | d.p1 == input.obj$ploidy.p1]
    out.dat <- subset_data(input.obj, select.mrk = mrk.names)
  }
  structure(list(mrk.names = mrk.names,
                 phases = phase,
                 pairwise = pairwise,
                 linkage.groups = linkage.groups,
                 order = list(mds = NA, genome = NA),
                 redundant = NULL,
                 data = out.dat),
            class = "mappoly2.sequence")
}

