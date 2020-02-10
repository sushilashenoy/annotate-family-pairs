#!/usr/bin/env Rscript

start <- proc.time()
shh <- suppressPackageStartupMessages

if ( ! shh(require('getopt')) ) {
  stop('This script requires the R package getopt')
}
options(stringsAsFactors=FALSE)

# ======== Get options ========
spec <- matrix(c(
  'help',     'h', 0, 'logical',   'print usage information',
  'input',    'i', 1, 'character', 'input file (.fam) or prefix',
  'output',   'o', 2, 'character', 'output file (default = ./basename(input).pairs)',
  'nthreads', 'n', 2, 'integer',   'number of threads to use for calculation',
  'multi', 'm', 0, 'logical', 'Speeds up computation when the input .fam contains multiple families with the same family structure (e.g. output from ped-sim) to speed up calculations',
  'extra', 'x', 0, 'logical', 'Include extra columns (d1/d2/a) in output',
  'redo-fam', 'r', 0, 'logical', 'Ignore the existing fam IDs & redefine families based on individuals connected via paternal or maternal ID',
  'parents', 'p', 0, 'logical', 'Include individuals who only appear as mothers or fathers in the .fam file'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)
if ( interactive() ) opt <- list(input = 'err.fam', parents = TRUE)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}



if ( is.null(opt$input) ) {
  cat(getopt(spec, usage = TRUE))
  stop('input file is required.')
} else {
  if ( endsWith(opt$input, '.fam') ) {
    fam.file <- opt$input
    file.prefix <- substring(fam.file, 1, nchar(fam.file)-4)
  } else {
    file.prefix <- opt$input
    fam.file <- paste0(file.prefix, '.fam')
  }
}
message(' - Input file: ', fam.file)

if ( file.access(fam.file) == -1 ) {
  stop('File ', fam.file, ' does not exist.')
}

if ( is.null(opt$output) ) {
  output.file <- paste0(basename(file.prefix), '.pairs')
} else {
  output.file <- paste0(opt$output)
}
message(' - Output file: ', output.file)

if ( is.null(opt$nthreads) ) {
  use.cores <- 1
} else {
  use.cores <- opt$nthreads
}
if ( use.cores > 1 ) {
  message(' - Using ', use.cores, ' threads for parallel calculations')
}

# ======= Other options ========
collapse.types <- FALSE
if ( !is.null(opt$multi) ) {
  collapse.types <- TRUE
  message(' - Assuming multiple families with identical structure (--multi)')
}

add.parents <- FALSE
if ( !is.null(opt$parents) ) {
  add.parents <- TRUE
  message(' - Including individuals who only appear as parents (--parents)')
}

add.dda <- FALSE
if ( !is.null(opt$extra) ) {
  add.dda <- TRUE
  message(' - Including additional columns (d1/d2/a) in output (--extra)')
}

construct.fam <- FALSE
if ( !is.null(opt$`redo-fam`) ) {
  construct.fam <- TRUE
  message(' - Ignoring existing family IDs (--redo-fam)')
}

message('Loading required R packages...')

if ( ! shh(require('parallel')) ) {
  use.cores <- 1
  mclapply <- function(mc.cores, ...) lapply(...)
}

# =========== Functions for identifying relationships from pedigree =========
shh(library('igraph'))

get.rel.from.fam <- function (fam) {
  # To find common ancestors build a graph with edges from each individual to their parents
  p.edges <- as.character(as.vector(t(cbind(fam$IID, fam$PID)[fam$PID != '0', ])))
  m.edges <- as.character(as.vector(t(cbind(fam$IID, fam$MID)[fam$MID != '0', ])))
  ig <- make_directed_graph(c(p.edges, m.edges))
  
  # Measure distance from every individual to their ancestors
  ddi.graph <- distances(ig, mode="in")
  ddi <- ddi.graph[match(fam$IID, rownames(ddi.graph)), ][, match(fam$IID, colnames(ddi.graph))]
  rownames(ddi) <- colnames(ddi) <- fam$IID
  ddi[is.na(ddi)] <- Inf
  
  get.dda <- function(i, j) {
    # If i or j are vectors (length > 1 ), vectorize
    if ( length(i) > 1 | length(j) > 1 ) {
      ij.df <- data.frame(i=i, j=j)
      return ( do.call(rbind, lapply(seq_len(nrow(ij.df)), function (k) get.dda(ij.df$i[k], ij.df$j[k]))) )
    }
    # If we get to this point we have one i and one j
    cas <- ddi[is.finite(ddi[, i]) & is.finite(ddi[, j]), c(i, j), drop=FALSE]
    if ( nrow(cas) == 0 ) {
      data.frame(d1=NA, d2=NA, a=0)
    } else {
      # Select only most recent common ancestors
      mrcas <- cas[cas[, 1]==min(cas[, 1]) & cas[, 2]==min(cas[, 2]), , drop=FALSE]
      data.frame(d1=mrcas[1, 1], d2=mrcas[1, 2], a=nrow(mrcas))
    }
  }
  
  get.rel <- function(i, j) {
    ij.df <- data.frame(i=i, j=j)
    dda.df <- get.dda(ij.df$i, ij.df$j)
    return ( get.relationship(dda.df$d1, dda.df$d2, dda.df$a) )
  }
  
  return( get.rel )
}


# ordinal number words
ord <- function(n) {
  ordinals <- c('first', 'second', 'third', 'fourth', 'fifth')
  if ( n %in% 1:5 ) {
    ordinals[n]
  } else {
    paste0(n, 'th')
  } 
}

# multiple adverb words
madv <- function(n) {
  ordinals <- c('once', 'twice', 'thrice')
  if ( n %in% 1:3 ) {
    ordinals[n]
  } else {
    paste0(n, 'x')
  } 
}

# description of relationship from
# meiotic distances from common ancestor(s) (d1, d2)
# and number of common ancestors (a)
get.relationship <- function(d1, d2, a) {
  # If any of d1, d2, or a are vectors (length > 1 ), vectorize
  if ( length(d1) > 1 | length(d2) > 1 | length(a) > 1 ) {
    dda.df <- data.frame(d1, d2, a)
    return ( do.call(rbind, lapply(seq_len(nrow(dda.df)),
                                   function (i) get.relationship(dda.df$d1[i],
                                                                 dda.df$d2[i],
                                                                 dda.df$a[i]))) )
  }
  
  # If we get to this point we have one value for each of d1, d2, a
  if ( a == 0 ) {
    # if a == 0 individuals are not related
    degree <- Inf
    rel <- 'unrelated'
    short <- 'UN'
  } else {
    
    degree <- d1 + d2 - log2(a)
    
    if ( degree < 0 | a > 2^min(d1, d2) ) {
      stop('Impossible relationship: ', d1, ' ', d2, ' ', a)
    }
    
    # Ensure dA >= dB
    if ( d1 >= d2 ) {
      dA <- d1
      dB <- d2
    } else {
      dA <- d2
      dB <- d1
    }
    
    if ( degree == 0 ) {
      # if degree == 0, MZ twins
      rel <- 'MZ_twin'
      short <- 'MZ'
    } else if ( dA == 1 & dB == 1 ) {
      rel <- 'sibling'
      short <- 'S'
    } else if ( dB == 0 | dB == 1 ) {
      if ( dB == 0 ) {
        rel <- 'parent/offspring'
        short <- 'PO'
        ngreat <- degree - 1
      } else {
        rel <- 'avuncular'
        short <- 'AV'
        ngreat <- degree - 2
        # half avuncular!
        if ( a == 1 ) {
          ngreat <- degree - 3
        }
      }
      if ( ngreat > 0 ) {
        rel <- paste0('grand', rel)
        short <- paste0('G', short)
        rel <- paste0(paste0(rep('great-', max(0, ngreat-1)), collapse=''), rel)
        short <- paste0(paste0(rep('G', max(0, ngreat-1)), collapse=''), short)
      }
      
    } else {
      rel <- paste(ord(dB-1), 'cousins', sep='_')
      short <- paste0(dB-1, 'C')
      
      if ( dA != dB ) {
        rel <- paste(rel, madv(dA-dB), 'removed', sep='_')
        short <- paste0(short, dA-dB, 'R')
      }
    }
    
    # If this isn't a parent-child relationship (in which case a = 1), then a=1
    # indicates a half prefix should be added
    if ( dB > 0  ) {
      if ( a == 1 ) {
        rel <- paste('half', rel, sep='_')
        short <- paste0('H', short)
      } else if ( a == 4 ) {
        rel <- paste('double', rel, sep='_')
        short <- paste0('D', short)
      } else if ( a == 2 ) {
        if ( dB != 1 | ( dB == 1 & dA == 1 ) ) {
          # For siblings and cousins also add a prefix for full
          # in cousins this prefix only appears in the short name
          if ( dA == 1)
            rel <- paste('full', rel, sep='_')
          short <- paste0('F', short)
        }
      } else {
        rel <- paste('??', rel, sep='_')
        short <- paste0('??', short)
      }
    }
  }
  
  data.frame(d1=d1, d2=d2, a=a, degree=degree, relationship=rel, short=short)
}


# =========== Load fam data =========
message('Loading fam data...')
fam <- read.table(fam.file, col.names=c('FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'))

# Force IDs to character type to avoid incorrect interpretation of numeric IDs as row indices
fam$FID <- as.character(fam$FID)
fam$IID <- as.character(fam$IID)
fam$PID <- as.character(fam$PID)
fam$MID <- as.character(fam$MID)

if ( construct.fam ) {
  message('Identifying families...')
  p.edges <- as.vector(t(cbind(fam$IID, fam$PID)[fam$PID != '0', ]))
  m.edges <- as.vector(t(cbind(fam$IID, fam$MID)[fam$MID != '0', ]))
  ig <- make_directed_graph(c(p.edges, m.edges))
  cfam <- components(ig)
  fam$FID <- cfam$membership[fam$IID]
  fam$FID[is.na(fam$FID)] <- 1:sum(is.na(fam$FID))+sum(!is.na(fam$FID))
}

# helper function to enumerate pairs of elements in x
make.pairs <- function(x) {
  len <- length(x)
  id <- diag(len) # matrix of size len) x len
  row.idx <- row(id)[upper.tri(id)]
  col.idx <- col(id)[upper.tri(id)]
  cbind(x[row.idx], x[col.idx])
}


if ( collapse.types ) {
  fam$type <- gsub('[1234567890]+$', '', fam$FID)
} else {
  fam$type <- fam$FID
}

families <- unique(fam[, c('FID', 'type')])
family.types <- unique(fam$type)
family.size.by.type <- lapply(family.types, function (x) unique(table(fam$FID[fam$type==x])))
# Confirm that families of the same type have the same size
if ( any(size.not.consistent <- sapply(family.size.by.type, length) > 1) ) {
  stop('Family type(s) ', family.types[size.not.consistent], ' have different sizes. (run without -m/--multi ?)')
}

message('Found ', nrow(families), ' families, ',
        length(family.types), ' family type(s).')

# enumerate pairs of ids for one family type
make.family.pairs <- function(ft, fam) {
  famt <- fam[fam$type==ft, ]
  fid1 <- famt$FID[1]
  fam1 <- famt[famt$FID==fid1, c('IID', 'PID', 'MID')]
  
  if ( collapse.types ) {
    # Remove family name from ID to generalize for all families of this type
    fam1$IID <- substring(fam1$IID, nchar(fid1)+2)
    fam1$PID[fam1$PID != 0] <- substring(fam1$PID[fam1$PID != 0], nchar(fid1)+2)
    fam1$MID[fam1$MID != 0] <- substring(fam1$MID[fam1$MID != 0], nchar(fid1)+2)
  }

  parent.ids <- setdiff(union(fam1$PID, fam1$MID), '0')
  missing.par <- !parent.ids %in% fam1$IID
  if ( any(missing.par) ) {
    fam1p <- rbind(fam1, data.frame(IID=parent.ids[missing.par], PID=0, MID=0))
  } else {
    fam1p <- fam1
  }

  if ( nrow(fam1) < 2 ) {
    return ( NULL )
  }
  
  # Enumerate relationships to appear in output
  if ( add.parents ) {
    pairs <- make.pairs(fam1p$IID)
  } else {
    pairs <- make.pairs(fam1$IID)
  }
  
  # Always include "missing parents" for this function to avoid errors
  get.rel <- get.rel.from.fam(fam1p)
  rels <- get.rel(pairs[, 1], pairs[, 2])
  
  message('Family ', ft, ' has ', nrow(fam1), ' individuals (', nrow(pairs), ' pairs).')
  
  pairs.df <- data.frame(fam.id=ft, id1=pairs[, 1], id2=pairs[, 2],
                         degree=rels$degree, relationship=rels$relationship, shortrel=rels$short)
  if ( add.dda ) {
    cbind(pairs.df, data.frame(d1=rels$d1, d2=rels$d2, a=rels$a))
  } else {
    pairs.df
  }
  
}

# enumerate pairs of ids for all family types
fpairs <- mclapply(family.types, make.family.pairs, fam, mc.cores=use.cores)



fam.pair.size <- unique(range(sapply(fpairs, nrow)))

if ( collapse.types ) {
  # Make a combined list of all pairs for all familes
  all.pairs <- do.call(rbind, mclapply(seq_len(nrow(families)), function (i) {
    fid <- families$FID[i]
    fp <- fpairs[[match(families$type[i], family.types)]]
    fp$fam.id <- fid
    fp$id1 <- paste0(fid, '_', fp$id1)
    fp$id2 <- paste0(fid, '_', fp$id2)
    fp
  }, mc.cores=use.cores))
} else {
  all.pairs <- do.call(rbind, fpairs)
}
message('Identified ', nrow(all.pairs), ' within-family pairs (',
        ifelse(length(fam.pair.size)==1, fam.pair.size,
               paste0(fam.pair.size[1], ' to ', fam.pair.size[2])),
        ' pairs per family).')


message('Writing family pair data to ', output.file)
write.table(all.pairs, output.file, sep='\t', quote=FALSE, row.names=FALSE)

# message('Summary: pairs by degree')
# print(table(all.pairs$degree))
# message('\n')
# 
# message('Summary: pairs by relationship')
# print(table(all.pairs$shortrel))
# message('\n')

elapsed.time <- proc.time()[3] - start[3]
message('Done! (run time ', elapsed.time, ' seconds)')



