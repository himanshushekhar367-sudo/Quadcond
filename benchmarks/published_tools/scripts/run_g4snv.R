.libPaths(c("/home/claude/Rlib", .libPaths()))
suppressMessages({library(G4SNVHunter); library(Biostrings); library(GenomicRanges)})
P <- read.csv("pairs_g4_tm.csv", stringsAsFactors=FALSE)
P$pair_id <- seq_len(nrow(P)) - 1
flank <- strrep("T", 30)
out <- list()
for (st in seq(1, nrow(P), by = 100)) {
  p <- P[st:min(nrow(P), st + 99), ]
  seqs <- DNAStringSet(paste0(flank, p$ref, flank)); names(seqs) <- paste0("pair", p$pair_id)
  g4 <- G4HunterDetect(seqs, threshold = 1, window_size = 25)
  v <- GRanges(seqnames = names(seqs), ranges = IRanges(start = 30 + p$pos + 1, width = 1),
               REF = p$ref_base, ALT = p$alt_base, pair_id = p$pair_id)
  res <- tryCatch(G4VarImpact(G4 = g4, variants = v, ref_col = "REF", alt_col = "ALT"), error = function(e) NULL)
  if (!is.null(res)) {
    df <- as.data.frame(mcols(res))
    out[[length(out) + 1]] <- data.frame(pair_id = df$pair_id, ref_score = df$max_score, mut_score = df$mutated.max_score, delta = df$score.diff)
  }
  cat(st, "\n")
}
o <- do.call(rbind, out)
o <- aggregate(. ~ pair_id, data = o, FUN = function(x) x[which.max(abs(x))])
write.csv(o, "scores_g4snvhunter.csv", row.names = FALSE)
cat(nrow(o), "pairs scored of", nrow(P), "\n")
