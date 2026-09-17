.libPaths(c("/home/claude/Rlib", .libPaths()))
suppressMessages({library(pqsfinder); library(Biostrings)})
d <- read.csv("allseqs.csv", stringsAsFactors=FALSE)
sc <- numeric(nrow(d)); sc_def <- numeric(nrow(d))
for (i in seq_len(nrow(d))) {
  s <- DNAString(d$sequence[i])
  p <- tryCatch(pqsfinder(s, min_score=10, strand="*", verbose=FALSE), error=function(e) NULL)
  sc[i] <- if (is.null(p) || length(p)==0) 0 else max(score(p))
  if (i %% 2000 == 0) cat(i, "\n")
}
d$pqsfinder <- sc
write.csv(d[,c("id","pqsfinder")], "scores_pqsfinder.csv", row.names=FALSE)
