test.match.mz.no.y.col.in.output <- function() {
	# Test that no 'y' column is added in output when running an MZ match
	call.search.mz(c('-m pos', '-i mz-input-small.tsv', '-o mz-output.tsv', '--same-rows'))
	df <- read.table('mz-output.tsv', header = TRUE)
	checkTrue( ! 'y' %in% colnames(df))
}
