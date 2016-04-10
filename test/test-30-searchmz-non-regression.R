test.match.mz.no.y.col.in.output <- function() {
	# Test that no 'y' column is added in output when running an MZ match
	call.search.mz(c('-m pos', '-i mz-input-small.tsv', '-o mz-output.tsv', '--same-rows'))
	df <- read.table('mz-output.tsv', header = TRUE)
	checkTrue( ! 'y' %in% colnames(df))
}

test.searchmz.ticket.2016031010000034 <- function() {
	call.search.mz(c('-d', 'file', '--url', 'res/ticket-2016031010000034-database_CEA_test_2_utf8.tsv', '-m', 'pos', '-i', 'res/ticket-2016031010000034-input_file_for_db_test_2.tsv', '--input-col-names', 'mz=mzmed,rt=rtmed', '-o', 'ticket-2016031010000034-output.tsv', '--peak-output-file', 'ticket-2016031010000034-output-peaks.tsv', '--html-output-file', 'ticket-2016031010000034-output.html'), use.global.conn.flags = FALSE)
}
