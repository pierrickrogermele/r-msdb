test.searchmz.ticket.2016031010000034 <- function() {
	call.search.mz(c('-d', 'file', '--url', 'res/ticket-2016031010000034-database_CEA_test_2_utf8.tsv', '-m', 'pos', '-i', 'res/ticket-2016031010000034-input_file_for_db_test_2.tsv', '--input-col-names', 'mz=mzmed,rt=rtmed', '-o', 'ticket-2016031010000034-output.tsv', '--peak-output-file', 'ticket-2016031010000034-output-peaks.tsv', '--html-output-file', 'ticket-2016031010000034-output.html'), use.global.conn.flags = FALSE)
}

test.searchmz.peakforest.estelle.20170314 <- function() {
	call.search.mz(c('-d', 'peakforest', '--url', 'https://peakforest-alpha.inra.fr/rest', '--db-token', ENV[['RMSBD_PEAKFOREST_TOKEN']], '-m', 'pos', '-p', '5', '-s', '0', '-i', 'res/peakforest.estelle.20170314-input.tsv', '--input-col-names', 'mz=mzmed,rt=rtmed', '-o', 'peakforest.estelle.20170314-output.tsv', '--peak-output-file', 'peakforest.estelle.20170314-output-peaks.tsv', '--html-output-file', 'peakforest.estelle.20170314-output.html'), use.global.conn.flags = FALSE)
}
