test.match.mz.no.y.col.in.output <- function() {
	# Test that no 'y' column is added in output when running an MZ match
	call.search.mz(c('-m pos', '-i', 'mz-input-small.tsv', '-o', 'mz-output.tsv', '--same-rows'))
	df <- read.table(get.res.path('mz-output.tsv'), header = TRUE)
	checkTrue( ! 'y' %in% colnames(df))
}

# Test bug replacing values of columns by integers
test.2017.01.26.w4m.sacurine.phenomenal.demo <- function() {

	res.name <- '2017-01-26-w4m-sacurine-phenomenal-demo'
	res.dir <- file.path(dirname(script.path), 'res', res.name)
	main.output <- file.path(dirname(script.path), paste(res.name, 'main.tsv', sep = '-'))
	peak.output <- file.path(dirname(script.path), paste(res.name, 'peak.tsv', sep = '-'))
	html.output <- file.path(dirname(script.path), paste(res.name, 'html', sep = '.'))

	call.search.mz(c('-d', 'file',
					 '--url', file.path(res.dir, 'massbank-neg-ms1-peaks.tsv'),
					 '--db-fields', 'mztheo=peak.mz,chromcolrt=chromcolrt,compoundid=accession,chromcol=chromcol,msmode=msmode,peakcomp=peak.formula,fullnames=name,compoundmass=mass,compoundcomp=formula,inchi=inchi,inchikey=inchikey,pubchemid=pubchemcompid,chebiid=chebiid,hmdbid=hmdbid,keggid=keggid',
					 '--db-ms-modes', 'pos=pos,neg=neg',
					 '-i', file.path(res.dir, 'Biosigner_variableMetadata.tsv'),
					 '--input-col-names', 'mz=mass_to_charge,rt=retention_time',
					 '-m', 'neg',
					 '-p', '10.0',
					 '-s', '0.0',
					 '--same-rows',
					 '--same-cols',
					 '-o', main.output,
					 '--peak-output-file', peak.output,
					 '--html-output-file', html.output
					 ), use.global.conn.flags = FALSE)
}
