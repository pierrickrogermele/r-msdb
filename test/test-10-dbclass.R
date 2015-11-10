source(file.path(dirname(script.path), '..', 'MsDbInputDataFrameStream.R'), chdir = TRUE)
source(file.path(dirname(script.path), '..', 'MsDbOutputDataFrameStream.R'), chdir = TRUE)
source(file.path(dirname(script.path), '..', 'common.R'), chdir = TRUE)

test.10.all.mol.ids <- function() {
	molids <- get.db()$getMoleculeIds()
	checkTrue(length(molids) > 0)
	checkTrue( all( ! duplicated(molids)))
}

test.15.nb.mols <- function() {
	checkTrue(get.db()$getNbMolecules() > 0)
	checkTrue(get.db()$getNbMolecules() == length(get.db()$getMoleculeIds()))
}

test.mol.names <- function() {
	molids <- get.db()$getMoleculeIds()
	checkTrue(is.na(get.db()$getMoleculeName(NA_character_)))
	checkTrue(is.na(get.db()$findByName(NA_character_)[[1]]))
	submolids <- molids[1:10] # get subset of molids
	checkTrue(length(get.db()$getMoleculeName(submolids)) == length(submolids))
	submolids <- c(molids[1:10], NA_character_) # get subset of molids
	checkTrue(length(get.db()$getMoleculeName(submolids)) == length(submolids))
	for (molid in molids[1:20]) {
		name <- get.db()$getMoleculeName(molid)
		checkTrue( ( ! is.na(molid) && ! is.na(name)) || (is.na(molid) && is.na(name)) )
		if ( ! is.na(molid))
			checkTrue(molid %in% get.db()$findByName(name)[[1]])
	}
}

test.unknown.mol <- function() {
	molids <- get.db()$getMoleculeIds()
	for (molid in 1:20)
		if ( ! molid %in% molids)
			checkTrue(is.na(get.db()$getMoleculeName(molid)))
}

test.findbyname <- function() {
	checkTrue(is.na(get.db()$findByName(NULL)))

	molids <- get.db()$getMoleculeIds()
	molids <- molids[1:10]
	names <- get.db()$getMoleculeName(molids)
	ids.of.names <- get.db()$findByName(names)
	checkTrue(length(molids) == length(ids.of.names))
	for (i in seq(molids))
		checkTrue(molids[[i]] %in% ids.of.names[[i]])
}

test.columns <- function() {
	checkTrue(class(get.db()$getChromCol()) == 'character')
	checkTrue(length(get.db()$getChromCol()) > 0)
	molids <- get.db()$getMoleculeIds()
	checkTrue(length(get.db()$getChromCol(molids[1:10])) > 0)
}

test.peaks <- function() {
	checkTrue(get.db()$getNbPeaks() > 0)
	checkTrue(get.db()$getNbPeaks(type = MSDB.TAG.POS) > 0)
	checkTrue(get.db()$getNbPeaks(type = MSDB.TAG.NEG) > 0)
	molids <- get.db()$getMoleculeIds()
	submolids <- molids[1:10] # get subset of molids
	checkTrue(get.db()$getNbPeaks(submolids) > 0)
	checkTrue(get.db()$getNbPeaks(submolids, MSDB.TAG.POS) > 0)
	checkTrue(get.db()$getNbPeaks(submolids, MSDB.TAG.NEG) > 0)
}

test.rt <- function() {
	checkTrue( ! is.null(get.db()$getChromCol()))
	checkTrue(length(get.db()$getChromCol()) > 0)

	# Get all molecule ids
	molids <- get.db()$getMoleculeIds()
	badids <- (1:10000)[! 1:10000 %in% molids]

	for (badid in badids[1:10]) {
		if ( ! is.na(badid)) {
			checkTrue(length(get.db()$getRetentionTimes(badid)) == 0)
			checkTrue(length(get.db()$getChromCol(badid)) == 0)
		}
	}

	for (molid in molids[1:10]) {
		if ( ! is.na(molid)) {
			checkTrue( ! is.null(get.db()$getRetentionTimes(molid)))
			checkTrue( ! is.null(get.db()$getChromCol(molid)))
		}
	}
}

test.search.mz.no.result <- function() {
	mzvals <- get.db()$getMzValues(mode = MSDB.TAG.POS)
	mzvals <- sort(mzvals)
	r <- get.db()$searchForMzRtList(x = msdb.make.input.df(mz = mzvals[1] - 10), mode = MSDB.TAG.POS, prec = 5)
	checkTrue(all(c('id', 'mz') %in% colnames(r)))
	checkTrue(nrow(r) == 1)
	checkTrue(is.na(r[1, 'id']))
}

test.search.mz <- function() {

	mzvals <- get.db()$getMzValues(mode = MSDB.TAG.POS)
	mzvals <- sort(mzvals)

	for (mz in mzvals) {
		if ( ! is.na(mz)) {

			# Search
			r <- get.db()$searchForMzRtList(x = msdb.make.input.df(mz = mz), mode = MSDB.TAG.POS, prec = 5)
			checkTrue(nrow(r) >= 1)
			checkTrue( ! is.na(r[1, 'id']))

			checkTrue(nrow(get.db()$searchForMzRtList(msdb.make.input.df(mz = mz), mode = MSDB.TAG.POS, prec = 5)) >= 1)
			checkTrue(nrow(get.db()$searchForMzRtList(msdb.make.input.df(mz = mz), mode = MSDB.TAG.POS, prec = 100)) >= 1)
			checkTrue(nrow(get.db()$searchForMzRtList(msdb.make.input.df(mz = c(mz)), mode = MSDB.TAG.POS, prec = 100)) >= 1)

			break
		}
	}
}

test.search.mzrt <- function() {
	mzvals <- get.db()$getMzValues(mode = MSDB.TAG.POS)
	mzvals <- sort(mzvals)

	# Search with no result
	r <- get.db()$searchForMzRtList(x = msdb.make.input.df(mz = mzvals[1] - 10), mode = MSDB.TAG.POS, prec = 5)
	checkTrue(all(c('id', 'mz') %in% colnames(r)))
	checkTrue(nrow(r) == 1)
	checkTrue(is.na(r[1, 'id']))

	for (mz in mzvals) {
		if ( ! is.na(mz)) {

			# Search with mz only
			r <- get.db()$searchForMzRtList(x = msdb.make.input.df(mz = mz), mode = MSDB.TAG.POS, prec = 5)
			checkTrue(nrow(r) >= 1)
			checkTrue( ! is.na(r[1, 'id']))
			molid <- r[1, 'id']

			# Get retention times of molecule
			rts <- get.db()$getRetentionTimes(molid)

			# Loop on all columns
			for (col in names(rts))
				for (rt in rts[col])
					checkTrue(nrow(get.db()$searchForMzRtList(x = msdb.make.input.df(mz = mz, rt = rt), mode = MSDB.TAG.POS, prec = 5, col = col)) >= 1)

			break
		}
	}
}
