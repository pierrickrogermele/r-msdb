if ( ! exists('MsPeakForestDb')) { # Do not load again if already loaded

	library(methods)
	source('MsDb.R')
	source(file.path('..', 'r-lib', 'UrlRequestScheduler.R'))

	#####################
	# CLASS DECLARATION #
	#####################
	
	MsPeakForestDb <- setRefClass("MsPeakForestDb", contains = "MsDb", fields = list(.url = "character", .url.scheduler = "ANY"))
	
	###############
	# CONSTRUCTOR #
	###############
	
	MsPeakForestDb$methods( initialize = function(url = NA_character_, useragent = NA_character_, ...) {

		# Check URL
		if (is.null(url) || is.na(url))
		    stop("No URL defined for new MsPeakForestDb instance.")

		.url <<- url
		.url.scheduler <<- UrlRequestScheduler$new(n = 3, useragent = useragent)
		.self$.url.scheduler$setVerbose(1L)

		callSuper(...)
	})

	####################
	# GET MOLECULE IDS #
	####################
	
	MsPeakForestDb$methods( getMoleculeIds = function() {

		library(RJSONIO)

		json <- .self$.url.scheduler$getUrl(url = paste0(.self$.url, 'compounds/all/ids'))
		ids <- as.character(fromJSON(json))

		return(ids)
	})

	####################
	# GET NB MOLECULES #
	####################
	
	MsPeakForestDb$methods( getNbMolecules = function() {

		n <- .self$.url.scheduler$getUrl(url = paste0(.self$.url, 'compounds/all/count'))

		return(as.integer(n))
	})
	
	###############################
	# GET CHROMATOGRAPHIC COLUMNS #
	###############################
	
	MsPeakForestDb$methods( getChromCol = function(molid = NULL) {

		library(RJSONIO)

		# Set URL
		url <- paste0(.self$.url, 'metadata/lc/list-code-columns')
		params <- NULL
		if ( ! is.null(molid))
			params <- c(filter = paste(molid, collapse = ','))

		# Call webservice
		json <- .self$.url.scheduler$getUrl(url = url, params = params)
		wscols <- fromJSON(json)

		# Build data frame
		cols <- NULL
		for(id in names(wscols))
			cols <- rbind(cols, data.frame(id = id, title = wscols[[id]]$name, stringsAsFactors = FALSE))

		return(cols)
	})
	
	#######################
	# GET RETENTION TIMES #
	#######################
	
	MsPeakForestDb$methods( getRetentionTimes = function(molid, col = NA_character_) {

		if (is.null(molid) || is.na(molid) || length(molid)  != 1)
			stop("The parameter molid must consist only in a single value.")
			
# TODO
		return(list())
	})
	
	#####################
	# GET MOLECULE NAME #
	#####################

	MsPeakForestDb$methods( getMoleculeName = function(molid) {

		library(RJSONIO)

		if (is.null(molid))
			return(NA_character_)

		# Initialize names
		names <- as.character(molid)

		# Get non NA values
		non.na.molid <- molid[ ! is.na(molid)]

		if (length(non.na.molid) > 0) {
			# Set URL
			url <- paste0(.self$.url, 'compounds/all/names')
			params <- c(filter = paste(non.na.molid, collapse = ','))

			# Call webservice
			json <- .self$.url.scheduler$getUrl(url = url, params = params)
			names[ ! is.na(molid)] <- fromJSON(json)
		}

		return(names)
	})

	################
	# FIND BY NAME #
	################

	MsPeakForestDb$methods( findByName = function(name) {

		library(RJSONIO)

		if (is.null(name))
			return(NA_character_)

		ids <- list()

		for (n in name) {

			if (is.na(n))
				ids <- c(ids, NA_character_)

			else {
				url <- paste0(.self$.url, 'search/compounds/name/', curlEscape(n))
				json <- .self$.url.scheduler$getUrl(url = url)
				compounds <- fromJSON(json)$compoundNames
				ids <- c(ids, list(vapply(compounds, function(c) as.character(c$compound$id), FUN.VALUE = '')))
			}
		}

		return(ids)
	})

	#################
	# GET NB PEAKS #
	#################
	
	MsPeakForestDb$methods( getNbPeaks = function(molid = NA_integer_, type = NA_character_) {

		# Build URL
		url <- paste0(.self$.url, 'spectra/lcms/count-peaks')
		params <- NULL
		if ( ! is.na(type))
			params <- c(params, mode = if (type == MSDB.TAG.POS) 'pos' else 'neg')
		if ( ! is.null(molid) || length(molid) > 1 || ! is.na(molid))
			params <- c(params, molids = paste(molid, collapse = ','))

		# Run request
		n <- .self$.url.scheduler$getUrl(url, params = params)

		return(as.integer(n))
	})
	
	#################
	# GET MZ VALUES #
	#################
	
	MsPeakForestDb$methods( getMzValues = function(mode = NULL) {

		library(RJSONIO)

		# Build URL
		url <- paste0(.self$.url, 'spectra/lcms/peaks/list-mz')

		# Query params
		params <- NULL
		if ( ! is.null(mode))
			params <- c(params, mode = if (mode == MSDB.TAG.POS) 'positive' else 'negative')

		# Call service and get JSON
		json <- .self$.url.scheduler$getUrl(url = url, params = params)

		# Convert JSON to R object
		mz <- fromJSON(json)

		return(mz)
	})

	##############################
	# DO SEARCH FOR MZ RT BOUNDS #
	##############################

	MsPeakForestDb$methods( .do.search.for.mz.rt.bounds = function(mode, mz.low, mz.high, rt.low = NULL, rt.high = NULL, col = NULL, attribs = NULL, molids = NULL) {

		library(RJSONIO)

		# Build URL for mz search
		url <- paste0(.self$.url, 'spectra/lcms/peaks/get-range/', mz.low, '/', mz.high)

		# Call service and get JSON
		json <- .self$.url.scheduler$getUrl(url = url)

		# Convert JSON to R object
		spectra <- fromJSON(json, nullValue = NA)

		# Build result data frame
		results <- data.frame(id = vapply(spectra, function(x) as.character(x$id), FUN.VALUE = ''))
		results[[.self$.output.fields$mztheo]] <- vapply(spectra, function(x) as.numeric(x$theoricalMass), FUN.VALUE = 1.1)
		results[[.self$.output.fields$comp]] <- vapply(spectra, function(x) as.character(x$composition), FUN.VALUE = '')
		results[[.self$.output.fields$attr]] <- vapply(spectra, function(x) as.character(x$attribution), FUN.VALUE = '')

		# RT search
		if ( ! is.null(rt.low) && ! is.null(rt.high)) {
			# Build URL for rt search
			url <- paste0(.self$.url, 'spectra/lcms/range-rt-min/', rt.low, '/', rt.high)
			params <- c(columns = paste(col, collapse = ',')) # TODO XXX What are the chrom col IDs ?
		}

		return(results)
	})
}
