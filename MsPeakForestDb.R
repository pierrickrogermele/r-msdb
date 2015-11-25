if ( ! exists('MsPeakForestDb')) { # Do not load again if already loaded

	library(methods)
	library(RJSONIO)
	source('MsDb.R')
	source(file.path('..', 'r-lib', 'UrlRequestScheduler.R'))

	#####################
	# CLASS DECLARATION #
	#####################
	
	MsPeakForestDb <- setRefClass("MsPeakForestDb", contains = "MsDb", fields = list(.url = "character", .url.scheduler = "ANY"))
	
	###############
	# CONSTRUCTOR #
	###############
	
	MsPeakForestDb$methods( initialize = function(url = NA_character_, ...) {

		# Check URL
		if (is.null(url) || is.na(url))
		    stop("No URL defined for new MsPeakForestDb instance.")

		.url <<- url
		.url.scheduler <<- UrlRequestScheduler$new(n = 3, useragent = 'r-msdb ; pierrick.roger@gmail.com')

		callSuper(...)
	})

	####################
	# GET MOLECULE IDS #
	####################
	
	MsPeakForestDb$methods( getMoleculeIds = function() {

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

		# Set URL
		url <- paste0(.self$.url, 'metadata/lc/list-columns')
		params <- NULL
		if ( ! is.null(molid))
			params <- c(filter = paste(molid, collapse = ','))

		# Call webservice
		json <- .self$.url.scheduler$getUrl(url = url, params = params)
		cols <- fromJSON(json)

		# Get column names
		cols <- vapply(cols, function(c) paste(c$constructor, c$length, c$diameter, c$particuleSize, sep = '-'), FUN.VALUE = '')

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
}
