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
		.url.scheduler <<- UrlRequestScheduler$new(n = 3, useragent = 'r-msdb ; pierrick.roger@gmail.com', ssl.verifypeer = FALSE)

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
		if ( ! is.null(molid))
			url <- paste0(url, '?filter=', paste(molid, collapse = ','))

		# Call webservice
		json <- .self$.url.scheduler$getUrl(url = url)
		cols <- fromJSON(json)

		# Get column names
		cols <- vapply(cols, function(c) paste(c$constructor, c$length, c$diameter, c$particuleSize, sep = '-'), FUN.VALUE = '')

		return(cols)
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
			url <- paste0(.self$.url, 'compounds/all/names?filter=', paste(non.na.molid, collapse = ','))

			# Call webservice
			json <- .self$.url.scheduler$getUrl(url = url)
			names[ ! is.na(molid)] <- fromJSON(json)
		}

		return(names)
	})

	################
	# FIND BY NAME #
	################

	MsPeakForestDb$methods( findByName = function(name) {

		if (is.null(name))
			return(list())

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
}
