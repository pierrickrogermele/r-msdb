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

	###########
	# GET URL #
	###########

	MsPeakForestDb$methods( .get.url = function(url, params = NULL, ret.type = 'json') {

		res <- NULL

		content <- .self$.url.scheduler$getUrl(url = url, params = params)

		if (ret.type == 'json') {

			library(RJSONIO)

			res <- fromJSON(content)

			if (class(res) == 'list' && 'success' %in% names(res) && res$success == FALSE) {
				param.str <- if (is.null(params)) '' else paste('?', vapply(names(params), function(p) paste(p, params[[p]], sep = '='), FUN.VALUE = ''), collapse = '&', sep = '')
				stop(paste0("Failed to run web service. URL was \"", url, param.str, "\"."))
			}
		} else {
			if (ret.type == 'integer') {
				if (grepl('^[0-9]+$', content, perl = TRUE))
					res <- as.integer(content)
				else {
					library(RJSONIO)
					res <- fromJSON(content)
				}
			}
		}

		return(res)
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

		# Set URL
		url <- paste0(.self$.url, 'metadata/lc/list-code-columns')
		params <- NULL
		if ( ! is.null(molid))
			params <- list(molids = paste(molid, collapse = ','))

		# Call webservice
		wscols <- .self$.get.url(url = url, params = params)

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

		rt <- list()

		# Set URL
		url <- paste0(.self$.url, 'spectra/lcms/search')
		params <- NULL
		if ( ! is.null(molid))
			params <- list(molids = paste(molid, collapse = ','))

		# Call webservice
		spectra <- .self$.get.url(url = url, params = params)
		if (class(spectra) == 'list' && length(spectra) > 0) {
			for (s in spectra)
				if (is.na(col) || s$liquidChromatography$columnCode %in% col) {
					ret.time <- (s$RTmin + s$RTmax) / 2
					c <- s$liquidChromatography$columnCode
					if (c %in% names(rt)) {
						if ( ! ret.time %in% rt[[c]])
							rt[[c]] <- c(rt[[c]], ret.time)
					} else
						rt[[c]] <- ret.time
				}
		}

		return(rt)
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
			params <- c(molids = paste(non.na.molid, collapse = ','))

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
		if ( ! is.null(molid) && (length(molid) > 1 || ! is.na(molid)))
			params <- c(params, molids = paste(molid, collapse = ','))

		# Run request
		n <- .self$.get.url(url = url, params = params, ret.type = 'integer')

		return(sum(n))
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
		ids <- vapply(spectra, function(x) as.character(x$id), FUN.VALUE = '')
		results <- data.frame(id = ids, stringsAsFactors = FALSE)
		colnames(results) <- MSDB.TAG.MOLID
		results[[MSDB.TAG.MZTHEO]] <- vapply(spectra, function(x) as.numeric(x$theoricalMass), FUN.VALUE = 1.1)
		results[[MSDB.TAG.COMP]] <- vapply(spectra, function(x) as.character(x$composition), FUN.VALUE = '')
		results[[MSDB.TAG.ATTR]] <- vapply(spectra, function(x) as.character(x$attribution), FUN.VALUE = '')

		# RT search
		if ( ! is.null(rt.low) && ! is.null(rt.high)) {
			# Build URL for rt search
			url <- paste0(.self$.url, 'spectra/lcms/range-rt-min/', rt.low, '/', rt.high)
			params <- c(columns = paste(col, collapse = ',')) # TODO XXX What are the chrom col IDs ?
		}

		return(results)
	})
}
