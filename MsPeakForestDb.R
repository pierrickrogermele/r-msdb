if ( ! exists('MsPeakForestDb')) { # Do not load again if already loaded

	library(methods)
	source('MsDb.R')
	source(file.path('..', 'r-lib', 'UrlRequestScheduler.R'))

	#####################
	# CLASS DECLARATION #
	#####################
	
	MsPeakForestDb <- setRefClass("MsPeakForestDb", contains = "MsDb", fields = list(.url = "character", .url.scheduler = "ANY", .token = "character"))
	
	###############
	# CONSTRUCTOR #
	###############
	
	MsPeakForestDb$methods( initialize = function(url = NA_character_, useragent = NA_character_, token = NA_character_, ...) {

		# Check URL
		if (is.null(url) || is.na(url))
		    stop("No URL defined for new MsPeakForestDb instance.")

		.url <<- url
		.url.scheduler <<- UrlRequestScheduler$new(n = 3, useragent = useragent)
		.self$.url.scheduler$setVerbose(1L)
		.token <<- token

		callSuper(...)
	})

	###########
	# GET URL #
	###########

	MsPeakForestDb$methods( .get.url = function(url, params = NULL, ret.type = 'json') {

		res <- NULL

		# Add token
		if ( ! is.na(.self$.token))
			params <- c(params, token = .self$.token)
				param.str <- if (is.null(params)) '' else paste('?', vapply(names(params), function(p) paste(p, params[[p]], sep = '='), FUN.VALUE = ''), collapse = '&', sep = '')
			print(paste0('URL ', url, param.str))

		# Get URL
		content <- .self$.url.scheduler$getUrl(url = url, params = params)

		if (ret.type == 'json') {

			library(RJSONIO)

			res <- fromJSON(content, nullValue = NA)

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
					res <- fromJSON(content, nullValue = NA)
				}
			}
		}

		return(res)
	})

	####################
	# GET MOLECULE IDS #
	####################
	
	MsPeakForestDb$methods( getMoleculeIds = function() {

		ids <- as.character(.self$.get.url(url = paste0(.self$.url, 'compounds/all/ids')))

		return(ids)
	})

	####################
	# GET NB MOLECULES #
	####################
	
	MsPeakForestDb$methods( getNbMolecules = function() {

		n <- .self$.get.url(url = paste0(.self$.url, 'compounds/all/count'), ret.type = 'integer')

		return(n)
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
		cols <- data.frame(id = character(), title = character())
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
			names[ ! is.na(molid)] <- .self$.get.url(url = url, params = params)
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
				compounds <- .self$.get.url(url = url)$compoundNames
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

		# Build URL
		url <- paste0(.self$.url, 'spectra/lcms/peaks/list-mz')

		# Query params
		params <- NULL
		if ( ! is.null(mode))
			params <- c(params, mode = if (mode == MSDB.TAG.POS) 'positive' else 'negative')

		# Get MZ valuels
		mz <- .self$.get.url(url = url, params = params)

		return(mz)
	})

	##############################
	# DO SEARCH FOR MZ RT BOUNDS #
	##############################

	MsPeakForestDb$methods( .do.search.for.mz.rt.bounds = function(mode, mz.low, mz.high, rt.low = NULL, rt.high = NULL, col = NULL, attribs = NULL, molids = NULL) {

		# Build URL for mz search
		url <- paste0(.self$.url, 'spectra/lcms/peaks/get-range/', mz.low, '/', mz.high)

		# Get spectra
		spectra <- .self$.get.url(url = url)

		# Build result data frame
		results <- data.frame(MSDB.TAG.MOLID = character(), MSDB.TAG.MOLNAMES = character(), MSDB.TAG.MZTHEO = numeric(), MSDB.TAG.COMP = character(), MSDB.TAG.ATTR = character())
		for (x in spectra)
			results <- rbind(results, data.frame(MSDB.TAG.MOLID = vapply(x$source$listOfCompounds, function(c) as.character(c$id), FUN.VALUE = ''),
			                                     MSDB.TAG.MOLNAMES = vapply(x$source$listOfCompounds, function(c) paste(c$names, collapse = MSDB.MULTIVAL.FIELD.SEP), FUN.VALUE = ''),
												 MSDB.TAG.MZTHEO = as.numeric(x$theoricalMass),
												 MSDB.TAG.COMP = as.character(x$composition),
												 MSDB.TAG.ATTR = as.character(x$attribution),
												 stringsAsFactors = FALSE))

		# RT search
		if ( ! is.null(rt.low) && ! is.null(rt.high)) {

			rt.res <- data.frame(MSDB.TAG.MOLID = character(), MSDB.TAG.COL = character(), MSDB.TAG.COLRT = numeric())

			if (nrow(results) > 0) {
				# Build URL for rt search
				url <- paste0(.self$.url, 'spectra/lcms/range-rt-min/', rt.low, '/', rt.high)
				params <- NULL
				if ( ! is.null(col))
					params <- c(columns = paste(col, collapse = ','))

				# Run query
				rtspectra <- .self$.get.url(url = url, params = params)

				# Get compound/molecule IDs
				for (x in spectra)
					rt.res <- rbind(rt.res, data.frame(MSDB.TAG.MOLID = vapply(x$listOfCompounds, function(c) as.character(c$id), FUN.VALUE = ''),
				                                   	   MSDB.TAG.COL = as.character(x$liquidChromatography$columnCode),
				                                   	   MSDB.TAG.COLRT = (as.numeric(x$RTmin) + as.numeric(x$RTmax)) / 2,
					                                   	   stringsAsFactors = FALSE))
			}	

			# Add retention times and column info
			results <- merge(results, rt.res)
		}
		
		# Rename columns with proper names
		colnames(results) <- vapply(colnames(results), function(s) eval(parse(text=s)), FUN.VALUE = '')

		return(results)
	})
}
