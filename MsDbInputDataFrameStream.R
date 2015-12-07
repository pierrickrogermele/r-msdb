if ( ! exists('MsDbInputDataFrameStream')) { # Do not load again if already loaded

	library(methods)
	source('MsDbInputStream.R')

	#####################
	# CLASS DECLARATION #
	#####################
	
	MsDbInputDataFrameStream <- setRefClass("MsDbInputDataFrameStream", contains = 'MsDbInputStream', fields = list( .df = "ANY", .i = "integer"))
	
	###############
	# CONSTRUCTOR #
	###############
	
	MsDbInputDataFrameStream$methods( initialize = function(df = data.frame(), input.fields = msdb.get.dft.input.fields(), ...) {

		.df <<- df
		.i <<- 0L
		
		callSuper(input.fields = input.fields, ...)
	})

	##########
	# GET MZ #
	##########

	MsDbInputDataFrameStream$methods( getMz = function() {

		if (.self$.i > 0 && .self$.i <= nrow(.self$.df) && ! is.null(.self$.input.fields[[MSDB.TAG.MZ]]))
			return(.self$.df[.self$.i, .self$.input.fields[[MSDB.TAG.MZ]]])

		return(NULL)
	})

	##########
	# GET RT #
	##########

	MsDbInputDataFrameStream$methods( getRt = function() {

		if (.self$.i > 0 && .self$.i <= nrow(.self$.df) && ! is.null(.self$.input.fields[[MSDB.TAG.RT]]))
			return(.self$.df[.self$.i, .self$.input.fields[[MSDB.TAG.RT]]])

		return(NULL)
	})

	###########
	# GET ALL #
	###########

	MsDbInputDataFrameStream$methods( getAll = function(but = NULL) {

		if (.self$.i > 0 && .self$.i <= nrow(.self$.df)) {

			vals <- .self$.df[.self$.i, , drop = FALSE]

			if ( ! is.null(but))
				vals <- vals[, ! colnames(vals) %in% .self$.input.fields[but]]

			return(vals)
		}

		return(NULL)
	})

	###############
	# NEXT VALUES #
	###############

	MsDbInputDataFrameStream$methods( nextValues = function() {

		if (.self$.i <= nrow(.self$.df))
			.i <<- .self$.i + 1L
	})

	###################
	# HAS NEXT VALUES #
	###################

	MsDbInputDataFrameStream$methods( hasNextValues = function() {
		return(.self$.i < nrow(.self$.df))
	})

} # end of load safe guard
