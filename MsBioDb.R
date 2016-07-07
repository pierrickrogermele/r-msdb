if ( ! exists('MsBioDb')) { # Do not load again if already loaded

	library(methods)
	source('MsDb.R')
	source(file.path('..', 'r-biodb', 'BioDbConn.R'))

	#####################
	# CLASS DECLARATION #
	#####################
	
	MsBioDb <- setRefClass("MsBioDb", contains = "MsDb", fields = list(.biodb = "ANY"))
	
	###############
	# CONSTRUCTOR #
	###############
	
	MsBioDb$methods( initialize = function(biodb = NA_character_, ...) {

		# Check bio database
		! is.na(biodb) || stop("You must set a bio database.")
		inherits(biodb, "BiodbConn") || stop("The bio database must inherit from BiodbConn class.")
		.biodb <<- biodb

		callSuper(...)
	})

	####################
	# GET MOLECULE IDS #
	####################
	
	MsBioDb$methods( getMoleculeIds = function() {
		return(.self$.biodb$getEntryIds(type = RBIODB.COMPOUND))
	})

	####################
	# GET NB MOLECULES #
	####################
	
	MsBioDb$methods( getNbMolecules = function() {
		return(.self$.biodb$getNbEntries(type = RBIODB.COMPOUND))
	})

}
