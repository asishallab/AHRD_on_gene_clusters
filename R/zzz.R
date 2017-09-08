# Function ensured data-sets to be loaded, whenever this package is loaded.
.onLoad <- function( libname = find.package( "AhrdOnGeneClusters" ), pkgname = "AhrdOnGeneClusters" ) {
    data( "ipr_and_families", package = "AhrdOnGeneClusters" )
    message("\nMAKE SURE THAT YOUR PROTEIN-INTERPRO-ANNOTATIONS ARE UNIQUE BEFORE YOU USE AhrdOnGeneClusters\n")
}
