
# SLICER is no longer on CRAN.
# This module will install the package
# if it does not exist and then attach it
# with `box::use()`



# see if `remotes` is installed
tryCatch(invisible(find.package("remotes")),
         # if not, install it
         error = function(cnd) install.packages("remotes"),
         finally = {function(){
           # a dependency `lle` was also kicked off CRAN
           # thus we must manually install as well.
           tryCatch(invisible(find.package("lle")),
                    error = function(cnd) {
                      remotes::install_github("cran/lle")
                    })

           # see if SLICER is installed
           tryCatch(invisible(find.package("SLICER")),
                    error = function(cnd) {
                      #if not install it via `remotes`
                      remotes::install_github("jw156605/SLICER")
                    })
         }}())


#' export all objects from the SLICER package
#' @export
box::use(SLICER[...])
