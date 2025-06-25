# pwSEM 1.0.0

I have corrected the following problems:
1.  I have added a reference for this package.  However,
it is for the 3rd edition of a book and the book has not
yet been published.

2.  I have writen TRUE and FALSE instead of T and F everywhere.

3. I have a \value-tag for every .Rd-file containing info about the structure of the output (class) and also what the output means.
I have done this using #' @returns in roxygen2

4.  I have added the following lines in the View.paths function
oldpar <- par(no.readonly = TRUE) 
on.exit(par(oldpar))
