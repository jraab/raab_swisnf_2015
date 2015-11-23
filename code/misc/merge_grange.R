# Function to merge two genmoic ranges. 
# takes two genomic ranges and a column to merge on 
# drops strand information
# returns a file of peak_name, and any annottation info merged 
# ----------------------------------------------------------
merge_grange <- function(gr1, gr2, cols_to_keep) { 
    m <- merge(as.data.frame(gr1), as.data.frame(gr2), by = 'name') 
    m <- m[,cols_to_keep]
    return(m) 
} 
