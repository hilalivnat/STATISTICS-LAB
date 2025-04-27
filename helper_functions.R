bin_reads = function(input, size){
  output = numeric(ceiling(max(input)/size) )
  dat = ceiling(input / size)
  for (i in 1:length(dat)){
    output[dat[i]] = output[dat[i]] + 1
  }
  output
}

getReadLineFast <- function(locations, beg_region, end_region) {
  # Assumes 'locations' is a sorted numeric vector
  N <- end_region - beg_region + 1 # Size of the output vector
  
  # Filter locations efficiently
  in_region_indices <- which(locations >= beg_region & locations <= end_region)
  if(length(in_region_indices) == 0) {
    return(numeric(N)) # Return zeros if no reads in region
  }
  locs_in_region <- locations[in_region_indices]
  
  # Adjust to 1-based indices relative to beg_regio
  adjusted_locs <- locs_in_region - beg_region + 1
  
  # Use tabulate for fast counting
  read_line <- tabulate(adjusted_locs, nbins = N)
  
  return(read_line)
}