#' Fast coverage calculation for genomic regions
#'
#' @param locations Sorted numeric vector of read locations
#' @param beg_region Start position of region to analyze
#' @param end_region End position of region to analyze
#' @return Numeric vector of coverage counts
getReadLineFast <- function(locations, beg_region, end_region) {
  # <<< Add input validation (from DeepSeek suggestion) >>>
  if (!is.numeric(beg_region) || !is.numeric(end_region) || length(beg_region) != 1 || length(end_region) != 1) {
      stop("beg_region and end_region must be single numeric values.")
  }
  if (beg_region >= end_region) {
    stop("beg_region must be less than end_region.")
  }
  if (beg_region < 0 || end_region < 0) {
      stop("beg_region and end_region must be non-negative.")
  }
  # <<< End validation add >>>
  
  N <- end_region - beg_region + 1 # Size of the output vector
  
  # Filter locations efficiently using data.table potentially if locations is large
  # Or stick to base R which() which is also efficient
  in_region_indices <- which(locations >= beg_region & locations <= end_region)
  
  # Handle cases with no relevant reads
  if (length(in_region_indices) == 0) {
    # Return vector of zeros of the correct length
    return(numeric(N)) 
  }
  locs_in_region <- locations[in_region_indices]
  
  # Adjust to 1-based indices relative to beg_region
  adjusted_locs <- locs_in_region - beg_region + 1
  
  # Use tabulate for fast counting
  # Ensure adjusted_locs are positive integers; tabulate ignores 0s and non-integers
  # Max value in adjusted_locs should not exceed N theoretically
  if (max(adjusted_locs) > N) {
       warning("Some read locations seem to extend beyond the defined end_region boundary after adjustment.")
       # Decide how to handle: clip or error. Clipping might hide issues.
       # Let's allow tabulate to handle it; it counts based on integer values.
  }
  read_line <- tabulate(adjusted_locs, nbins = N)
  
  return(read_line)
}

# -----------------------------------------------------

#' Efficiently calculate binned coverage across a chromosome
#' (Replaces the need for get_binned_coverage from DeepSeek fix)
#' 
#' @param locations Numeric vector of all read start locations on the chromosome.
#' @param bin_size The size of each bin.
#' @param chromosome_length Estimated or known length of the chromosome.
#' @return Numeric vector of read counts per bin for the entire chromosome.
getBinnedCoverage <- function(locations, bin_size, chromosome_length) {
  if (bin_size < 1) stop("bin_size must be >= 1")
  if (length(locations) == 0) return(numeric(0))
  
  # Define breaks for the entire chromosome
  breaks <- seq(0, chromosome_length + bin_size, by = bin_size)
  num_bins <- length(breaks) - 1
  
  # Use cut to assign each location to a bin index
  # right=FALSE means bins are [start, end), typical for genomics
  # include.lowest=TRUE ensures the very first position is included
  bin_indices <- cut(locations, breaks = breaks, labels = FALSE, include.lowest = TRUE, right = FALSE) 
  
  # Handle potential NAs if locations exceed chromosome_length estimate slightly
  bin_indices_valid <- na.omit(bin_indices) 
  if (length(bin_indices_valid) != length(bin_indices)) {
      warning(sprintf("%d read locations were outside the defined chromosome length range [0, %d) and were ignored during binning.", 
                      length(bin_indices) - length(bin_indices_valid), chromosome_length))
  }
  
  # Use tabulate for efficient counting of bin assignments
  # Need to ensure we cover all possible bins from 1 to num_bins
  # tabulate(factor(...)) would work but table(factor(...)) is often clearer
  bin_counts_table <- table(factor(bin_indices_valid, levels = 1:num_bins))
  
  # Return as numeric vector
  return(as.numeric(bin_counts_table))
}


# -----------------------------------------------------

#' Calculate GC count per bin for a given sequence region
#' (Used for Part 2)
#'
#' @param sequence_vector Character vector representing the DNA sequence.
#' @param seq_start The starting position of the sequence_vector in the chromosome coordinates.
#' @param region_start The starting position of the analysis region.
#' @param region_end The ending position of the analysis region.
#' @param bin_size The size of each bin.
#' @return Numeric vector of GC counts per bin.
calculate_gc_per_bin <- function(sequence_vector, seq_start, region_start, region_end, bin_size) {
  
  # Ensure the requested region is within the provided sequence
  seq_end <- seq_start + length(sequence_vector) - 1
  if (region_start < seq_start || region_end > seq_end ) {
      stop(sprintf("Requested region [%d, %d] is outside the bounds of the provided sequence vector [%d, %d].",
                   region_start, region_end, seq_start, seq_end))
  }
   if (bin_size < 1) stop("bin_size must be >= 1")
   
  # Adjust region coordinates to be relative to the sequence_vector
  relative_start <- region_start - seq_start + 1
  relative_end <- region_end - seq_start + 1
  
  # Calculate the number of bins needed for the specified region
  region_length <- region_end - region_start + 1
  n_bins <- ceiling(region_length / bin_size)
  gc_counts <- numeric(n_bins) # Pre-allocate
  
  cat(sprintf("Calculating GC content for %d bins in region %d-%d...\n", n_bins, region_start, region_end))
  
  for (i in 1:n_bins) {
    # Define bin boundaries within the relative coordinates of the sequence_vector
    bin_relative_start <- relative_start + (i - 1) * bin_size
    bin_relative_end <- min(relative_start + i * bin_size - 1, relative_end)
    
    # Ensure indices are within the bounds of sequence_vector (1 to length(sequence_vector))
    # This also handles the case where the last bin might be smaller
    actual_start_index <- max(1, bin_relative_start)
    actual_end_index <- min(length(sequence_vector), bin_relative_end)

    # If the calculated start is already past the end, skip (can happen with ceiling)
    if (actual_start_index > actual_end_index) {
      gc_counts[i] <- 0 # Assign 0 if bin is empty or outside bounds
      next 
    }
    
    # Extract the sequence segment for the current bin
    segment <- sequence_vector[actual_start_index:actual_end_index]
        
    # Count G and C, case-insensitive and ignoring NAs or other characters
    gc_counts[i] <- sum(segment %in% c("G", "C", "g", "c"), na.rm = TRUE)

  }
  cat(" Done.\n")
  
  return(gc_counts)
}
