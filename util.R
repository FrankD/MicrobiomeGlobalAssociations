# Utility functions
require(phyloseq)
require(readr)

# Load in phyloseq data from individual samples
load_in_phyloseq <- function(data.dir='~/Work/Projects/IBDProject/Data/') {
  sample.list = list.dirs(paste0(data.dir, 'counts_biom/'), 
                          full.names=FALSE, recursive=FALSE)
  
  # Get mapping from A1 etc. labels to sample IDs
  temp.data = read_csv(file=paste0(data.dir, 'level_2_phylum_mb.csv'))
  mapping = temp.data[1,2:dim(temp.data)[2]]
  control.status = temp.data[2,2:dim(temp.data)[2]]
  names(control.status) = mapping
  
  # Retrieve only samples in mapping
  sample.list = sample.list[sample.list %in% names(mapping)]
  
  # Need to set sample IDs
  samples = lapply(sample.list, function(sample.dir) {
    sample.phyloseq = import_biom(paste0(data.dir, 'counts_biom/', 
                                         sample.dir, '/otu_table.biom'))
    sample_names(sample.phyloseq) = mapping[sample.dir]
    return(sample.phyloseq)
  })
  
  return(list(data=do.call(merge_phyloseq, samples), status=control.status))
  
}

# Place root at longest individual branch
root.by.outgroup <- function(tree.unrooted,max.length=0.25) {
  #longest.edge <- which.max(tree.unrooted$edge.length)
  #long.nodes <- tree.unrooted$edge[longest.edge,]     #this should usually include one tip
  long.nodes <- tree.unrooted$edge[tree.unrooted$edge.length>max.length,]  # hacky way to find a tip since longest edge not a tip
  new.outgroup <- long.nodes[long.nodes < ape::Ntip(tree.unrooted)]
  tree.rooted <- ape::root(tree.unrooted, outgroup=new.outgroup, resolve.root=T)
  tree.rooted
}

# Load in phyloseq data from individual samples
load_in_phyloseq_reanalysis <- function(data.dir='~/Work/Projects/IBDProject/Data/') {
  samples.phyloseq = import_biom(paste0(data.dir, 'otu_table_mc2_w_tax_no_pynast_failures.biom'),
                                treefilename=paste0(data.dir, 'rep_set.tre'))
  
  # Place tree root at longest individual branch
  phy_tree(samples.phyloseq) = root.by.outgroup(phy_tree(samples.phyloseq))
  
  # Get mapping from A1 etc. labels to sample IDs
  temp.data = read_csv(file=paste0(data.dir, 'level_2_phylum_mb.csv'))
  mapping = temp.data[1,2:dim(temp.data)[2]]
  control.status = temp.data[2,2:dim(temp.data)[2]]
  names(control.status) = mapping
  
  # Swap sample names for more informative IDs
  sample.list = sample_names(samples.phyloseq)
  sample_names(samples.phyloseq) = mapping[sample.list]
  
  return(list(data=samples.phyloseq, status=control.status))
  
}