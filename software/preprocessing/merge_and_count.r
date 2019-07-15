###function for merging peaksets and counting the reads in them
#projectDir is where the project/data is
#name is the joint name of the folder
#control is the name of the control 
#toCompare is the names of things to compare to the control
#overwrite dictates whether to overwrite old files or read them in if they exist


merge_and_count = function(projectDir, name, control, toCompare) {
    
    top_dir = getwd()
    #this code itself is happiest in the 'out' folder of the project
    setwd(paste0(projectDir, "/out"))
    
    #This pulls all the technical replicates for a the inputs 
    #Assuming they have the format {name}_x
    samples = list.files(pattern = paste0(control,'_'))
    names(samples) = rep(control, length(samples))
    for(n in toCompare) {
        samples = c(samples, list.files(pattern = paste0('^', n, '_')))
        names(samples)[names(samples) == ''] = n
    }
    
    #run code to merge peaksets
    command = paste0(top_dir, '/software/peak_merge.sh ', name, ' ', paste(samples, collapse=' '))
    print(command)
    system(command)

    #read peak union
    peaks = read.table(paste0(name, '/', name, '_union.bed'))
    peaks = cbind(paste0('peak', 1:nrow(peaks)), peaks, rep('-', nrow(peaks)))
    names(peaks) = c("GeneID", 'Chr', 'Start', 'End', 'Strand')
    
    #make feature counts for each sample
    reads = c()
    for(b in samples) {
        reads = cbind(reads, Rsubread::featureCounts(
            paste0(b,'/',b,'_sorted.bam'), annot.ext=peaks)$counts)
    }
    colnames(reads) = samples
    
    write.csv(reads, paste0(name, '/read_table.csv'), quote=F)
    
    setwd(top_dir)
    return(list(reads=reads, peaks=peaks))
}
