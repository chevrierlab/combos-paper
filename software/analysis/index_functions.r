### Some helper functions for isserlis-related computation

### values from key builds out the 1,2,3,12,13,23 and 123 components of a dataframe
### it uses the key to get this information for each row of data, from data
### name indicates which column to pull from in data
### match indicates whether to match on a 'Screen_rep_ID' column as well as Treatment
### if match is false, it matches on just Treatment

values_from_key = function(data, name, key, match) { 
    
    colnames = str_c(name, '_', as.character(c(1,2,3,12,13,23,123)))
    keynames = str_c('C_', c(1,2,3,12,13,23,123))
    
    values = matrix(nrow = nrow(data), ncol = length(colnames), dimnames = list(NULL, colnames))
    
    # for each row in data
    for (i in 1:nrow(data)) {
        # for each component that needs getting
        for (j in 1:length(colnames)) { 
            # get the name from the key
            key_name = unlist(filter(key, Treatment==data$Treatment[i])[keynames[j]])
            if (is.na(key_name)) { 
                # if it the component doesn't exist for this row
                values[i,j] = NA
            } else if (match) {
                # if it does and we're matching
                values[i,j] = 
                    unlist(filter(data, Treatment==key_name, Screen_rep_ID == data$Screen_rep_ID[i])[name]);
            } else {
                # if it does and we're not matching
                values[i,j] = unlist(filter(data, Treatment==key_name)[name])
            }
        }
    }
    
    return(values)
}

### split_values_from_key() works like values_from_key() but has been optimized for bootstrapping
### it takes a different reference data frame (ref_data) from which data is pulled
### so that one can reconstruct only a certain portion of the data (like triplets) 
### it is optimized in various other ways to speed it up, it can take multiple names at once (as a char vector)
### although it is not currently optimized for match=FALSE


split_values_from_key = function(data, ref_data, name, key, match) {
    colnames = unlist(map(name, function(n) str_c(n, '_',  as.character(c(1,2,3,12,13,23,123)))))
    keynames = str_c('C_', c(1,2,3,12,13,23,123))
    
    values = matrix(nrow = nrow(data), ncol = length(colnames), dimnames = list(NULL, colnames))
    n = length(name)
    m = length(colnames)/n
    
    for (i in 1:nrow(data)) {
        for (j in 1:m) { 
            key_name = unlist(key[key$Treatment == data$Treatment[i], j + 1])
            if (is.na(key_name)) { 
                values[i,j] = NA
            } else if (match) {
                # if multiple names, they are equally spaced
                values[i,seq(j, n*m, m)] = 
                    unlist(ref_data[ref_data$Treatment == key_name & ref_data$Screen_rep_ID == data$Screen_rep_ID[i], name]);
            } else {
                values[i,j] = unlist(filter(ref_data, Treatment==key_name)[name])
            }
        }
    }
    
    return(values)
}

### quick function to calculate standard error ###
se = function(v) {
    return(sd(v, na.rm=T) / sqrt(length(v)))
}

### function for calculating the isserlis prediction of a metrix (prefix) on a dataframe (data) ###
isserlis = function(data, prefix) {
    return(data[, str_c(prefix, '_1')] * data[, str_c(prefix, '_23')] + 
               data[, str_c(prefix, '_2')] * data[, str_c(prefix, '_13')] + 
               data[, str_c(prefix, '_3')] * data[, str_c(prefix, '_12')] - 
               2 * data[, str_c(prefix, '_1')] * data[, str_c(prefix, '_2')] * data[, str_c(prefix, '_3')])
}

### calculates interaction of PI values ###
make_interaction = function(PIs) {
    return(PIs %>%
               mutate(I_12 = PI_norm_12 - PI_norm_1 * PI_norm_2,
                      I_13 = PI_norm_13 - PI_norm_1 * PI_norm_3,
                      I_23 = PI_norm_23 - PI_norm_2 * PI_norm_3,
                      I_123 = PI_norm_123 - PI_norm_1*PI_norm_2*PI_norm_3,
                      I_iss = PI_iss - PI_norm_1*PI_norm_2*PI_norm_3))
}

### returns the value in a vector with the minimum absolute value (returns actual value) ###
abs_min = function(x) {
    x[which.min(abs(x))]
}