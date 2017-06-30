divergencePenaltyForUnaligned = function(pwm, unaligned_at_left,
                                         aligned_length,
                                         unaligned_at_right, divergence,
                                         base_distribution) {
    shift_divergence = 0
    if (unaligned_at_left > 0) {
       shift_divergence = shift_divergence +
                          pwmDivergence(    matrix( pwm[,1:unaligned_at_left], nrow=nrow(pwm)),
                                            baseDistributionPwm( unaligned_at_left, nrow(pwm), base_distribution = base_distribution), 
                                            divergence
                                        );
    }
    if (unaligned_at_right > 0) {
       aligned_length = ncol(pwm) - unaligned_at_right
       shift_divergence = shift_divergence +
                          pwmDivergence(    matrix(pwm[,(aligned_length+1):ncol(pwm)], nrow=nrow(pwm)), 
                                            baseDistributionPwm(ncol(pwm)-aligned_length, nrow(pwm), base_distribution = base_distribution), 
                                            divergence
                                        );
    }
    return(shift_divergence)
}

findBestShiftForPwms = function(static_pwm, shifted_pwm, divergence,  unaligned_penalty, base_distribution, length_normalization) {
    best_divergence = Inf
    best_shift = 0
    for (shift in 0:(ncol(static_pwm)-1)){
        intersection_length = min(ncol(static_pwm) - shift, ncol(shifted_pwm))
        shift_divergence = pwmDivergence(
                                matrix( static_pwm[, (shift+1):(shift+intersection_length)], nrow=nrow(static_pwm)),
                                matrix( shifted_pwm[,1:intersection_length], nrow=nrow(shifted_pwm)),
                                divergence=divergence)
        shift_divergence = shift_divergence+unaligned_penalty(
                                      static_pwm,
                                      unaligned_at_left = shift,
                                      intersection_length,
                                      unaligned_at_right = max(0, ncol(static_pwm) - (shift+intersection_length)),
                                      divergence,
                                      base_distribution);
        shift_divergence = shift_divergence+unaligned_penalty(
                                      shifted_pwm,
                                      unaligned_at_left = 0,
                                      intersection_length,
                                      unaligned_at_right = max(0, ncol(shifted_pwm) - intersection_length),
                                      divergence,
                                      base_distribution);
        if (length_normalization) {
           shift_divergence = shift_divergence / max(ncol(static_pwm), shift + ncol(shifted_pwm));
        }
        if (shift_divergence < best_divergence) {
           best_divergence = shift_divergence
           best_shift = shift
        }
    }
    return(list('divergence' = best_divergence, 'shift' = best_shift));
}

##' Finds best local alignment for two PWMs.
##'
##' @title Align pwms
##' @param pwm_left first PWM, a matrix of type matrix
##' @param pwm_right first PWM, a matrix of type matrix
##' @param divergence is a measure of difference between two pwm columns. Smaller is more similar. If you want to use non-uniform background distribution, provide your own function.
##' @param unaligned_penalty distance for unaligned columns at edges of matrixes. See divergencePenaltyForUnaligned as an example for providing your own function
##' @param try_reverse_complement If false the alignment will not be performed on reverse complements. If true, the input pwms should have column order of ACTG/ACGU.
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @param length_normalization If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
##' @return list of length two containing the alignment and the divergence
##' @author Lando Andrey
localPwmAlignment = function(pwm_left, pwm_right, divergence=shannonDivergence, unaligned_penalty=divergencePenaltyForUnaligned, try_reverse_complement=TRUE, base_distribution=NULL, length_normalization = FALSE) {
    no_change = list("shift"=0, "direction"="forward")
    result_alignment_vector = list()
    best_divergence = Inf

    alignment = findBestShiftForPwms(pwm_left, pwm_right, divergence=divergence,
                                     unaligned_penalty=unaligned_penalty,
                                     base_distribution,
                                     length_normalization=length_normalization)
    if (alignment$divergence < best_divergence) {
       result_alignment_vector[[1]] = no_change
       result_alignment_vector[[2]] = list("shift" = alignment$shift, "direction"="forward")
       result_divergence = alignment$divergence
       best_divergence = alignment$divergence
    }

    alignment = findBestShiftForPwms(pwm_right, pwm_left, divergence=divergence,
                                     unaligned_penalty=unaligned_penalty,
                                     base_distribution=base_distribution,
                                     length_normalization=length_normalization)
    if (alignment$divergence < best_divergence) {
        result_alignment_vector[[1]] = list("shift" = alignment$shift, "direction"="forward")
        result_alignment_vector[[2]] = no_change
        result_divergence = alignment$divergence
        best_divergence = alignment$divergence
    }
    
    if (try_reverse_complement) {
        alignment = findBestShiftForPwms(pwm_left, revCompPwm(pwm_right),
                                         divergence=divergence,
                                         unaligned_penalty=unaligned_penalty,
                                         base_distribution=base_distribution,
                                         length_normalization=length_normalization)
        if (alignment$divergence < best_divergence) {
            result_alignment_vector[[1]] = no_change
            result_alignment_vector[[2]] = list("shift" = alignment$shift, "direction"="reverse")
            result_divergence = alignment$divergence
            best_divergence = alignment$divergence
        }

        alignment = findBestShiftForPwms(revCompPwm(pwm_right), pwm_left,
                                         divergence=divergence,
                                         unaligned_penalty=unaligned_penalty,
                                         base_distribution=base_distribution,
                                         length_normalization=length_normalization)
        if (alignment$divergence < best_divergence) {
            result_alignment_vector[[1]] = list("shift" = alignment$shift, "direction"="forward")
            result_alignment_vector[[2]] = list("shift" = 0, "direction"="reverse")
            result_divergence = alignment$divergence
            best_divergence = alignment$divergence
        }
    }
    return(list("vector"=result_alignment_vector, "divergence"=result_divergence))
}

addBefore = function (matrix, to_add_length, base_distribution=NULL) {
    return(cbind(baseDistributionPwm(to_add_length, nrow(matrix), base_distribution), matrix))
}

addAfter = function (matrix, to_add_length, base_distribution=NULL) {
    return(cbind( matrix, baseDistributionPwm(to_add_length, nrow(matrix), base_distribution)))
}

##' Extends pwms by adding base_distribution to both sides, so that they keep aligned, but have equal length.
##'
##' @title Extend pwms with respect to alignment
##' @param pwms is a list of matrixes
##' @param alignment_vector is a list of shifts ($shift) and orientations ($direction)
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @return extended pwms
##' @export
##' @author Lando Andrey
##' @examples 
##' file1 = system.file("extdata/homer/Max.motif",   package = "DiffLogo")
##' file2 = system.file("extdata/homer/c-Myc.motif", package = "DiffLogo")
##' pwm1 = getPwmFromFile(file1)
##' pwm2 = getPwmFromFile(file2)
##' 
##' pwms <- list(pwm1, pwm2)
##' multiple_pwms_alignment = multipleLocalPwmsAlignment(pwms)
##' aligned_pwms = extendPwmsFromAlignmentVector(pwms, multiple_pwms_alignment$alignment$vector)
extendPwmsFromAlignmentVector = function(pwms, alignment_vector, base_distribution=NULL) {
    stopifnot(length(pwms)==length(alignment_vector))
    max_length = 0
    for (i in 1:length(pwms)) {
        if (alignment_vector[[i]]$direction == 'reverse') {
            pwms[[i]] = revCompPwm(pwms[[i]])
        }
        pwms[[i]] = addBefore(pwms[[i]], alignment_vector[[i]]$shift, base_distribution)
        max_length = max(max_length, ncol(pwms[[i]]))
    }
    for (i in 1:length(pwms)) {
        pwms[[i]] = addAfter(pwms[[i]], max_length - ncol(pwms[[i]]), base_distribution)
    }
    return(pwms)
}

##' Computes average pwm divergence from alignment vector for two datasets. This equals to average divergence between all pairs where one pwm comes from left set, and other comes from right
##'
##' @title Average divergence between two sets.
##' @param left_pwms_list is a list of matrixes
##' @param left_pwms_alignment is a list of shifts ($shift) and orientations ($direction)
##' @param right_pwms_list is a list of matrixes
##' @param right_pwms_alignment is a list of shifts ($shift) and orientations ($direction)
##' @param divergence divergence measure.
##' @return float
##' @author Lando Andrey
twoSetsAveragePwmDivergenceFromAlignmentVector = function(left_pwms_list,
                                                          left_pwms_alignment,
                                                          right_pwms_list,
                                                          right_pwms_alignment,
                                                          divergence = shannonDivergence){
    stopifnot(length(left_pwms_list)==length(left_pwms_alignment))
    stopifnot(length(right_pwms_list)==length(right_pwms_alignment))
    pwms_divergences = c()
    for (left_pwm_index in 1:length(left_pwms_list)) {
        for (right_pwm_index in 1:length(right_pwms_list)) {
            stopifnot(!is.null(left_pwms_list[left_pwm_index]))
            stopifnot(!is.null(right_pwms_list[right_pwm_index]))

            stopifnot(!is.null(left_pwms_alignment[left_pwm_index]))
            stopifnot(!is.null(right_pwms_alignment[right_pwm_index]))
            pwms = extendPwmsFromAlignmentVector(c(left_pwms_list[left_pwm_index], right_pwms_list[right_pwm_index]),
                                                 c(left_pwms_alignment[left_pwm_index], right_pwms_alignment[right_pwm_index]))
            pwms_divergences = c(pwms_divergences, pwmDivergence(pwms[[1]], pwms[[2]], divergence = divergence))
        }
    }
    return(mean(pwms_divergences))
}

findBestShiftForPwmSets = function(static_pwms_list, static_pwms_alignment,
                                   shifted_pwms_list, shifted_pwms_alignment,
                                   divergence = shannonDivergence,
                                   unaligned_penalty = divergencePenaltyForUnaligned) {
    stopifnot(length(static_pwms_list)==length(static_pwms_alignment$vector))
    stopifnot(length(shifted_pwms_list)==length(shifted_pwms_alignment$vector))
    best_divergence = Inf
    best_shift = 0
    max_possible_shift = 0
    result = list()
    for (i in 1:length(static_pwms_list)) {
        max_possible_shift = max(max_possible_shift,
                                 ncol(static_pwms_list[[i]])
                                     +static_pwms_alignment$vector[[i]]$shift)
    }
    for (shift in 0:max_possible_shift) {
        divergence_for_shift = twoSetsAveragePwmDivergenceFromAlignmentVector(
            static_pwms_list, static_pwms_alignment$vector,
            shifted_pwms_list, shifted_pwms_alignment$vector,
            divergence=divergence)
        if (divergence_for_shift < best_divergence) {
            best_divergence = divergence_for_shift
            result = list("vector" = c(static_pwms_alignment$vector, shifted_pwms_alignment$vector),
                          "divergence" = best_divergence)
        }
        for (i in 1:length(shifted_pwms_alignment$vector)) {
            shifted_pwms_alignment$vector[[i]]$shift = shifted_pwms_alignment$vector[[i]]$shift+1
        }
    }
    return(result)
}

##' Switches between 'forward' and 'reverse'
##' 
##' @title Switches between 'forward' and 'reverse'
##' @param direction either 'forward' or 'reverse'
##' @return either 'reverse' or 'forward'
switchDirection = function(direction) {
    if (direction=='forward') {
       return('reverse')
    } else if (direction=='reverse') {
       return('forward')
    }
    stopifnot(FALSE);
}

##' Returns alignment vector as if all pwm were reverted.
##'
##' @title Reverse for alignment vector
##' @param alignment_vector list of list which $shift and $orientation
##' @param pwms list of matrixes.
##' @return list - reversed alignment vector
##' @author Lando Andrey
reverseAlignmentVector = function(alignment_vector, pwms) {
    stopifnot(length(alignment_vector)==length(pwms));
    alignment_length = 0;
    for (i in 1:length(alignment_vector)) {
        alignment_vector[[i]]$direction = switchDirection(alignment_vector[[i]]$direction);
        stopifnot(alignment_vector[[i]]$shift >= 0);
        alignment_length = max(alignment_length, ncol(pwms[[i]]) + alignment_vector[[i]]$shift);
    }
    min_shift = 0;
    for (i in 1:length(alignment_vector)) {
        alignment_vector[[i]]$shift = alignment_length - ncol(pwms[[i]]) - alignment_vector[[i]]$shift;
        min_shift = min(min_shift, alignment_vector[[i]]$shift);
    }
    for (i in 1:length(alignment_vector)) {
        alignment_vector[[i]]$shift = -min_shift + alignment_vector[[i]]$shift;
    }
    return(alignment_vector)
}


##' Align two sets of pwms
##'
##' @title Multiple PWMs alignment
##' @param left_pwms_set list of pwms(matrixes)
##' @param left_alignment alignment of left_pwms_set. 
##' @param right_pwms_set list of pwms;
##' @param right_alignment alignment of right_pwms_set.
##' @param try_reverse_complement if true(default), also try reverse complement.
##' @return list - alignment of concatination of left_pwms_set and right_pwms_set
##' @author Lando Andrey
alignPwmSets = function(left_pwms_set, left_alignment, right_pwms_set, right_alignment, try_reverse_complement) {
    stopifnot(length(left_pwms_set)==length(left_alignment$vector))
    stopifnot(length(right_pwms_set)==length(right_alignment$vector))
    result = list()
    best_divergence = Inf
    alignment = findBestShiftForPwmSets(left_pwms_set, left_alignment, right_pwms_set, right_alignment)
    if (alignment$divergence < best_divergence) {
       result = alignment
       best_divergence = alignment$divergence
    }

    alignment = findBestShiftForPwmSets(right_pwms_set, right_alignment, left_pwms_set, left_alignment)
    if (alignment$divergence < best_divergence) {
       best_divergence = alignment$divergence
       result$vector = c(alignment$vector[(length(right_pwms_set)+1):(length(alignment$vector))],
                            alignment$vector[1:length(right_pwms_set)])
       result$divergence = best_divergence
    }

    if(try_reverse_complement) {
      right_alignment$vector = reverseAlignmentVector(right_alignment$vector, right_pwms_set)
      alignment = findBestShiftForPwmSets(left_pwms_set, left_alignment, right_pwms_set, right_alignment)
      if (alignment$divergence < best_divergence) {
         result = alignment
         best_divergence = alignment$divergence
      }
    }

    alignment = findBestShiftForPwmSets(right_pwms_set, right_alignment, left_pwms_set, left_alignment)
    if (alignment$divergence < best_divergence) {
       best_divergence = alignment$divergence
       result$vector = c(alignment$vector[(length(right_pwms_set)+1):(length(alignment$vector))], alignment$vector[1:length(right_pwms_set)])
       result$divergence = best_divergence
    }
    return(result)
}

##' Creates a distance matrix for pwms
##'
##' @title Multiple PWMs alignment
##' @param pwms list of pwms
##' @param diagonal_value value to put on diagonal.
##' @param bottom_default_value value to put on bottom triangle. Set to NULL to get symmetric distance matrix.
##' @param divergence divergence measure.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @param length_normalization is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @return list
##' @author Lando Andrey
pwmsDistanceMatrix = function(pwms,
                              diagonal_value = 0,
                              bottom_default_value = NULL,
                              divergence=shannonDivergence,
                              unaligned_penalty=divergencePenaltyForUnaligned,
                              try_reverse_complement=TRUE, base_distribution=NULL,
                              length_normalization = FALSE){
    distance_matrix = matrix(0, ncol=length(pwms), nrow=length(pwms))
    for (i in 1:(length(pwms)-1)) {
        distance_matrix[[i, i]] = diagonal_value
        for (j in (i+1):length(pwms)) {
            distance_matrix[[i, j]] = localPwmAlignment(
                pwms[[i]], pwms[[j]],
                divergence=divergence,
                unaligned_penalty=unaligned_penalty,
                try_reverse_complement=try_reverse_complement,
                base_distribution=base_distribution,
                length_normalization = length_normalization)$divergence
            if (is.null(bottom_default_value)) {
                distance_matrix[[j, i]] = distance_matrix[[i, j]]
            } else {
                distance_matrix[[j, i]] = bottom_default_value
            }
        }
    }
    distance_matrix[[length(pwms), length(pwms)]] = diagonal_value
    return(distance_matrix)
}

minimumMatrixElementIndexes = function(matrix) {
    minimum_element = which.min(matrix)
    col_of_minimum = round(minimum_element/ncol(matrix)-0.5)+1
    row_of_minimum = minimum_element - ((col_of_minimum-1) * (ncol(matrix)))
    first_index = min(col_of_minimum, row_of_minimum)
    second_index = max(col_of_minimum, row_of_minimum)
    return(list(first_index, second_index))
}

addLastTreeNodeToDistanceMatrix = function(distance_matrix, tree_nodes, try_reverse_complement) {
    distance_matrix = cbind(distance_matrix, rep(Inf, ncol(distance_matrix)))
    distance_matrix = rbind(distance_matrix, rep(Inf, ncol(distance_matrix)))
    last = ncol(distance_matrix)
    if (last > 0) {
        for (i in 1:(last-1)) {
            distance_matrix[[i, last]] = alignPwmSets( tree_nodes[[last]]$pwms, tree_nodes[[last]]$pwms_alignment,
                tree_nodes[[i]]$pwms,tree_nodes[[i]]$pwms_alignment, try_reverse_complement)$divergence
        }
    }
    return(distance_matrix)
}

joinTwoNodesInAlignmentTree = function(distance_matrix, tree_nodes, join_counter, try_reverse_complement) {
    min_pwm_row_col = minimumMatrixElementIndexes(distance_matrix)
    first_pwm_index = min_pwm_row_col[[1]]
    second_pwm_index = min_pwm_row_col[[2]]
    pwms_alignment = alignPwmSets(tree_nodes[[first_pwm_index]]$pwms,
                                  tree_nodes[[first_pwm_index]]$pwms_alignment,
                                  tree_nodes[[second_pwm_index]]$pwms,
                                  tree_nodes[[second_pwm_index]]$pwms_alignment,
                                  try_reverse_complement)
    stopifnot(length(c(tree_nodes[[first_pwm_index]]$pwms, tree_nodes[[second_pwm_index]]$pwms)) == length(pwms_alignment$vector))
    tree_nodes = c(tree_nodes,
                   list(
                       list("pwms"=c(tree_nodes[[first_pwm_index]]$pwms, tree_nodes[[second_pwm_index]]$pwms),
                            "left_child"=tree_nodes[[first_pwm_index]],
                            "right_child"=tree_nodes[[second_pwm_index]],
                            "num_leafs"=tree_nodes[[first_pwm_index]]$num_leaf + tree_nodes[[second_pwm_index]]$num_leaf,
                            "node_number"=join_counter,
                            "pwms_alignment" = pwms_alignment
                          )
                    )
                );
    join_counter = join_counter+1
    tree_nodes[[first_pwm_index]] = NULL
    tree_nodes[[second_pwm_index-1]] = NULL
    distance_matrix = matrix(distance_matrix[-c(first_pwm_index, second_pwm_index), -c(first_pwm_index, second_pwm_index)],
        nrow = nrow(distance_matrix)-2, ncol = ncol(distance_matrix)-2)
    if (ncol(distance_matrix) != 0) {
        distance_matrix = addLastTreeNodeToDistanceMatrix(distance_matrix, tree_nodes, try_reverse_complement)
    }
    return(list("distance_matrix"=distance_matrix, "tree_nodes"=tree_nodes))
}

alignmentTree = function(pwms, distance_matrix, try_reverse_complement) {
    tree_nodes = list()
    for (i in 1:length(pwms)) {
        tree_nodes = c(tree_nodes,
                       list(
                            list(   "pwms"=list(pwms[[i]]),
                                    "left_child"=NULL,
                                    "right_child"=NULL,
                                    "num_leafs"=1,
                                    "pwm_index"=i,
                                    "pwms_alignment"=list("vector"=list(list("direction"="forward", "shift"=0)))
                            )
                        )
                    );
    }
    join_counter = 1
    for (step in 1:(length(pwms)-1)) {
        joined = joinTwoNodesInAlignmentTree(distance_matrix, tree_nodes, join_counter, try_reverse_complement)
        distance_matrix = joined$distance_matrix
        tree_nodes = joined$tree_nodes
        join_counter = join_counter+1
    }
    return(tree_nodes[[1]])
}

alignmentTreeLeftRightTriversal = function(node) {
    merge = list()
    order = list()
    height = list()
    leftRightTreeTriversalHelper = function(node){
        if (is.null(node$left_child)) {
            order <<- c(order, node$pwm_index)
            return(-(node$pwm_index))
        } else {
            left_index = leftRightTreeTriversalHelper(node$left_child)
            right_index = leftRightTreeTriversalHelper(node$right_child)
            merge <<- c(merge, list(c(left_index, right_index)))
            height <<- c(height, list(node$pwms_alignment$divergence))
            #height <<- c(height, list(node$node_number))
        }
        return(length(merge))
    }
    leftRightTreeTriversalHelper(node)
    return(list("merge"=merge, "order"=unlist(order), "height"=height))
}


##' Creates a multiple alignment of pwms
##'
##' @title Multiple PWMs alignment
##' @param pwms list of pwms
##' @param divergence Divergence measure.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @param length_normalization If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
##' @export
##' @return list
##' @author Lando Andrey
##' @examples 
##' file1 = system.file("extdata/homer/Max.motif",   package = "DiffLogo")
##' file2 = system.file("extdata/homer/c-Myc.motif", package = "DiffLogo")
##' pwm1 = getPwmFromFile(file1)
##' pwm2 = getPwmFromFile(file2)
##' 
##' multiple_pwms_alignment = multipleLocalPwmsAlignment(list(pwm1, pwm2))
multipleLocalPwmsAlignment = function(
  pwms,
  divergence=shannonDivergence,
  unaligned_penalty=divergencePenaltyForUnaligned,
  try_reverse_complement=TRUE,
  base_distribution=NULL,
  length_normalization = FALSE)
{
    distance_matrix = pwmsDistanceMatrix(
        pwms,
        bottom_default_value = Inf,
        diagonal_value = Inf,
        divergence=divergence,
        unaligned_penalty=unaligned_penalty,
        try_reverse_complement=try_reverse_complement,
        base_distribution=base_distribution,
        length_normalization = length_normalization);
    alignment_tree_nodes = alignmentTree(pwms, distance_matrix, try_reverse_complement)
    traversal_result = alignmentTreeLeftRightTriversal(alignment_tree_nodes)
    alignment_tree_nodes$pwms_alignment$vector = alignment_tree_nodes$pwms_alignment$vector[ sort.list( unlist( traversal_result$order))];
    alignment_length = 0;
    for (i in 1:length(pwms)) {
        alignment_length = max( alignment_length, ncol(pwms[[i]]) + alignment_tree_nodes$pwms_alignment$vector[[i]]$shift);
    }
    return(list("alignment"=alignment_tree_nodes$pwms_alignment,
                "order"=unlist(traversal_result$order),
                "merge"=t(matrix(unlist(traversal_result$merge), 2, byrow=FALSE)),
                "height"=unlist(traversal_result$height),
                "raw_tree"=list(alignment_tree_nodes),
                "distance_matrix"=distance_matrix,
                "alignment_length"=alignment_length
                ))
}
