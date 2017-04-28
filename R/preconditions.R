preconditionPWM = function(pwm) {
  if (any(abs(1 - apply(pwm, 2, sum)) > 0.01)) { # check if the sum of each columns is near 1
    stop("Columns of PWM must add up to 1.0")
  }
}
preconditionProbabilityVector = function(vec) {
  if (abs(1 - sum(vec)) > 0.01) { # check if the sum of each columns is near 1
    stop("Vector must add up to 1.0")
  }
}

preconditionTransformPWM = function(pwm, alphabet) {
    if (class(pwm) == "pwm") {
        return (pwm@pwm)
    } else if (class(pwm) == "data.frame") {
        return (as.matrix(pwm))
    } else if (class(pwm) != "matrix") {
        print("pwm must be of class matrix or data.frame. Trying to convert")
        return (matrix(pwm,alphabet$size,length(pwm)/alphabet$size))
    }
    return(pwm)
}

preconditionPWMSameSize = function(pwm1, pwm2) {
  if(ncol(pwm1) != ncol(pwm2) || nrow(pwm1) != nrow(pwm2)) {  # check if the two PWMs have the same length
    stop("The two given PWMs must have the same dimension");
  }
}
preconditionVectorSameSize = function(vec1, vec2) {
  if(length(vec1) != length(vec2)) {  # check if the two vectors have the same length
    stop("The two given vectors must have the same length");
  }
}

preconditionStackHeight = function(h) {
    if(is.na(h$height)) {
        stop("Height is not allowed to be NA.");
    }
}

preconditionBaseDistribution = function(distr) {
    if(any(is.na(distr))) {
        stop("A value in the given distribution is NA.");
    }
}