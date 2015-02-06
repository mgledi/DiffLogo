preconditionPWM = function(pwm) {
    if (any(abs(1 - apply(pwm, 2, sum)) > 0.01)) { # check if the sum of each columns is near 1
        stop("Columns of PWM must add up to 1.0")
    }
}

preconditionTransformPWM = function(pwm) {
    if (class(pwm) == "pwm") {
        return (pwm@pwm)
    } else if (class(pwm) == "data.frame") {
        return (as.matrix(pwm))
    } else if (class(pwm) != "matrix") {
        print("pwm must be of class matrix or data.frame. Trying to convert")
        return (matrix(pwm,4,length(pwm)/4))
    }
    return(pwm)
}

preconditionPWMSameSize = function(pwm1, pwm2) {
    if(ncol(pwm1) != ncol(pwm2)) {  # check if the two PWMs have the same length
        stop("The two given PWMs must have the same dimension");
    }
}
