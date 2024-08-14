ProduceOutput <- function(){
  Summary.Output.Mean[k,2] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,2]))
  Summary.Output.Mean[k,3] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,3]))
  Summary.Output.Mean[k,4] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,4]))
  Summary.Output.Mean[k,5] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,5]))
  Summary.Output.Mean[k,6] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,6]))
  Summary.Output.Mean[k,7] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,7]))
  Summary.Output.Mean[k,8] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,8]))
  Summary.Output.Mean[k,9] <<- abs(as.numeric(Output[2,10]) - as.numeric(Output[2,9]))
  Summary.Output.Mean[k,10] <<- NA
  Summary.Output.CI[k,2] <<- as.numeric(Output[4,2]) - as.numeric(Output[3,2])
  Summary.Output.CI[k,3] <<- as.numeric(Output[4,3]) - as.numeric(Output[3,3])
  Summary.Output.CI[k,4] <<- as.numeric(Output[4,4]) - as.numeric(Output[3,4])
  Summary.Output.CI[k,5] <<- as.numeric(Output[4,5]) - as.numeric(Output[3,5])
  Summary.Output.CI[k,6] <<- as.numeric(Output[4,6]) - as.numeric(Output[3,6])
  Summary.Output.CI[k,7] <<- as.numeric(Output[4,7]) - as.numeric(Output[3,7])
  Summary.Output.CI[k,8] <<- as.numeric(Output[4,8]) - as.numeric(Output[3,8])
  Summary.Output.CI[k,9] <<- as.numeric(Output[4,9]) - as.numeric(Output[3,9])
  Summary.Output.CI[k,10] <<- as.numeric(Output[4,10]) - as.numeric(Output[3,10])
}

