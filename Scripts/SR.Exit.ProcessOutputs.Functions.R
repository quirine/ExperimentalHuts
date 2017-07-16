
read.and.process.chains <- function(output.file){
  chain = read.csv(output.file,header=TRUE,colClasses = 'numeric')
  chain = chain[-1,]
  chain = mcmc(chain)
  return(chain)
}