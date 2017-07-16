# User defined parameters for SR.Exit exercise

# Parameters used for data simulation -------------------------------------
movement.tied.to.exit = TRUE            # set to true when one movement rate is defined and r determines the proportion of movement into exit traps rather than to neighboring huts
Qs <- c(0.005, 0.005, 0.005)            # q_AorE, q_BorD, q_C: movement rates between huts
Ps <- c(0.7, 0.5)                       # p_BorD, p_C (p_AorE is always 1 and thus not defined) : probability of moving away from SR-hut (C)
Rs <- c(0.1, 0.1, 0.1)
# Rs <- c(0.002, 0.002, 0.002) #c(0.003, 0.003, 0.003) #c(0.001, 0.001, 0.001)  # r_AorE, r_BorD, r_C : exit rates
Ks <- c(0.0001, 0.0001, 0.0001)         # k_AorE, k_BorD, k_C : knockdown rates
u <- 0.0007
delta = 0.1                             # stepsize for ode
Interval <- 30                          # measurement interval in minutes (60 minutes as used for kd-monitoring in experimental hut studies)
Duration <- 750                         # 12.5 hours (from 5.30 am to 6 pm)    
Mosquitoes <- 1000                      # Number of mosquitoes released per hut
kd.interval.factor <- 2 
first.measurement <- 30
num.exp.days <- 5
# Parameters for MCMC -----------------------------------------------------
# we need eleven parameters to fit this model. Start of with using the simulated parameters as starting values
defaults <- c( q_AorE = Qs[1], q_BorD = Qs[2], q_C = Qs[3], 
               p_BorD = Ps[1], p_C = Ps[2], 
               r_AorE = Rs[1], r_BorD = Rs[2], r_C = Rs[3], 
               k_AorE = Ks[1], k_BorD = Ks[2], k_C = Ks[3],
               u = u) ## initial condition
startvalue <- defaults
# fix = c(1,2,3,5,6,7,8,9,10,11)        # vector that allows you to fix specific parameters. The order being the same as seen in 'defaults'
fix = NULL
Qfix = FALSE
iter = 50000                          # number of iterations
wr.freq = round(iter/1000)              # how frequently would you like the parameters to be printed to file?
mcmc.output.file = 'output.ltfu.csv'       # file that prints MCMC results 
sd = c(rep(0.005,3),rep(0.2,2),rep(0.2,3),rep(0.0001,3),0.004)
draw.probs = c(rep(1/11,4),0,rep(1/11,7))
relative.sds = FALSE
adaptive = FALSE
adaptive.paramaters = NULL
LTFU = TRUE
p.prior = 'Flat'
q.prior = 'Flat'

