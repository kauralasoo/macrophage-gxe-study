import peer
import scipy
import pylab

#Import data
expr = scipy.loadtxt("results/SL1344/cond_A_exprs.peer.txt", delimiter = ",")

#Create model object
model = peer.PEER()
model.setPhenoMean(expr)
model.setAdd_mean(True)
model.getPhenoMean().shape

#Set the number of hidden confounders
model.setNk(10)
model.update()

#Observe results
factors = model.getX()
weights = model.getW()
precision = model.getAlpha()
residuals = model.getResiduals()

#Save results to disk
scipy.savetxt("results/SL1344/cond_A.peer_factors.txt", factors, delimiter = ",")
scipy.savetxt("results/SL1344/cond_A.peer_weights.txt", weights, delimiter = ",")
scipy.savetxt("results/SL1344/cond_A.peer_precision.txt", precision, delimiter = ",")
scipy.savetxt("results/SL1344/cond_A.peer_residuals.txt", residuals, delimiter = ",")

#Plot precision
pylab.plot(precision)
pylab.show()
