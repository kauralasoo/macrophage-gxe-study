import peer
import scipy
import argparse
import os

parser = argparse.ArgumentParser(description = "Apply PEER to a gene expression matrix.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input", help = "Path to the input matrix.")
parser.add_argument("--outdir", help = "Directory of the output files.")
parser.add_argument("--n_factors", help = "Number of hidden factors to detect.", type=int)

args = parser.parse_args()

#Make output directory
if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

#Import data
expr = scipy.loadtxt(args.input, delimiter = ",")

#Create model object
model = peer.PEER()
model.setPhenoMean(expr)
model.setAdd_mean(True)

#Set the number of hidden factors
model.setNk(args.n_factors)
model.update()

#Observe results
factors = model.getX()
weights = model.getW()
precision = model.getAlpha()
residuals = model.getResiduals()

#Save results to disk
scipy.savetxt(os.path.join(args.outdir, "factors.txt"), factors, delimiter = ",")
scipy.savetxt(os.path.join(args.outdir, "weights.txt"), weights, delimiter = ",")
scipy.savetxt(os.path.join(args.outdir, "precision.txt"), precision, delimiter = ",")
scipy.savetxt(os.path.join(args.outdir, "residuals.txt"), residuals, delimiter = ",")

