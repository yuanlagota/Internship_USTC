import argparse

from flare import model, analysis



parser = argparse.ArgumentParser(prog="flare find_resonances", description="Scan equilibrium q-propfile for position of resonances")
parser.add_argument("n", type=int, help="toroidal mode number")
args = parser.parse_args()



model.load()
resonances = analysis.find_resonances(args.n)
for m in resonances:
    print(m, *resonances[m])
