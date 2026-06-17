import argparse

from ..analysis import footprint_parameters



parser = argparse.ArgumentParser(prog="flare footprint_parameters", description="Evaluate characteristic parameters of magnetic footprint.")
parser.add_argument("filename")
parser.add_argument("-raw", action="store_true", help="Print data values in one line")
parser.add_argument("-separatrix", default=1.0, type=float, help="Alternative psiN threshold for footprint domain.")
args = parser.parse_args()


d, A, avR, minR, avL = footprint_parameters(args.filename, args.separatrix)
if args.raw:
    print(d, A, avR, minR, avL)
else:
    print("Footprint width [cm]:              {}".format(d))
    print("Footprint area [m²]:               {}".format(A))
    print("Average radial connection:         {}".format(avR))
    print("Deepest radial connection:         {}".format(minR))
    print("Average connection length [p.t.]:  {}".format(avL))
