import numpy as np

from moose.core.txtio import TxtIO
from moose.grids import Mesh1d
from moose.data import Metadata, Data



class ConnectionHistogram(TxtIO):
    def __init__(self, psiN, turns, nsamples, counts):
        self.psiN = psiN
        self.turns = turns
        self.nsamples = nsamples
        self.counts = counts
        if not counts.shape == (2, psiN.size, turns.size+1):
            raise(ValueError("cdf has unexpected shape {}, should be (2, {}, {})".format(counts.shape, psiN.size, turns.size+1)))


    @classmethod
    def _readtxt(cls, f, psin: int, turns: int, nsamples: int):
        n, m = psin, turns
        psiN = np.fromfile(f, dtype=float, count=n, sep=' ')
        turns = np.fromfile(f, dtype=float, count=m, sep=' ')
        counts = np.fromfile(f, dtype=int, count=n*(1+m)*2, sep=' ').reshape(2,n,(1+m))
        return cls(psiN, turns, nsamples, counts)


    @property
    def cdf_bwd(self):
        return self.counts[0,...] / self.nsamples


    @property
    def cdf_fwd(self):
        return self.counts[1,...] / self.nsamples


    @staticmethod
    def __select_direction(direction, ybwd, yfwd):
        if direction == "bwd":
            return ybwd
        elif direction == "fwd":
            return yfwd
        elif direction == "both":
            return (ybwd + yfwd) / 2
        else:
            raise(ValueError("invalid direction = '{}'".format(direction)))


    def profile(self, iturns, direction="both", raw_data=False):
        """Radial profile of field line connections within selected length (turns[iturns])."""

        y = self.__select_direction(direction, self.cdf_bwd[:,iturns], self.cdf_fwd[:,iturns])

        if raw_data: return self.psiN, y
        grid1d = Mesh1d(self.psiN, "Normalized Poloidal Flux")
        return Data(y, grid1d, Metadata("{Fraction of connecting field lines}", "Field line connections"))


    def pdf(self, ipsiN, direction="both", raw_data=False):
        """Probability density of field line connections from ipsiN-th flux surface."""

        turns   = np.hstack(([0], self.turns))
        pdf_bwd = (self.cdf_bwd[ipsiN,1:] - self.cdf_bwd[ipsiN,0:-1]) / (turns[1:] - turns[:-1])
        pdf_fwd = (self.cdf_fwd[ipsiN,1:] - self.cdf_fwd[ipsiN,0:-1]) / (turns[1:] - turns[:-1])
        pdf     = self.__select_direction(direction, pdf_bwd, pdf_fwd)
        x       = (turns[1:] + turns[:-1])/2 * np.pi

        if raw_data: return x, pdf
        grid1d = Mesh1d(x, "Toroidal turns")
        return Data(pdf, grid1d, Metadata("p_connect", "Probability density of field line connections"))


    def cdf(self, ipsiN, direction="both", raw_data=False):
        """Distribution of field line connections from ipsiN-th flux surface."""

        y = self.__select_direction(direction, self.cdf_bwd[ipsiN,:], self.cdf_fwd[ipsiN,:])
        x = np.hstack(([0], self.turns))

        if raw_data: return x, y
        grid1d = Mesh1d(x, "Toroidal turns")
        return Data(y, grid1d, Metadata("{Fraction of connecting field lines}", "Distribution of field line connections"))



if __name__ == "__main__":
    import argparse
    import matplotlib.pyplot as plt
    import sys

    # 0. main parser
    parser     = argparse.ArgumentParser(description="Visualize connection length statistics.")
    parser.add_argument("filename")
    subparsers = parser.add_subparsers(dest="plot")

    # 1. radial profile
    profile = subparsers.add_parser("profile", help="Radial profile of field line connections.")
    profile.add_argument("nturns", type=int, nargs='+', help="Max. length for field line connection [toroidal turns].")
    profile.add_argument("-d", "--direction", default="both", choices=["bwd", "fwd", "both"])

    # 2. pdf
    pdf = subparsers.add_parser("pdf", help="Probability density of field line connections.")
    pdf.add_argument("ipsiN", type=int, nargs='+', help="Radial index.")
    pdf.add_argument("-d", "--direction", default="both", choices=["bwd", "fwd", "both"])

    # 3. cdf
    cdf = subparsers.add_parser("cdf", help="Distribution of field line connections.")
    cdf.add_argument("ipsiN", type=int, nargs='+', help="Radial index.")
    cdf.add_argument("-d", "--direction", default="both", choices=["bwd", "fwd", "both"])


    #-------------------------------------------------------------------
    args = parser.parse_args()
    if not args.plot:
        parser.print_help()
        sys.exit()
    esc = ConnectionHistogram.loadtxt(args.filename)

    # 1. radial profile
    if args.plot == "profile":
        for iturns in args.nturns:
            esc.profile(iturns, args.direction).plot(label="{} turns".format(iturns))
    # 2. pdf
    elif args.plot == "pdf":
        for ipsiN in args.ipsiN:
            esc.pdf(ipsiN, args.direction).plot(label=esc.psiN[ipsiN])
    # 3. cdf
    elif args.plot == "cdf":
        for ipsiN in args.ipsiN:
            esc.cdf(ipsiN, args.direction).plot(label=esc.psiN[ipsiN])
    plt.legend()
    plt.show()
