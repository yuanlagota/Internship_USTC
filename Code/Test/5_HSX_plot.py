import matplotlib.pyplot as plt
from moose.data import Dataset
from flare import model
from flare.analysis import boundary

####################################################################################################### 
                                            # BOUNDARY 1 # 
#######################################################################################################   

base = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/"
plots = "../../Data/Plots"
# --- Connection-length map on the launch flux surface (.dat) ---
connectionlength = Dataset.loadtxt(base + "lc_0.05.dat")     # grid rebuilt from the .grid in the header
connectionlength["Lc"].plot()                                 # Lc = Lc_neg + Lc_pos  [m]
plt.savefig(plots + "lc_map.png", dpi=200)
plt.show()

# --- Strike-point density on the PFC (.nc) ---
strike_points = Dataset.loadnc(base + "strike_point_density_0.05.nc")
strike_points["p"].plot()                                  # p = strike-point density [m^-2]
plt.savefig(plots + "strike_density.png", dpi=200)
plt.show()


# gather all strike points from 1st boundary and plot their location
# strike points are converted from mesh (u1, u2) coordinates to surface (phi, v) coordinates
for direction, color in zip(["fwd", "bwd"], ['r', 'b']):
    x, y = boundary.strike_points(strike_points, direction, 1)
    plt.scatter(x, y, color=color, s=0.01, alpha=0.2)

plt.xlabel("Toroidal angle [deg]")
plt.ylabel("Poloidal angle [deg]")
plt.show()
