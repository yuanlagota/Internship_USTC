from dataclasses import dataclass
import netCDF4
import numpy as np
from os.path import basename, dirname, join

from moose.grids import Uqmesh
from moose.geometry import Axisurf, Torosurf, Hypersurf3d, Quad

from .. import f2py


cm_to_m = 0.01



class Rzmesh(Uqmesh):
    """
    Cross-section of Mmesh with index map for field lines and flux tubes.
    """

    def __init__(self, iphi, iside, nodes, quads, next_cell, aux_nodes, iline, itube):
        super().__init__(nodes, quads, next_cell, aux_nodes, x1label="r [m]", x2label="z [m]")
        self.iphi = iphi
        self.iside = iside
        self.iline = iline
        self.itube = itube


    @classmethod
    def _encoded_type(cls):
        return "uqmesh"


    @property
    def _metadata(self):
        return super()._metadata | {"iphi": self.iphi, "iside": self.iside}



@dataclass
class Mmesh:
    """Unstructured magnetic flux tube mesh"""

    symmetry: int   #: Toroidal symmetry of domain (< 0 for half-field period with stellarator symmetry).

    nphi: int       #: Number of planes in toroidal direction.

    nzones: int     #: Number of zones.

    nnodes: int     #: Number of mesh nodes.

    nlines: int     #: Number of field lines along mesh nodes.

    ntubes: int     #: Number of flux tubes.

    nbsect: int     #: Number of virtual field lines from bisection of flux tube sides.

    nxmaps: int     #: Number of flux tube side surfaces with non-default neighbor relations.


    @classmethod
    def allocate(cls, symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps):
        """Allocate arrays for mesh."""
        self = cls(symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps)
        self.phi = np.zeros(self.nphi)
        self.iphi_zone = np.zeros((self.nzones, 2), dtype=int)
        self.x = np.zeros((self.nnodes, 2))
        self.g = np.zeros((self.nnodes, 2))
        self.b = np.zeros(self.nnodes)
        self.corner = np.zeros((self.ntubes, 4), dtype=int)
        self.next_tube = np.zeros((self.ntubes, 4, 2), dtype=int)
        self.rparam_tmap = np.zeros((self.ntubes, 2, 16))
        self.iparam_tmap = np.zeros((self.ntubes, 2, 2), dtype=int)
        self.izone_tube = np.zeros(self.ntubes, dtype=int)
        self.izone_line = np.zeros(self.nlines, dtype=int)
        self.inode_offset = np.zeros(self.nlines+self.nbsect, dtype=int)
        self.bsect = np.zeros((self.nbsect, 2), dtype=int)
        self.iparam_xmap = np.zeros((self.nxmaps, 2), dtype=int)
        self.tmap = 0
        return self


    @classmethod
    def loadnc(cls, filename):
        with netCDF4.Dataset(filename, 'r') as ncfile:
            return cls.readnc(ncfile)


    @classmethod
    def readnc(cls, ncfile):
        symmetry = ncfile.symmetry
        nphi = ncfile.dimensions['nphi'].size
        nzones = ncfile.dimensions['nzones'].size
        nnodes = ncfile.dimensions['nnodes'].size
        nlines = ncfile.dimensions['nlines'].size
        ntubes = ncfile.dimensions['ntubes'].size
        nbsect = ncfile.dimensions['nbsect'].size
        nxmaps = ncfile.dimensions['nxmaps'].size
        self = cls(symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps)

        # toroidal domain
        self.phi = ncfile['phi'][:]
        self.iphi_zone = ncfile['iphi_zone'][:]

        # field lines
        self.x = ncfile['x'][:]
        self.g = ncfile['g'][:]
        self.b = ncfile['b'][:]
        self.izone_line = ncfile['izone_line'][:] - 1
        self.inode_offset = ncfile['inode_offset'][:]
        if nbsect > 0:
            self.bsect = ncfile['bsect'][:]
            self.inode_offset = np.roll(self.inode_offset, -nbsect)
        else:
            self.bsect = np.zeros(0, dtype=int)


        # flux tubes
        self.corner = ncfile['corner'][:]
        self.izone_tube = ncfile['izone_tube'][:] - 1
        self.next_tube = ncfile['next_tube'][:]
        self.tmap = ncfile.tmap
        if nxmaps > 0:
            self.iparam_xmap = ncfile['iparam_xmap'][:]
        else:
            self.iparam_xmap = np.zeros(0, dtype=int)


        # torosf_map
        if self.tmap == 1:
            self.rparam_tmap = ncfile['rparam_tmap'][:]
            self.iparam_tmap = ncfile['iparam_tmap'][:]
        else:
            self.rparam_tmap = np.zeros((self.ntubes, 2, 16))
            self.iparam_tmap = np.zeros((self.ntubes, 2, 2), dtype=int)
        return self


    def lbound_line(self, iline):
        if iline >= 0:
            return self.iphi_zone[self.izone_line[iline], 0]
        else:
            k = self.bsect[-iline-1,:]
            return max(self.lbound_line(k[0]), self.lbound_line(k[1]))


    def ubound_line(self, iline):
        if iline >= 0:
            return self.iphi_zone[self.izone_line[iline], 1]
        else:
            k = self.bsect[-iline-1,:]
            return min(self.ubound_line(k[0]), self.ubound_line(k[1]))


    def lbound_tube(self, itube):
        return self.iphi_zone[self.izone_tube[itube], 0]


    def ubound_tube(self, itube):
        return self.iphi_zone[self.izone_tube[itube], 1]


    def node_index(self, iline, iphi):
        return self.inode_offset[iline] + iphi - self.lbound_line(iline)


    def rzcoords(self, iline, iphi):
        if iline >= 0:
            return self.x[self.node_index(iline, iphi), :]
        else:
            k = self.bsect[-iline+1,:]
            return (self.rzcoords(k[0], iphi) + self.rzcoords(k[1], iphi)) / 2


    def xsect(self, itube, iphi):
        x1 = self.rzcoords(self.corner[itube, 0], iphi)
        x2 = self.rzcoords(self.corner[itube, 1], iphi)
        x3 = self.rzcoords(self.corner[itube, 2], iphi)
        x4 = self.rzcoords(self.corner[itube, 3], iphi)
        return Quad(x1, x2, x3, x4)


    @classmethod
    def import_emc3(cls, g):
        """Import mesh from EMC3."""

        # 1. layout (map EMC3 zones to mmesh zones)
        zone_map = np.zeros(g.nzonet, dtype=int)
        zones = []
        nphi = 0
        nzones = 0
        nnodes = g.grid_p_os[-1]
        nlines = 0
        ntubes = 0
        iline_offset = np.zeros(g.nzonet, dtype=int)
        itube_offset = np.zeros(g.nzonet, dtype=int)
        for iz in range(g.nzonet):
            iline_offset[iz] = nlines
            itube_offset[iz] = ntubes
            nlines += g.srf_radi[iz] * g.srf_polo[iz]
            ntubes += g.zon_radi[iz] * g.zon_polo[iz]

            phi = tuple(g.phi_plane[g.phi_pl_os[iz]:g.phi_pl_os[iz+1]])
            if not phi in zones:
                zone_map[iz] = len(zones)
                zones.append(phi)
                nzones += 1
                nphi = nphi + g.srf_toro[iz]
            else:
                zone_map[iz] = zones.index(phi)
        nphi -= nzones - 1
        self = cls.allocate(g.symmetry, nphi, nzones, nnodes, nlines, ntubes, 0, 0)

        # 2. toroidal domain (map EMC3-phi to mmesh-phi)
        i = 0
        for iz in range(nzones):
            n = len(zones[iz])
            self.phi[i:i+n] = zones[iz]
            self.iphi_zone[iz,0] = i
            self.iphi_zone[iz,1] = i + n - 1
            i += n - 1
        if any(i >= j for i, j in zip(self.phi, self.phi[1:])):
            print(self.phi)
            print(self.iphi_zone)
            raise(ValueError("phi array is not strictly increasing"))

        # 3. field lines
        iline = 0
        inode = 0
        for iz in range(g.nzonet):
            for ir in range(g.srf_radi[iz]):
                for ip in range(g.srf_polo[iz]):
                    self.inode_offset[iline] = inode
                    self.x[inode:inode+g.srf_toro[iz],0] = g.zone[iz].r[:,ip,ir] * cm_to_m
                    self.x[inode:inode+g.srf_toro[iz],1] = g.zone[iz].z[:,ip,ir] * cm_to_m
                    self.izone_line[iline] = zone_map[iz]
                    iline += 1
                    inode += g.srf_toro[iz]

        # 4. flux tubes
        # 4.1. corners, side surfaces
        def map_tube(imap):
            ipair = g.map_convect_r[0,imap]
            iz = g.map_convect_r[1,ipair]
            ir = g.map_convect_r[2,ipair]
            ip = g.map_convect_r[3,ipair]
            return itube_offset[iz] + ir * g.zon_polo[iz] + ip

        itube = 0
        for iz in range(g.nzonet):
            icore = 2 if g.idsurp[g.nps_off[iz]] == 1 else 0
            for ir in range(g.zon_radi[iz]):
                for ip in range(g.zon_polo[iz]):
                    self.corner[itube,0] = iline_offset[iz] + ir * g.srf_polo[iz] + ip
                    self.corner[itube,1] = self.corner[itube,0] + 1
                    self.corner[itube,2] = self.corner[itube,1] + g.srf_polo[iz]
                    self.corner[itube,3] = self.corner[itube,0] + g.srf_polo[iz]

                    self.next_tube[itube,:,0] = 1
                    self.next_tube[itube,0,1] = itube - g.zon_polo[iz]
                    self.next_tube[itube,1,1] = itube + 1
                    self.next_tube[itube,2,1] = itube + g.zon_polo[iz]
                    self.next_tube[itube,3,1] = itube - 1

                    # lower radial boundary
                    if ir == 0:
                        isurf = ir + ip * g.srf_radi[iz] + g.nrs_off[iz]
                        istat = g.idsurr[isurf]
                        if istat > 0:
                            self.next_tube[itube,0,1] = map_tube(istat // 10)
                        else:
                            self.next_tube[itube,0,0] = 0
                            self.next_tube[itube,0,1] = 1 + iz*10

                    elif ir == g.r_surf_pl_trans_range[0,iz]:
                        self.next_tube[itube,0,0] = 0
                        self.next_tube[itube,0,1] = 5 + icore + iz*10

                    # upper radial boundary
                    if ir == g.zon_radi[iz]-1:
                        isurf = ir + 1 + ip * g.srf_radi[iz] + g.nrs_off[iz]
                        istat = g.idsurr[isurf]
                        if istat > 0:
                            self.next_tube[itube,2,1] = map_tube(istat // 10)
                        else:
                            self.next_tube[itube,2,0] = 0
                            self.next_tube[itube,2,1] = 3 + iz*10

                    elif ir == g.r_surf_pl_trans_range[1,iz]-1:
                        self.next_tube[itube,2,0] = 0
                        self.next_tube[itube,2,1] = 6 + icore + iz*10

                    # lower poloidal boundary
                    if ip == 0:
                        isurf = ir + ip * g.zon_radi[iz] + g.nps_off[iz]
                        istat = g.idsurp[isurf]
                        if istat == 1:
                            self.next_tube[itube,3,1] += g.zon_polo[iz]
                        else:
                            self.next_tube[itube,3,0] = 0
                            self.next_tube[itube,3,1] = 4 + iz*10

                    # upper poloidal boundary
                    if ip == g.zon_polo[iz]-1:
                        isurf = ir + (ip + 1) * g.zon_radi[iz] + g.nps_off[iz]
                        istat = g.idsurp[isurf]
                        if istat == 1:
                            self.next_tube[itube,1,1] -= g.zon_polo[iz]
                        else:
                            self.next_tube[itube,1,0] = 0
                            self.next_tube[itube,1,1] = 2 + iz*10

                    self.izone_tube[itube] = zone_map[iz]
                    itube += 1

        # 4.2. up/down symmetric surfaces
        nupdown = 0
        for iz in range(g.nzonet):
            for ip in range(g.zon_polo[iz]):
                for ir in range(g.zon_radi[iz]):
                    itube = itube_offset[iz] + ir * g.zon_polo[iz] + ip
                    for iside, it in zip(range(2), [0, g.zon_toro[iz]]):
                        isurf = ir + (ip + it * g.zon_polo[iz]) * g.zon_radi[iz] + g.nts_off[iz]
                        istat = g.idsurt[isurf]
                        if istat == 2:
                            nupdown += 1
                            iphi = self.iphi_zone[iz,0] + it
                            xsect = self.xsect(itube, iphi)
                            c = xsect.interp_params
                            w = xsect.inverse_params(c)
                            s = xsect.xstep_params(c, w)
                            self.rparam_tmap[itube, iside, 0:8] = c.flatten()
                            self.rparam_tmap[itube, iside, 8:14] = s.flatten()
                            self.iparam_tmap[itube, iside, 0] = itube + g.zon_polo[iz] - 1 - 2*ip

        if not nupdown + g.wks_map == ntubes*2:
            print("ntubes = ", ntubes)
            print("wks_map = ", g.wks_map)
            raise(RuntimeError("unexpected number of mapping surfaces"))

        # 4.3. toroidal mapping surfaces
        for imap in range(g.total_map_sf_t):
            iz = g.zone_nr_map_t[imap]
            for ip in range(g.map_srf_np2_t[imap]):
                for ir in range(g.map_srf_nr1_t[imap]):
                    ipl = ir + ip * g.map_srf_nr1_t[imap] + g.ic_map_offset[imap]
                    itube = itube_offset[iz] + ir * g.zon_polo[iz] + ip
                    iside = 1 if g.t_sf_nr_map_t[imap] > 0 else 0

                    imap_pair = g.index_trans[0,ipl] - 1
                    jr = g.index_trans[1,ipl]
                    jp = g.index_trans[2,ipl]
                    jz = g.zone_nr_map_t[imap_pair]
                    # interp_params
                    self.rparam_tmap[itube, iside,  0] = g.coord_trans[ 1, ipl] * cm_to_m
                    self.rparam_tmap[itube, iside,  1] = g.coord_trans[ 5, ipl] * cm_to_m
                    self.rparam_tmap[itube, iside,  2] = g.coord_trans[ 3, ipl] * cm_to_m
                    self.rparam_tmap[itube, iside,  3] = g.coord_trans[ 7, ipl] * cm_to_m
                    self.rparam_tmap[itube, iside,  4] = g.coord_trans[ 2, ipl] * cm_to_m
                    self.rparam_tmap[itube, iside,  5] = g.coord_trans[ 6, ipl] * cm_to_m
                    self.rparam_tmap[itube, iside,  6] = g.coord_trans[ 4, ipl] * cm_to_m
                    self.rparam_tmap[itube, iside,  7] = g.coord_trans[ 8, ipl] * cm_to_m
                    # xstep_params
                    self.rparam_tmap[itube, iside,  8] = g.coord_trans[12, ipl] / cm_to_m
                    self.rparam_tmap[itube, iside,  9] =-g.coord_trans[11, ipl] / cm_to_m
                    self.rparam_tmap[itube, iside, 10] =-g.coord_trans[10, ipl] / cm_to_m
                    self.rparam_tmap[itube, iside, 11] = g.coord_trans[ 9, ipl] / cm_to_m
                    self.rparam_tmap[itube, iside, 12] = g.coord_trans[14, ipl]
                    self.rparam_tmap[itube, iside, 13] = g.coord_trans[13, ipl]
                    # xi-map
                    self.rparam_tmap[itube, iside, 14] = g.coord_trans[16, ipl]
                    self.rparam_tmap[itube, iside, 15] = g.coord_trans[15, ipl]
                    self.iparam_tmap[itube, iside, 0] = itube_offset[jz] + jr * g.zon_polo[jz] + jp
        self.tmap = 1

        return self


    def savenc(self, filename):
        """Save mesh as netcdf file."""
        ncfile = netCDF4.Dataset(filename, 'w')
        ncfile.symmetry = self.symmetry
        ncfile.tmap = self.tmap

        # supporting dimensions
        dim = {
            2: ncfile.createDimension('dim_0002', 2),
            4: ncfile.createDimension('dim_0004', 4),
            16: ncfile.createDimension('dim_0016', 16)
            }

        # named dimensions
        dim['nphi'] = ncfile.createDimension('nphi', self.nphi)
        dim['nzones'] = ncfile.createDimension('nzones', self.nzones)
        dim['nnodes'] = ncfile.createDimension('nnodes', self.nnodes)
        dim['nlines'] = ncfile.createDimension('nlines', self.nlines)
        dim['ntubes'] = ncfile.createDimension('ntubes', self.ntubes)
        dim['nbsect'] = ncfile.createDimension('nbsect', self.nbsect)
        dim['nxmaps'] = ncfile.createDimension('nxmaps', self.nxmaps)
        dim['noffsets'] = ncfile.createDimension('noffsets', self.nlines+self.nbsect)

        # toroidal domain variables
        ncfile.createVariable('phi', np.float64, ('nphi',))
        ncfile.createVariable('iphi_zone', np.int32, ('nzones', 'dim_0002'))

        # field line variables
        ncfile.createVariable('x', np.float64, ('nnodes', 'dim_0002'))
        ncfile.createVariable('g', np.float64, ('nnodes', 'dim_0002'))
        ncfile.createVariable('b', np.float64, ('nnodes',))
        ncfile.createVariable('izone_line', np.int32, ('nlines',))
        ncfile.createVariable('inode_offset', np.int32, ('noffsets',))
        if self.nbsect > 0:
            ncfile.createVariable('bsect', np.int32, ('nbsect', 'dim_0002'))

        # flux tube variables
        ncfile.createVariable('corner', np.int32, ('ntubes', 'dim_0004'))
        ncfile.createVariable('next_tube', np.int32, ('ntubes', 'dim_0004', 'dim_0002'))
        ncfile.createVariable('izone_tube', np.int32, ('ntubes',))
        if self.nxmaps > 0:
            ncfile.createVariable('iparam_xmap', np.int32, ('nxmaps', 'dim_0002'))
        if self.tmap > 0:
            ncfile.createVariable('rparam_tmap', np.float64, ('ntubes', 'dim_0002', 'dim_0016'))
            ncfile.createVariable('iparam_tmap', np.int32, ('ntubes', 'dim_0002', 'dim_0002'))


        # write data
        ncfile['phi'][:] = self.phi
        ncfile['iphi_zone'][:] = self.iphi_zone
        ncfile['x'][:] = self.x
        ncfile['g'][:] = self.g
        ncfile['b'][:] = self.b
        ncfile['izone_line'][:] = self.izone_line + 1
        ncfile['inode_offset'][:] = self.inode_offset
        if self.nbsect > 0:
            ncfile['bsect'][:] = self.bsect
        ncfile['corner'][:] = self.corner
        ncfile['next_tube'][:] = self.next_tube
        ncfile['izone_tube'][:] = self.izone_tube + 1
        if self.nxmaps > 0:
            ncfile['iparam_xmap'][:] = self.iparam_xmap
        if self.tmap > 0:
            ncfile['rparam_tmap'][:] = self.rparam_tmap
            ncfile['iparam_tmap'][:] = self.iparam_tmap
        ncfile.close()


    def rzmesh(self, iphi, iside):
        """
        Cross-section of mesh at toroidal index *iphi*

        iside = 0: flux tubes on right side (corresponds to coordinates (iphi, 0.0))
                1: flux tubes on left side (corresponds to coordinates (iphi-1, 1.0))
        """

        def transposed_arrays(*args):
            return [arr.T for arr in args]

        f2py.mmesh.import_unstructured_mmesh(self.symmetry, *transposed_arrays(
            self.phi, self.iphi_zone, self.x, self.g, self.b,
            self.izone_line + 1, self.inode_offset[self.nlines:], self.inode_offset[:self.nlines],
            self.bsect, self.corner, self.next_tube, self.izone_tube + 1,
            self.rparam_tmap, self.iparam_tmap, self.iparam_xmap
            ))

        rzmesh_dims = f2py.mmesh.mmesh_rzmesh(iphi, iside)
        return Rzmesh(iphi, iside, *transposed_arrays(*f2py.mmesh.get_rzmesh(*rzmesh_dims)))


    def inner_boundary(self, incr=1):
        """Export inner boundary as *Torosurf* object. Optional: *incr* > 1 for lower poloidal resolution."""

        # determine poloidal resolution along inner boundary
        for itube in range(self.ntubes):
            if self.next_tube[itube, 1, 1] == 0:
                nv = itube // incr + 1
                break

        # construct Torosurf representation of inner boundary
        inner_boundary = Torosurf.new(self.nphi, nv+1, abs(self.symmetry))
        inner_boundary.phi = self.phi
        for j in range(nv):
            k0 = self.inode_offset[self.corner[j*incr, 0]]
            k1 = k0 + self.nphi
            inner_boundary.rz[:, j, :] = self.x[k0:k1, :].T
        inner_boundary.rz[:, -1, :] = inner_boundary.rz[:, 0, :]
        return inner_boundary
