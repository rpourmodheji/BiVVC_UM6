class Getvolumes(object):
    """docstring for Getvolumes."""

    def __init__(self, arguments):
        self.volumes = arguments["Volumes"]

    def gettotalvolume(self):

        V_sa = self.volumes["sa"]
        V_sv = self.volumes["sv"]
        V_RA = self.volumes["RA"]
        V_RV = self.volumes["RV"]
        V_PA = self.volumes["PA"]
        V_pc = self.volumes["pc"]
        V_pv = self.volumes["pv"]
        V_LA = self.volumes["LA"]
        V_LV = self.volumes["LV"]
        V_AO = self.volumes["AO"]
        V_total = \
                V_sa +\
                V_sv +\
                V_RA +\
                V_RV +\
                V_PA +\
                V_pc +\
                V_pv +\
                V_LA +\
                V_LV +\
                V_AO
        return V_total
