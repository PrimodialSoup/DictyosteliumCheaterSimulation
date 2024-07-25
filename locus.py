class Locus:
    def __init__(self, cheat_v, res_v):
        self.cheater_allele = cheat_v
        self.resistor_allele = res_v

    def both(self):
        return int(self.cheater_allele + self.resistor_allele == 2) # 0 or 1

    def clone(self):
        return Locus(self.cheater_allele, self.resistor_allele)