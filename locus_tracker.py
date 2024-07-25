class LocusTracker:
    def __init__(self, name, ch_count, res_count, both_count):
        self.name = name
        self.ch_count = ch_count
        self.res_count = res_count
        self.both_count = both_count

    def clone(self):
        return LocusTracker(self.name, self.ch_count, self.res_count, self.both_count)