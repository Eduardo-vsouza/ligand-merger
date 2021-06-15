import os
import sys

import pandas as pd


class AffinityWrapper(object):
    def __init__(self, folder):
        self.folder = folder
        self.files = os.listdir(folder)
        self.topLigands = []

    def unify(self):
        modes = []
        affinities = []
        rmsd_1 = []
        rmsd_2 = []
        ligands = []
        for file in self.files:
            name = self.__get_names(file)
            start = 10000
            with open(f'{self.folder}/{file}', 'r') as aff:
                for i, line in enumerate(aff, 1):
                    line = line.rstrip()
                    if line.startswith('mode |'):
                        start = i
                    if i == start+3:
                        info = line.split("   ")
                        modes.append(info[1])
                        affinities.append(float(info[4]))
                        rmsd_1.append(info[6])
                        rmsd_2.append(info[8])
                        ligands.append(name)
                        self.topLigands.append(f'{name}{line}')
        df = pd.DataFrame(data={'ligand': ligands, 'mode': modes, 'affinity': affinities, 'rmsd l.b.': rmsd_1,
                                'rmsd u.b.': rmsd_2})
        df = df.sort_values(by='affinity')
        self.topLigands = df

    @staticmethod
    def __get_names(file):
        end = file.find("pdb")
        name = file[:end-1]
        return name

    def save(self, output='bestligands.txt'):
        self.topLigands.to_csv(output, sep='\t', index=False)
        return self


data = AffinityWrapper(folder=sys.argv[1])
tops = data.unify()
# print(tops)
if sys.argv[2] is not None:
    data.save(output=sys.argv[2])
else:
    data.save()