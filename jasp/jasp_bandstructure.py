from jasp import *
import os
import matplotlib.pyplot as plt
from ase.dft import DOS
def get_bandstructure(self,
                      kpts_path=None,
                      kpts_nintersections=10):

    kpts = [k[1] for k in kpts_path]
    labels = [k[0] for k in kpts_path]

    dos = DOS(self, width=0.2)
    d = dos.get_dos()
    e = dos.get_energies()

    ef = self.get_fermi_level()

    # run in non-selfconsistent directory
    cwd = os.getcwd()
    base, end = os.path.split(cwd)
    wd = cwd + '/bandstructure'
    self.clone(wd)

    with jasp(wd,
              kpts=kpts,
              kpts_nintersections=kpts_nintersections,
              reciprocal=True,
              nsw=0, # no ionic updates required
              isif=None,
              ibrion=None,
              debug=logging.DEBUG,
              icharg=11) as calc:

        calc.calculate()

        fig = plt.figure()
        with open('EIGENVAL') as f:
            line1 = f.readline()
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline()
            comment = f.readline()
            unknown, npoints, nbands = [int(x) for x in f.readline().split()]

            blankline = f.readline()

            band_energies = [[] for i in range(nbands)]

            for i in range(npoints):
                x,y,z, weight = [float(x) for x in f.readline().split()]

                for j in range(nbands):
                    fields = f.readline().split()
                    id, energy = int(fields[0]), float(fields[1])
                    band_energies[id-1].append(energy)
                blankline = f.readline()
            f.close()
            ax1 = plt.subplot(121)
            for i in range(nbands):
                plt.plot(range(npoints), np.array(band_energies[i]) - ef)

            ax = plt.gca()
            ax.set_xticks([]) # no tick marks
            plt.xlabel('k-vector')
            plt.ylabel('Energy (eV)')

            nticks = len(labels)/2 + 1
            ax.set_xticks(np.linspace(0,npoints,nticks))
            L = []
            L.append(labels[0])
            for i in range(2,len(labels)):
                if i % 2 == 0:
                    L.append(labels[i])
                else:
                    pass
            L.append(labels[-1])
            ax.set_xticklabels(L)
            plt.axhline(0,c='r')

    plt.subplot(122, sharey=ax1)
    plt.plot(d,e)
    plt.axhline(0, c='r')
    plt.ylabel('energy (eV)')
    plt.xlabel('DOS')

    plt.subplots_adjust(wspace=0.26)
    plt.show()
    return (npoints, band_energies, fig)

Vasp.get_bandstructure = get_bandstructure
