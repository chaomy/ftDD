#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-04-13 22:53:34
# @Last Modified by:   chaomy
# @Last Modified time: 2019-01-08 12:44:33

import plt_drv
import dd_data
import numpy as np
from itertools import cycle
from optparse import OptionParser
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt


class ftDDtool(plt_drv.plt_drv):

    def __init__(self):
        self.figsize1 = (10, 7)
        plt_drv.plt_drv.__init__(self)
        self.set_keys()

    def loop_plt(self):
        for kk in ["Alue4", "Alue5", "Alue6"]:
            self.plot_check(kk)

        # for kk in ["AluPrec1", "AluPrec2", "AluPrec3", "AluPrec4",
        #            "AluPrec5"]:
        #     self.plot_check(kk)

        # for kk in ["AluPrecA", "AluPrecB", "AluPrecC", "AluPrecD",
        #            "AluPrecE"]:
        #     self.plot_check(kk)

        # for kk in ["AluPrec5e4", "AluPrec5e6"]:
        #     self.plot_check(kk)

    def plot_fitv_total(self):
        self.set_611plt((13, 11))
        # self.ax.get_xaxis().get_major_formatter().set_useOffset(False)
        
        for kk, ax in zip(["Alue5", "AluPrec1", "AluPrec2",
                           "AluPrec3", "AluPrec4", "AluPrec5"], self.axls):
            data = np.loadtxt('chk' + kk + '.txt')
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0f'))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(5.0))
            ax.plot(1e6 * data[:, 0], data[:, -2],
                    label=dd_data.labtag[
                        kk + "_results"], color=next(self.coloriter),
                    lw=2.5)
            ax.plot(1e6 * data[:, 0], data[:, -1], lw=2.0, color='k')
            if (kk != "AluPrec5"):
                plt.setp(ax.get_xticklabels(), visible=False)

        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.add_title(r"Velocity [m/s] vs. Time [$\mu$s]", self.axls[0])
        # self.add_x_labels(cycle([r"Time [$\mu$s]"]), self.axls[-1])
        # self.add_y_labels(cycle(["Velcocity [m/s]"]), self.axls[2])
        self.fig.savefig("FIG_VLE_{}.png".format(kk), **self.figsave)

    def plot_check(self, kk):
        data = np.loadtxt('chk' + kk + '.txt')
        self.set_111plt(self.figsize1)
        # self.ax.get_xaxis().get_major_formatter().set_useOffset(False)
        self.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.0f'))
        # self.ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))

        self.ax.plot(1e6 * data[:, 0], data[:, -2], label="DD {}".format(kk))
        self.ax.plot(1e6 * data[:, 0], data[:, -1], label="FIT")

        self.add_x_labels(cycle([r"Time [$\mu$s]"]), *self.axls)
        self.add_y_labels(cycle(["Velcocity [m/s]"]), *self.axls)

        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig("FIG_VLE_{}.png".format(kk), **self.figsave)

    def plot_beta_den(self):
        for f, k in zip(["beta.tprec.txt", "beta.theta.txt"], ["Ttype", "Theta"]):
            data = np.loadtxt(f)
            self.set_111plt()
            self.ax.plot(1e6 * data[:, 1], data[:, 2],
                         label=k,  **next(self.keysiter))
            p = np.polyfit(data[:, 1], data[:, 2], 1)
            self.ax.plot(1e6 * data[:, 1], data[:, 1] * p[0],
                         label="k = {:.2e}".format(p[0]),  **next(self.keysiter))
            self.add_x_labels(cycle([r"Rho [1e-6/m^2]"]), *self.axls)
            self.add_y_labels(cycle(["Beta"]), *self.axls)
            self.add_legends(*self.axls)
            self.set_tick_size(*self.axls)
            self.fig.savefig("FIG_BETA_{}.png".format(k), **self.figsave)

    def plot_beta(self):
        data = np.loadtxt("beta.tprec.txt")
        self.set_111plt()
        self.ax.plot(data[:, 0], data[:, 2],
                     label="Ttype",  **next(self.keysiter))
        p = np.polyfit(data[:, 0], data[:, 2], 1)
        m = np.sqrt(data[:, 0] * p[0])
        self.ax.plot(data[:, 0], data[:, 0] * p[0],
                     label="k = {:.2e}".format(p[0]),  **next(self.keysiter))
        for i in range(5):
            fnm = "ft.t{:02d}".format(i + 1)
            with open(fnm, 'w') as fid:
                fid.write("ptype b\n")
                fid.write("burg 2.86e-10\n")
                fid.write("taylor 0.40824829046386307\n")
                fid.write("mu 2.7e10\n")
                fid.write("mob 2384.59\n")
                fid.write("alp 0.23709\n")
                fid.write("bta {}\n".format(m[i + 1]))

        data = np.loadtxt("beta.theta.txt")
        p = np.polyfit(data[:, 0], data[:, 2], 1)
        self.ax.plot(data[:, 0], data[:, 2],
                     label="Theta",  **next(self.keysiter))
        self.ax.plot(data[:, 0], data[:, 0] * p[0],
                     label="k = {:.2e}".format(p[0]),  **next(self.keysiter))
        m = np.sqrt(data[:, 0] * p[0])
        for i in range(5):
            fnm = "ft.a{:02d}".format(i + 1)
            with open(fnm, 'w') as fid:
                fid.write("ptype b\n")
                fid.write("burg 2.86e-10\n")
                fid.write("taylor 0.40824829046386307\n")
                fid.write("mu 2.7e10\n")
                fid.write("mob 2384.59\n")
                fid.write("alp 0.23709\n")
                fid.write("bta {}\n".format(m[i + 1]))

        self.add_x_labels(cycle([r"Volume fraction"]), *self.axls)
        self.add_y_labels(cycle(["Beta"]), *self.axls)
        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig("FIG_BETA_T.png", **self.figsave)

    def plot_beta_erate(self):
        data = np.loadtxt("beta.tprec.erate.txt")
        self.set_111plt()
        self.ax.plot(data[:, 0], data[:, 1],
                     label="Beta-Erate",  **next(self.keysiter))
        p = np.polyfit(data[:, 0], data[:, 1], 1)
        print(p)
        self.add_x_labels(cycle([r"Erate"]), *self.axls)
        self.add_y_labels(cycle(["Beta"]), *self.axls)
        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig("FIG_BETA_E.png", **self.figsave)

if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = ftDDtool()
    dispatcher = {'plt': drv.loop_plt,
                  'beta': drv.plot_beta,
                  'beta2': drv.plot_beta_erate,
                  'beta3': drv.plot_beta_den,
                  'fitall': drv.plot_fitv_total}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
