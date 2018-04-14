#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-04-13 22:53:34
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-13 23:57:33

import plt_drv
import numpy as np
from optparse import OptionParser


class ftDDtool(plt_drv.plt_drv):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)
        self.set_keys()

    def loop_plt(self):
        for kk in ["Alue4", "Alue5", "Alue6"]:
            self.plot_check(kk)

    def plot_check(self, kk):
        data = np.loadtxt('chk' + kk + '.txt')
        self.set_111plt()
        self.ax.plot(data[:-2, 0], data[:-2, 1])
        self.fig.savefig("FIG_FRE_{}.png".format(kk), **self.figsave)

        # self.set_111plt()
        # self.ax.plot(data[:, 0], data[:, 2])
        # self.fig.savefig("FIG_VLE_{}.png".format(kk), **self.figsave)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = ftDDtool()
    dispatcher = {'plt': drv.loop_plt}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
