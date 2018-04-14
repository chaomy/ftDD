#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-04-13 22:53:34
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-13 23:00:02

import plt_drv
import numpy as np
from optparse import OptionParser


class ftDDtool(plt_drv.plt_drv):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)

    def plot_velocity(self):
        self.set_keys()
        self.set_111plt()
        data = np.loadtxt('chk.txt')
        self.ax.plot(data[:, 0], data[:, 1])
        self.fig.savefig("FIG_VEL.png", **self.figsave)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = ftDDtool()
    dispatcher = {'pltv': drv.plot_velocity}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
