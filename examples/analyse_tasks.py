#!/usr/bin/env python
"""
Usage:
    analsyse_tasks.py [options] input.dat

where input.dat is a thread info file for a step.  Use the '-y interval' flag
of the swift command to create these.

The output is an analysis of the task timings, including deadtime per thread
and step, total amount of time spent for each task type, for the whole step
and per thread and the minimum and maximum times spent per task type.

This file is part of SWIFT.
Copyright (c) 2017 Peter W. Draper (p.w.draper@durham.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.collections as collections
import matplotlib.ticker as plticker
import pylab as pl
import sys
import argparse

#  Handle the command line.
parser = argparse.ArgumentParser(description="Analyse task dumps")

parser.add_argument("input", help="Thread data file (-y output)")
parser.add_argument("-v", "--verbose", dest="verbose",
                    help="Verbose output (default: False)",
                    default=False, action="store_true")

args = parser.parse_args()
infile = args.input

#  Tasks and subtypes. Indexed as in tasks.h.
TASKTYPES = ["none", "sort", "self", "pair", "sub_self", "sub_pair",
             "init_grav", "ghost", "extra_ghost", "drift_part",
             "drift_gpart", "kick1", "kick2", "timestep", "send", "recv",
             "grav_top_level", "grav_long_range", "grav_mm", "grav_down",
             "cooling", "sourceterms", "count"]

SUBTYPES = ["none", "density", "gradient", "force", "grav", "external_grav",
            "tend", "xv", "rho", "gpart", "multipole", "spart", "count"]

#  Read input.
data = pl.loadtxt( infile )

maxthread = int(max(data[:,0])) + 1
print "# Maximum thread id:", maxthread

#  Recover the start and end time
full_step = data[0,:]
tic_step = int(full_step[4])
toc_step = int(full_step[5])
CPU_CLOCK = float(full_step[-1]) / 1000.0
data = data[1:,:]
if args.verbose:
    print "CPU frequency:", CPU_CLOCK * 1000.0

#  Avoid start and end times of zero.
data = data[data[:,4] != 0]
data = data[data[:,5] != 0]

#  Calculate the time range.
total_t = (toc_step - tic_step)/ CPU_CLOCK
print "# Data range: ", total_t, "ms"

#  Correct times to relative values.
start_t = float(tic_step)
data[:,4] -= start_t
data[:,5] -= start_t

tasks = {}
tasks[-1] = []
for i in range(maxthread):
    tasks[i] = []

#  Gather into by thread data.
num_lines = pl.size(data) / 10
for line in range(num_lines):
    thread = int(data[line,0])
    tic = int(data[line,4]) / CPU_CLOCK
    toc = int(data[line,5]) / CPU_CLOCK
    tasktype = int(data[line,1])
    subtype = int(data[line,2])

    tasks[thread].append([tic,toc,tasktype,subtype])

#  Sort by tic and gather used thread ids.
threadids = []
for i in range(maxthread):
    if len(tasks[i]) > 0:
        tasks[i] = sorted(tasks[i], key=lambda task: task[0])
        threadids.append(i)

#  Times per task.
print "# Task times:"
print "# {0:<16s}: {1:>7s} {2:>9s} {3:>9s} {4:>9s} {5:>9s} {6:>9s}"\
      .format("type/subtype", "count","minimum", "maximum",
              "sum", "mean", "percent")
alltasktimes = {}
for i in threadids:
    tasktimes = {}
    for task in tasks[i]:
        key = TASKTYPES[task[2]] + "/" + SUBTYPES[task[3]]
        dt = task[1] - task[0]
        if not key in tasktimes:
            tasktimes[key] = []
        tasktimes[key].append(dt)

        if not key in alltasktimes:
            alltasktimes[key] = []
        alltasktimes[key].append(dt)

    print "# Thread : ", i
    for key in sorted(tasktimes.keys()):
        taskmin = min(tasktimes[key])
        taskmax = max(tasktimes[key])
        tasksum = sum(tasktimes[key])
        print "{0:18s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}"\
              .format(key, len(tasktimes[key]), taskmin, taskmax, tasksum,
                      tasksum / len(tasktimes[key]), tasksum / total_t * 100.0)
    print

print "# All threads : "
for key in sorted(alltasktimes.keys()):
    taskmin = min(alltasktimes[key])
    taskmax = max(alltasktimes[key])
    tasksum = sum(alltasktimes[key])
    print "{0:18s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}"\
          .format(key, len(alltasktimes[key]), taskmin, taskmax, tasksum,
                  tasksum / len(alltasktimes[key]),
                  tasksum / (len(threadids) * total_t) * 100.0)
print

#  Dead times.
print "# Deadtimes:"
print "# no.    : {0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:>9s} {5:>9s}"\
      .format("count", "minimum", "maximum", "sum", "mean", "percent")
alldeadtimes = []
for i in threadids:
    deadtimes = []
    last = 0
    for task in tasks[i]:
        dt = task[0] - last
        deadtimes.append(dt)
        last = task[1]
    dt = total_t - last
    deadtimes.append(dt)

    deadmin = min(deadtimes)
    deadmax = max(deadtimes)
    deadsum = sum(deadtimes)
    print "thread {0:2d}: {1:9d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}"\
          .format(i, len(deadtimes), deadmin, deadmax, deadsum,
                  deadsum / len(deadtimes), deadsum / total_t * 100.0)
    alldeadtimes.extend(deadtimes)

deadmin = min(alldeadtimes)
deadmax = max(alldeadtimes)
deadsum = sum(alldeadtimes)
print "all      : {0:9d} {1:9.4f} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.2f}"\
      .format(len(alldeadtimes), deadmin, deadmax, deadsum,
              deadsum / len(alldeadtimes),
              deadsum / (len(threadids) * total_t ) * 100.0)
print


sys.exit(0)
