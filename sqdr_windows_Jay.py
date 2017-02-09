#!/usr/bin/env python

import sys
import math
import numpy as np
import os
import argparse
import time
# python 2.7, 3.x compatible, maybe
try:
    import tkinter as tk
    from tkinter import filedialog
except:
    import Tkinter as tk
    import tkFileDialog as filedialog
import matplotlib
from matplotlib import pyplot, ticker
plt = pyplot
from scipy.optimize import curve_fit
from IPython import embed
import threading

# factor to convert from scaled voltage x3 to EMU
# may want to check this later
squid_factor = 1.09589

# limits for rso_response fitting
x1_min = None
x1_max = None
x2_min = None
x2_max = None
x3_min = None
x3_max = None
x4_min = None
x4_max = None

def within_lims(x1, x2, x3, x4):
    if x1_min is not None:
        if x1 < x1_min:
            return False
    if x2_min is not None:
        if x2 < x2_min:
            return False
    if x3_min is not None:
        if x3 < x3_min:
            return False
    if x4_min is not None:
        if x4 < x4_min:
            return False
    if x1_max is not None:
        if x1 > x1_max:
            return False
    if x2_max is not None:
        if x2 > x2_max:
            return False
    if x3_max is not None:
        if x3 > x3_max:
            return False
    if x4_max is not None:
        if x4 > x4_max:
            return False
    return True

# column names in raw squid file
class colnames:
    time = 'Time'
    start_temp = 'Start_Temperature_K'
    end_temp = 'End_Temperature_K'
    field = 'Field_Oe'
    pos = 'Position_cm'
    scaled_vol = 'Long_Scaled_Response'
    raw_vol = 'Long_Voltage'
    demeaned_vol = 'Long_Demeaned_Voltage'
    demeaned_vol_fit = 'Long_Demeaned_Fit'

# parse and read in raw squid data
def read_sqd(fname):
    sqd_data = np.genfromtxt(fname, delimiter = ',',
            names = True, skip_header = 30)
    # ignore empty columns
    columnnames = [n for n in sqd_data.dtype.names if not np.isnan(sqd_data[n][0])]
    return sqd_data[columnnames]

# separate complete raw sqd array into separate measurements
def split_sqd(sqd_data):
    unique_times = sorted(set(sqd_data[colnames.time]))
    measurements = {}
    for i, d in enumerate(sqd_data):
        if not d[colnames.time] in measurements:
            measurements[d[colnames.time]] = [i]
        else:
            measurements[d[colnames.time]].append(i)
    return [sqd_data[index] for time, index in sorted(measurements.items())]

# use scientific format for graph y axis
def set_sci_format(*fignums):
    if fignums == []:
        fignums = plt.get_fignums()
    for i in fignums:
        f = plt.figure(i)
        formatter = ticker.ScalarFormatter(useMathText = True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-3, 3))
        for ax in f.get_axes():
            if ax.get_yscale() == 'linear':
                ax.yaxis.set_major_formatter(formatter)

# window for controlling stuff
class Window(threading.Thread):
    def __init__(self,):
        threading.Thread.__init__(self)
        self.start()

    def callback(self):
        self.root.quit()

    def run(self):
        self.root = tk.Tk()
        self.root.protocol("WM_DELETE_WINDOW", self.callback)
        self.create_widgets()
        time.sleep(0.5)
        self.root.mainloop()

    def load_squid(self, squid):
        self.squid = squid
        self.total_scans["text"] = "/ {}".format(len(self.squid.scans) - 1)
        self.scan_selector(0)

    def create_widgets(self):
       # Frame for current scan info
        self.current_scan = tk.LabelFrame(self.root, text = "Current scan", labelanchor = "n")

        self.scan_select = tk.Frame(self.current_scan)
        self.scan_entry = tk.Entry(self.scan_select, width = 4)
        self.scan_entry_contents = tk.StringVar()
        self.scan_entry_contents.set("0")
        self.scan_entry["textvariable"] = self.scan_entry_contents
        get_entry = lambda event : self.scan_selector(int(self.scan_entry_contents.get()))
        self.scan_entry.bind('<Return>', get_entry)
        self.decrement = tk.Button(self.scan_select, text = "<", command = self.scan_decrement)
        self.total_scans = tk.Label(self.scan_select, text = "")
        self.increment = tk.Button(self.scan_select, text = ">", command = self.scan_increment)

        self.decrement.pack(side = "left", padx = 5)
        self.scan_entry.pack(side = "left", padx =5)
        self.total_scans.pack(side = "left", padx =5)
        self.increment.pack(side = "left", padx =5)
        self.scan_select.pack(pady = 5)

        self.scan_temperature = tk.Frame(self.current_scan)
        self.temp_label = tk.Label(self.scan_temperature, text = "Temperature:", anchor = 'w')
        self.temp_value = tk.Label(self.scan_temperature, text = "", anchor = 'e')
        self.temp_label.pack(side = "left")
        self.temp_value.pack(side = "right")
        self.scan_temperature.pack(fill = 'both', padx = 5, pady = 5)

        self.scan_field = tk.Frame(self.current_scan)
        self.field_label = tk.Label(self.scan_field, text = "Applied field:")
        self.field_value = tk.Label(self.scan_field, text = "")
        self.field_label.pack(side = "left")
        self.field_value.pack(side = "right")
        self.scan_field.pack(fill = "both", padx = 5, pady = 5)

        self.scan_moment = tk.Frame(self.current_scan)
        self.moment_label = tk.Label(self.scan_moment, text = "Calculated moment:")
        self.moment_value = tk.Label(self.scan_moment, text = "")
        self.moment_label.pack(side = "left")
        self.moment_value.pack(side = "right")
        self.scan_moment.pack(fill = "both", padx = 5, pady = 5)

        def rs():
            self.scan.fit(center = True, print_new_fit = False)
            self.update_scan_details()
            self.scan.plot_fit()
        def ras():
            self.squid.refit_centers()
            self.update_scan_details()
            self.scan.plot_fit()
            self.squid.plot_dependence()
            self.squid.plot_params()
        def rss():
            self.squid.selected_scan.original_fit(center = True,)
            self.update_scan_details()
            self.squid.selected_scan.plot_fit()
        def rsa():
            for scan in self.squid.scans:
                scan.original_fit(center = True,)
            self.update_scan_details()
            self.scan.plot_fit()

        self.refitting = tk.LabelFrame(self.current_scan, text = 'Refitting', labelanchor = "n")

        self.refit_scan = tk.Button(self.refitting, text = "Refit center of\ncurrent scan", command = rs)
        self.refit_all_scans = tk.Button(self.refitting, text = "Refit center of\nall scans", command = ras)
        self.refit_squid_scan = tk.Button(self.refitting, text = "Recalculate SQUID's fit\nof current scan", command = rss)
        self.refit_squid_all = tk.Button(self.refitting, text = "Recalculate SQUID's fit\nof all scans", command = rsa)

        self.refit_scan.grid(row = 0, column = 0, padx = 5, pady = 5)
        self.refit_all_scans.grid(row = 1, column = 0, padx = 5, pady = 5)
        self.refit_squid_scan.grid(row = 0, column = 1, padx = 5, pady = 5)
        self.refit_squid_all.grid(row = 1, column = 1, padx = 5, pady = 5)

        self.refitting.pack(padx = 5, pady = 5)

        self.current_scan.grid(row = 0, column = 0, sticky = "wens", padx = 5, pady = 5)

       # Frame for fit details
        self.fit_details = tk.LabelFrame(self.root, text = "Fit details", labelanchor = "n")
        self.x1 = tk.Label(self.fit_details, text = "x1")
        self.x2 = tk.Label(self.fit_details, text = "x2")
        self.x3 = tk.Label(self.fit_details, text = "x3")
        self.x4 = tk.Label(self.fit_details, text = "x4")
        self.chi2 = tk.Label(self.fit_details, text = "chi2")
        
        self.best_column = tk.Label(self.fit_details, text = "Best fit")
        self.x1_best_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x2_best_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x3_best_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x4_best_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.chi2_best_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x1_best = tk.Label(self.x1_best_f, text = "", width = 9)
        self.x2_best = tk.Label(self.x2_best_f, text = "", width = 9)
        self.x3_best = tk.Label(self.x3_best_f, text = "", width = 9)
        self.x4_best = tk.Label(self.x4_best_f, text = "", width = 9)
        self.chi2_best = tk.Label(self.chi2_best_f, text = "", width = 9)
        self.x1_best.pack()
        self.x2_best.pack()
        self.x3_best.pack()
        self.x4_best.pack()
        self.chi2_best.pack()

        self.squid_column = tk.Label(self.fit_details, text = "Squid fit")
        self.x1_squid_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x2_squid_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x3_squid_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x4_squid_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.chi2_squid_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x1_squid = tk.Label(self.x1_squid_f, text = "", width = 9)
        self.x2_squid = tk.Label(self.x2_squid_f, text = "", width = 9)
        self.x3_squid = tk.Label(self.x3_squid_f, text = "", width = 9)
        self.x4_squid = tk.Label(self.x4_squid_f, text = "", width = 9)
        self.chi2_squid = tk.Label(self.chi2_squid_f, text = "", width = 9)
        self.x1_squid.pack()
        self.x2_squid.pack()
        self.x3_squid.pack()
        self.x4_squid.pack()
        self.chi2_squid.pack()

        self.fit_column = tk.Label(self.fit_details, text = "Initial fitting\nparameters")
        self.x1_fit = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x2_fit = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x3_fit = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x4_fit = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x1_fit_contents = tk.StringVar()
        self.x2_fit_contents = tk.StringVar()
        self.x3_fit_contents = tk.StringVar()
        self.x4_fit_contents = tk.StringVar()
        self.x1_fit_contents.set("0")
        self.x2_fit_contents.set("0")
        self.x3_fit_contents.set("0")
        self.x4_fit_contents.set("0")
        self.x1_fit["textvariable"] = self.x1_fit_contents
        self.x2_fit["textvariable"] = self.x2_fit_contents
        self.x3_fit["textvariable"] = self.x3_fit_contents
        self.x4_fit["textvariable"] = self.x4_fit_contents

        self.min_column = tk.Label(self.fit_details, text = 'Lower limit')
        self.x1_min = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x2_min = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x3_min = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x4_min = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x1_min_contents = tk.StringVar()
        self.x2_min_contents = tk.StringVar()
        self.x3_min_contents = tk.StringVar()
        self.x4_min_contents = tk.StringVar()
        self.x1_min_contents.set("")
        self.x2_min_contents.set("")
        self.x3_min_contents.set("")
        self.x4_min_contents.set("")
        self.x1_min["textvariable"] = self.x1_min_contents
        self.x2_min["textvariable"] = self.x2_min_contents
        self.x3_min["textvariable"] = self.x3_min_contents
        self.x4_min["textvariable"] = self.x4_min_contents

        self.max_column = tk.Label(self.fit_details, text = 'Upper limit')
        self.x1_max = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x2_max = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x3_max = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x4_max = tk.Entry(self.fit_details, width = 10, relief = 'sunken',)
        self.x1_max_contents = tk.StringVar()
        self.x2_max_contents = tk.StringVar()
        self.x3_max_contents = tk.StringVar()
        self.x4_max_contents = tk.StringVar()
        self.x1_max_contents.set("")
        self.x2_max_contents.set("")
        self.x3_max_contents.set("")
        self.x4_max_contents.set("")
        self.x1_max["textvariable"] = self.x1_max_contents
        self.x2_max["textvariable"] = self.x2_max_contents
        self.x3_max["textvariable"] = self.x3_max_contents
        self.x4_max["textvariable"] = self.x4_max_contents

        self.opt_column = tk.Label(self.fit_details, text = "Optimised fit")
        self.x1_opt_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x2_opt_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x3_opt_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x4_opt_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.chi2_opt_f = tk.LabelFrame(self.fit_details, relief = 'sunken', padx = 5)
        self.x1_opt = tk.Label(self.x1_opt_f, text = "", width = 9)
        self.x2_opt = tk.Label(self.x2_opt_f, text = "", width = 9)
        self.x3_opt = tk.Label(self.x3_opt_f, text = "", width = 9)
        self.x4_opt = tk.Label(self.x4_opt_f, text = "", width = 9)
        self.chi2_opt = tk.Label(self.chi2_opt_f, text = "", width = 9)
        self.x1_opt.pack()
        self.x2_opt.pack()
        self.x3_opt.pack()
        self.x4_opt.pack()
        self.chi2_opt.pack()


        cbf = lambda : self.insert_initial_fit('best')
        csf = lambda : self.insert_initial_fit('squid')
        
        self.copy_best_fit = tk.Button(self.fit_details, text = 'Copy params', command = cbf)
        self.copy_squid_fit = tk.Button(self.fit_details, text = 'Copy params', command = csf)
        self.optimise_params = tk.Button(self.fit_details, text = 'Optimise', command = self.optimise)
        self.apply_lims = tk.Button(self.fit_details, text = "Apply fitting\nlimits", command = self.update_fitting_lims)
        self.clear_lims = tk.Button(self.fit_details, text = "Clear all\n limits", command = self.clear_fitting_lims)

        self.x1.grid(row = 1, column = 0, padx = 5, pady = 5)
        self.x2.grid(row = 2, column = 0, padx = 5, pady = 5)
        self.x3.grid(row = 3, column = 0, padx = 5, pady = 5)
        self.x4.grid(row = 4, column = 0, padx = 5, pady = 5)
        self.chi2.grid(row = 5, column = 0, padx = 5, pady = 5)

        self.best_column.grid(row = 0, column = 1, padx = 5, pady = 5)
        self.x1_best_f.grid(row = 1, column = 1, padx = 5, pady = 5)
        self.x2_best_f.grid(row = 2, column = 1, padx = 5, pady = 5)
        self.x3_best_f.grid(row = 3, column = 1, padx = 5, pady = 5)
        self.x4_best_f.grid(row = 4, column = 1, padx = 5, pady = 5)
        self.chi2_best_f.grid(row = 5, column = 1, padx = 5, pady = 5)
        self.copy_best_fit.grid(row = 6, column = 1, padx = 5, pady = 5)

        self.squid_column.grid(row = 0, column = 2, padx = 5, pady = 5)
        self.x1_squid_f.grid(row = 1, column = 2, padx = 5, pady = 5)
        self.x2_squid_f.grid(row = 2, column = 2, padx = 5, pady = 5)
        self.x3_squid_f.grid(row = 3, column = 2, padx = 5, pady = 5)
        self.x4_squid_f.grid(row = 4, column = 2, padx = 5, pady = 5)
        self.chi2_squid_f.grid(row = 5, column = 2, padx = 5, pady = 5)
        self.copy_squid_fit.grid(row = 6, column = 2, padx = 5, pady = 5)

        self.fit_column.grid(row = 0, column = 3, padx = 5, pady = 5)
        self.x1_fit.grid(row = 1, column = 3, padx = 5, pady = 5)
        self.x2_fit.grid(row = 2, column = 3, padx = 5, pady = 5)
        self.x3_fit.grid(row = 3, column = 3, padx = 5, pady = 5)
        self.x4_fit.grid(row = 4, column = 3, padx = 5, pady = 5)
        self.optimise_params.grid(row = 6, column = 3, padx = 5, pady = 5)

        self.min_column.grid(row = 0, column = 4, padx = 5, pady = 5)
        self.x1_min.grid(row = 1, column = 4, padx = 5, pady = 5)
        self.x2_min.grid(row = 2, column = 4, padx = 5, pady = 5)
        self.x3_min.grid(row = 3, column = 4, padx = 5, pady = 5)
        self.x4_min.grid(row = 4, column = 4, padx = 5, pady = 5)
        self.apply_lims.grid(row = 6, column = 4, padx = 5, pady = 5)

        self.max_column.grid(row = 0, column = 5, padx = 5, pady = 5)
        self.x1_max.grid(row = 1, column = 5, padx = 5, pady = 5)
        self.x2_max.grid(row = 2, column = 5, padx = 5, pady = 5)
        self.x3_max.grid(row = 3, column = 5, padx = 5, pady = 5)
        self.x4_max.grid(row = 4, column = 5, padx = 5, pady = 5)
        self.clear_lims.grid(row = 6, column = 5, padx = 5, pady = 5)

        self.opt_column.grid(row = 0, column = 6, padx = 5, pady = 5)
        self.x1_opt_f.grid(row = 1, column = 6, padx = 5, pady = 5)
        self.x2_opt_f.grid(row = 2, column = 6, padx = 5, pady = 5)
        self.x3_opt_f.grid(row = 3, column = 6, padx = 5, pady = 5)
        self.x4_opt_f.grid(row = 4, column = 6, padx = 5, pady = 5)
        self.chi2_opt_f.grid(row = 5, column = 6, padx = 5, pady = 5)


        self.fit_details.grid(row = 1, column = 0, columnspan = 2, sticky = 'wens', padx = 5, pady = 5)

       # Frame for general commands
        self.gen_commands = tk.LabelFrame(self.root, text = "Commands", labelanchor = "n")
        self.refit_all_button = tk.Button(self.gen_commands, text = "Refit all scans")
        self.save_to = tk.Button(self.gen_commands, text = "Save as...", command = self.save_fits)
        self.dependence = tk.LabelFrame(self.gen_commands, text = "Independent variable", labelanchor = 'n')

       # switch dependence
        stt = lambda : self.switch_dependence('temperature')
        stf = lambda : self.switch_dependence('field')

        self.temp_dependence = tk.Button(self.dependence, text = "Temperature", relief = 'sunken', command = stt)
        self.field_dependence = tk.Button(self.dependence, text = "Field", command = stf)

        self.temp_dependence.pack(side = "left", padx = 5, pady = 5)
        self.field_dependence.pack(side = "right", padx = 5, pady = 5)
        self.replot = tk.LabelFrame(self.gen_commands, text = "Redraw plot", labelanchor = "n")

        dplot = lambda : self.squid.plot_dependence()
        dparam = lambda : self.squid.plot_params()
        iscan = lambda : self.squid.selected_scan.plot_fit()
        gp = lambda : self.squid.selected_scan.update_offset(plot_only = True)

        self.dependence_plot = tk.Button(self.replot, text = "Moment dependence", command = dplot)
        self.dependence_params = tk.Button(self.replot, text = "Fitting parameters", command = dparam)
        self.individual_scan = tk.Button(self.replot, text = "Selected scan", command = iscan)
        self.gradient_plot = tk.Button(self.replot, text = "Scan gradient", command = gp)
        self.save_to.pack()
        self.dependence.pack(padx = 5, pady = 5)
        self.dependence_plot.pack(padx = 5, pady = 5, fill = 'both')
        self.dependence_params.pack(padx = 5, pady = 5, fill = 'both')
        self.individual_scan.pack(padx = 5, pady = 5, fill = 'both')
        self.gradient_plot.pack(padx = 5, pady = 5, fill = 'both')
        self.replot.pack(padx = 5, pady = 5)
        self.gen_commands.grid(row = 0, column = 1, sticky = 'wens', padx = 5, pady = 5)

    def update_fitting_lims(self):
        global x1_min, x2_min, x3_min, x4_min
        global x1_max, x2_max, x3_max, x4_max
        try:
            x1_min = float(self.x1_min.get())
        except ValueError:
            x1_min = None
        try:
            x2_min = float(self.x2_min.get())
        except ValueError:
            x2_min = None
        try:
            x3_min = float(self.x3_min.get())
        except ValueError:
            x3_min = None
        try:
            x4_min = float(self.x4_min.get())
        except ValueError:
            x4_min = None
        try:
            x1_max = float(self.x1_max.get())
        except ValueError:
            x1_max = None
        try:
            x2_max = float(self.x2_max.get())
        except ValueError:
            x2_max = None
        try:
            x3_max = float(self.x3_max.get())
        except ValueError:
            x3_max = None
        try:
            x4_max = float(self.x4_max.get())
        except ValueError:
            x4_max = None
        print("Fitting limits applied\n"
                "x1: {}, {}\n"
                "x2: {}, {}\n"
                "x3: {}, {}\n"
                "x4: {}, {}".format(x1_min, x1_max, x2_min, x2_max,
                    x3_min, x3_max, x4_min, x4_max))

    def clear_fitting_lims(self):
        for x in [self.x1_min_contents, self.x2_min_contents, self.x3_min_contents, self.x4_min_contents,
                self.x1_max_contents, self.x2_max_contents, self.x3_max_contents, self.x4_max_contents]:
            x.set("")
        self.update_fitting_lims()

    def insert_initial_fit(self, type = None):
        if type is None:
            return
        scan = self.squid.selected_scan
        if type == 'best':
            x1, x2, x3, x4 = ['{:.3e}'.format(x) for x in scan.best_popt]
        elif type == 'squid':
            x1, x2, x3, x4 = ['{:.3e}'.format(x) for x in scan.squid_popt]
        else:
            # something went wrong
            return
        self.x1_fit_contents.set(x1)
        self.x2_fit_contents.set(x2)
        self.x3_fit_contents.set(x3)
        self.x4_fit_contents.set(x4)

    def get_scan_entry(self):
        try:
            return int(self.scan_entry_contents.get())
        except ValueError:
            print("Input index must be an integer")

    def scan_decrement(self):
        if self.get_scan_entry() < 1:
            return
        new_scan_index = self.get_scan_entry() - 1
        self.update_scan_details(new_scan_index)
        self.scan_selector(new_scan_index)

    def scan_increment(self):
        if self.get_scan_entry() >= len(self.squid.scans) - 1:
            return
        new_scan_index = self.get_scan_entry() + 1
        self.update_scan_details(new_scan_index)
        self.scan_selector(new_scan_index)

    def scan_selector(self, new_scan_index):
        if new_scan_index < 0 or new_scan_index >= len(self.squid.scans):
            return
        self.scan_entry_contents.set(str(new_scan_index))
        self.update_scan_details(new_scan_index)
        self.squid.select_scan(new_scan_index)
        self.scan = self.squid.selected_scan

    def update_scan_details(self, scan = None):
        if scan is None:
            scan = self.scan
        if type(scan) == int:
            scan = self.squid[scan]
        self.temp_value["text"] = '{:.3f} K'.format(scan.temperature)
        self.field_value["text"] = '{:.3f} Oe'.format(scan.field)
        self.moment_value["text"] = '{:.3e} EMU'.format(scan.get_best_moment())
        self.x1_best["text"] = '{:.3e}'.format(scan.best_popt[0])
        self.x2_best["text"] = '{:.3e}'.format(scan.best_popt[1])
        self.x3_best["text"] = '{:.3e}'.format(scan.best_popt[2])
        self.x4_best["text"] = '{:.3e}'.format(scan.best_popt[3])
        self.chi2_best["text"] = '{:.3e}'.format(scan.best_chi2)
        self.x1_squid["text"] = '{:.3e}'.format(scan.squid_popt[0])
        self.x2_squid["text"] = '{:.3e}'.format(scan.squid_popt[1])
        self.x3_squid["text"] = '{:.3e}'.format(scan.squid_popt[2])
        self.x4_squid["text"] = '{:.3e}'.format(scan.squid_popt[3])
        self.chi2_squid["text"] = '{:.3e}'.format(scan.squid_chi2)

    def optimise(self):
        try:
            x1 = float(self.x1_fit_contents.get())
            x2 = float(self.x2_fit_contents.get())
            x3 = float(self.x3_fit_contents.get())
            x4 = float(self.x4_fit_contents.get())
        except ValueError:
            print("Fitting parameters must be floats")
            return
        p0 = [x1, x2, x3, x4]
        scan = self.squid.selected_scan
        popt, pcov = curve_fit(rso_response, scan.position, scan.voltage, p0 = p0)
        self.x1_opt["text"] = '{:.3e}'.format(popt[0])
        self.x2_opt["text"] = '{:.3e}'.format(popt[1])
        self.x3_opt["text"] = '{:.3e}'.format(popt[2])
        self.x4_opt["text"] = '{:.3e}'.format(popt[3])

        sq_res_sum = sum((rso_response(scan.position, *popt, ignore_lims = True) - scan.voltage)**2)
        chi2 = sq_res_sum / (len(scan.position) - len(popt))
        self.chi2_opt["text"] = '{:.3e}'.format(chi2)
        ax1, ax2, ax3, ax4 = self.squid.rmf.get_axes()
        plt.sca(ax1)
        fitted_voltage = rso_response(scan.position, *popt, ignore_lims = True)
        plt.plot(scan.position, fitted_voltage, linestyle = '--', marker = '^')
        plt.sca(ax3)
        plt.plot(scan.position, fitted_voltage - scan.voltage, linestyle = '--', marker = '^')
        self.squid.rmf.canvas.draw()
        if chi2 < scan.best_chi2:
            scan.best_popt, scan.best_chi2 = popt, chi2
            self.update_scan_details(scan)

    def save_fits(self):
        dump_fit(self.squid)

    def switch_dependence(self, dependent):
        if dependent == 'temperature':
            self.temp_dependence["relief"] = 'sunken'
            self.field_dependence["relief"] = 'raised'
        else:
            self.temp_dependence["relief"] = 'raised'
            self.field_dependence["relief"] = 'sunken'
        if dependent == self.squid.dependent:
            return
        else:
            self.squid.toggle_dependence()

# class to contain raw squid data, methods...
class Squid:
    def __init__(self, fname, dependent = 'temperature'):
        self.fname = fname
        self.autoupdate_on_click = True
        self.dependent = dependent
        start_time = time.time()
        self.raw_data = read_sqd(fname)
        self.scans = []
        split_data = split_sqd(self.raw_data)
        # start fitting from coldest temperature; then sort by time of scan
        # if dependent variable is temperature
        if dependent == 'temperature':
            split_data.sort(key = lambda x : x[colnames.start_temp][0])
        self.scans.append(Squid_measurement(split_data[0], parent = self, dependent = dependent))
        for i, s in enumerate(split_data):
            if not i:
                continue
            self.scans.append(Squid_measurement(s, parent = self,
                dependent = dependent, p0 = self.scans[-1].best_popt))
        # clicked_scans contain events
        self.clicked_scans = []
        for s in self.scans:
            s.original_fit()
        self.figures(1, 2, 3, 4)
        if dependent == 'temperature':
            self.scans.sort(key = lambda x : x.time)
        self.select_scan(self.scans[0])
        self.plot_dependence()
        
        print("\n{:.3f} seconds to load {}".format(time.time() - start_time, fname))

    def __getitem__(self, k):
        return self.scans[k]

    def __len__(self):
        return len(self.scans)

    def toggle_dependence(self):
        self.dependent = {'temperature':'field', 'field':'temperature'}[self.dependent]
        for scan in self.scans:
            scan.dependent = {'temperature':scan.temperature, 'field':scan.field}[self.dependent]
        self.plot_dependence()
        self.plot_params()

    def toggle_autoupdate(self):
        self.autoupdate_on_click = not self.autoupdate_on_click

    def figures(self, *args):
        if args == ():
            args = (1,2,3,4)
        # functions for selecting scans off plots
        pe = lambda x : self.get_scans(x)
        bre = lambda x : self.closest_scan(x)
        # function for refitting scan by clicking starting center
        cfe = lambda x : self.click_fit(x)

        # figure for temperature (field) dependence (1)
        if 1 in args:
            self.temp_dependence_fig = plt.subplots(3, sharex = True)
            self.tdf = self.temp_dependence_fig[0]
            title = "{} dependence".format('Temperature' if self.dependent == 'temperature' else 'Field')
            plt.suptitle(title)
            self.tdf.canvas.mpl_connect('pick_event', pe)
            self.tdf.canvas.mpl_connect('button_release_event', bre)

        # figure to show fit parameters (temperature dependence)
        if 2 in args:
            self.param_fig = plt.subplots(4, 2, sharex = True, sharey = 'row')
            self.pf = self.param_fig[0]
            plt.suptitle('Fit parameters')
            self.plot_params()
            self.pf.canvas.mpl_connect('pick_event', pe)
            self.pf.canvas.mpl_connect('button_release_event', bre)
            set_sci_format(self.pf.number)

        # figure for offset correction
        if 3 in args:
            self.offset_fig = plt.subplots(2, sharex = True)
            self.of = self.offset_fig[0]
            plt.suptitle('Slope/offset correction')
        
        # figure for individual raw measurement data
        if 4 in args:
            self.raw_measurement_fig = plt.subplots(2, 2, sharex = 'col', sharey = 'row')
            self.rmf = self.raw_measurement_fig[0]
            plt.suptitle('Raw measurement')
            # fit by clicking
            # self.rmf.canvas.mpl_connect('button_release_event', cfe)

        set_sci_format()

    def get_scans(self, event):
        self.clicked_scans.append(event)

    def closest_scan(self, event):
        if not self.clicked_scans:
            return
        click_x, click_y = event.x, event.y
        dists = []
        for s in self.clicked_scans:
            xy = s.artist.get_xydata()
            pt_x, pt_y = s.artist.get_axes().transData.transform(xy)[0]
            dists.append(math.hypot(click_x - pt_x, click_y - pt_y))
        closest = self.clicked_scans[dists.index(min(dists))]
        self.clicked_scans = []
        self.select_scan(closest.artist.scan)

    def click_fit(self, event):
        scan = self.selected_scan
        if event.xdata is None:
            return
        center = -event.xdata
        v_range = max(scan.voltage) - min(scan.voltage)
        p0 = [0, 0, v_range/2.65, center] 
        popt, pcov = curve_fit(rso_response, scan.position, scan.voltage, p0 = p0)
        sq_res_sum = sum((rso_response(scan.position, *popt, ignore_lims = True) - scan.voltage)**2)
        chi2 = sq_res_sum / (len(scan.position) - len(popt))
        ax1, ax2, ax3, ax4 = self.rmf.get_axes()
        plt.sca(ax1)
        fitted_voltage = rso_response(scan.position, *popt, ignore_lims = True)
        plt.plot(scan.position, fitted_voltage, linestyle = '--', marker = '^')
        plt.sca(ax3)
        plt.plot(scan.position, fitted_voltage - scan.voltage, linestyle = '--', marker = '^')
        if chi2 < scan.best_chi2:
            print('Better fit found!')
            if not self.autoupdate_on_click:
                print('Not updated')
            else:
                scan.best_popt, scan.best_chi2 = popt, chi2

    def select_scan(self, scan):
        if type(scan) == int:
            scan = self.scans[scan]
        self.selected_scan = scan
        try:
            if plt.fignum_exists(self.parent.of.number):
                scan.update_offset(plot_only = True)
        except:
            pass
        scan.plot_fit()
        scan.update_offset(plot_only = True)
        scan_index = self.scans.index(self.selected_scan)
        new_title = 'Scan index {0}, {1} = {2:.2f}'
        title = new_title.format(scan_index, self.dependent, scan.dependent)
        self.rmf.texts[0].set_text(title)
        self.rmf.canvas.draw()
        self.of.texts[0].set_text(title)
        self.of.canvas.draw()

    def plot_dependence(self, chi2bound = 100):
        try:
            ax1, ax2, ax3 = self.tdf.get_axes()
        except:
            self.figures(1)
            ax1, ax2, ax3 = self.tdf.get_axes()
        selected_scans = [s for s in self.scans if s.best_chi2 < chi2bound]
        title = "{} dependence".format('Temperature' if self.dependent == 'temperature' else 'Field')
        self.tdf.texts[0].set_text(title)
        # ax1 for temp/field dependence of magnetization
        plt.sca(ax1)
        plt.cla()
        ax1.set_title('Magnetic moment dependence')
        ax1.set_ylabel('Moment (EMU)')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.best_popt[2], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
            squid_pt, = plt.plot(s.dependent, s.squid_popt[2], marker = 'o', color = 'g', picker = 10)
            squid_pt.scan = s
        plt.setp(ax1.get_xticklabels(), visible = False)
        # ax2 for inverse susceptibility, units 1/cm^3
        plt.sca(ax2)
        plt.cla()
        ax2.set_title('Inverse susceptibility')
        ax2.set_ylabel(r'1/$\chi$ (1/cm$^3$)')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.field/s.best_popt[2], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
            squid_pt, = plt.plot(s.dependent, s.field/s.squid_popt[2], marker = 'o', color = 'g', picker = 10)
            squid_pt.scan = s
        plt.setp(ax2.get_xticklabels(), visible = False)
        # ax3 for fit parameter
        plt.sca(ax3)
        plt.cla()
        ax3.set_title('Chi2 fitting parameter')
        ax3.set_ylabel('chi2')
        if self.dependent == 'temperature':
            ax3.set_xlabel('Temperature (K)')
        else:
            ax3.set_xlabel('Applied field (Oe)')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.best_chi2, marker = 'o', color = 'b', picker = 10)
            pt.scan = s
            squid_pt, = plt.plot(s.dependent, s.squid_chi2, marker = 'o', color = 'g', picker = 10)
            squid_pt.scan = s
        ax3.set_yscale('log')
        set_sci_format(self.tdf.number)
        self.tdf.canvas.draw()

    def plot_params(self, chi2bound = 100):
        try:
            ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8 = self.pf.get_axes()
        except:
            self.figures(2)
            ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8 = self.pf.get_axes()
        selected_scans = [s for s in self.scans if s.best_chi2 < chi2bound]
        plt.sca(ax1)
        plt.cla()
        ax1.set_ylabel('x1')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.best_popt[0], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.setp(ax1.get_xticklabels(), visible = False)
        plt.sca(ax2)
        plt.cla()
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.squid_popt[0], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.setp(ax2.get_xticklabels(), visible = False)
        plt.setp(ax2.get_yticklabels(), visible = False)
        plt.sca(ax3)
        plt.cla()
        ax3.set_ylabel('x2')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.best_popt[1], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.setp(ax3.get_xticklabels(), visible = False)
        plt.sca(ax4)
        plt.cla()
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.squid_popt[1], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.setp(ax4.get_xticklabels(), visible = False)
        plt.setp(ax4.get_yticklabels(), visible = False)
        plt.sca(ax5)
        plt.cla()
        ax5.set_ylabel('x3')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.best_popt[2], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.setp(ax5.get_xticklabels(), visible = False)
        plt.sca(ax6)
        plt.cla()
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.squid_popt[2], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.setp(ax6.get_xticklabels(), visible = False)
        plt.setp(ax6.get_yticklabels(), visible = False)
        plt.sca(ax7)
        plt.cla()
        if self.dependent == 'temperature':
            ax7.set_xlabel('Temperature (K)')
        else:
            ax7.set_xlabel('Applied field (Oe)')
        ax7.set_ylabel('x4')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.best_popt[3], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.sca(ax8)
        plt.cla()
        if self.dependent == 'temperature':
            ax8.set_xlabel('Temperature (K)')
        else:
            ax8.set_xlabel('Applied field (Oe)')
        for s in selected_scans:
            pt, = plt.plot(s.dependent, s.squid_popt[3], marker = 'o', color = 'b', picker = 10)
            pt.scan = s
        plt.setp(ax8.get_yticklabels(), visible = False)
        set_sci_format(self.pf.number)
        self.pf.canvas.draw()

    def bad_fits(self, chi2bound):
        return [s for s in self.scans if s.best_chi2 > chi2bound]

    def refit(self, scans = None, reverse = False, force_sign = 0):
        if scans is None:
            scans = self.scans
        if reverse:
            scans = scans[::-1]
        scans[0].fit(center = True)
        for i, scan in enumerate(scans):
            if not i:
                continue
            scan.fit(force_update = True, p0 = scans[i - 1].best_popt,
                    force_center = scans[i - 1].best_popt[3],
                    force_sign = force_sign)
            if i > 1:
                scan.fit(p0 = scans[i - 2].best_popt,
                        force_center = scans[i-2].best_popt[3])

    def refit_centers(self, force_sign = 0):
        start_time = time.time()
        refitted = 0
        for scan in self.scans:
            if scan.fit(center = True, force_sign = force_sign, print_new_fit = False):
                refitted += 1
        #self.plot_dependence()
        #self.plot_params()
        print("Refit {}/{} scans in {:.2f} s".format(refitted, len(self.scans), time.time()- start_time ))

class Squid_measurement:
    # class to contain individual measurements:
    def __init__(self, data, p0 = None, parent = None, dependent = 'temperature'):
        self.data = data
        self.time = data[colnames.time][0]
        self.start_temp = data[colnames.start_temp][0]
        self.end_temp = data[colnames.end_temp][0]
        self.delta_temp = abs(self.start_temp - self.end_temp)
        self.temperature = 0.5*(self.start_temp + self.end_temp)
        self.field = data[colnames.field][0]
        self.best_popt = None
        self.best_chi2 = None
        self.squid_popt = None
        # squid_chi2 for actual data, squid_fit_chi2 for fit to squid fit
        self.squid_chi2 = None
        self.squid_fit_chi2 = None
        # test for non zero:
        nzi = np.nonzero(data[colnames.scaled_vol])
        self.scaling_factor = np.average(data[colnames.scaled_vol][nzi]/data[colnames.raw_vol][nzi])
        self.position = data[colnames.pos]
        self.voltage = data[colnames.demeaned_vol] * self.scaling_factor
        self.fit_voltage = data[colnames.demeaned_vol_fit] * self.scaling_factor
        if p0 is None:
            self.fit(center = True, print_new_fit = False)
        else:
            self.fit(p0 = p0, print_new_fit = False)
        self.parent = parent
        self.dependent = {'temperature':self.temperature,
                'field':self.field}[dependent]

    def original_fit(self, center = False, force_center = None, p0 = None,
            force_sign = 0, force_update = False):
        # check if current fit is outside parameter limits
        if self.squid_popt is not None:
            if not within_lims(*self.squid_popt):
                force_update = True
        if p0 is None:
            # default fitting values
            # magnitude ~ voltage range / 2.5
            v_range = max(self.fit_voltage) - min(self.fit_voltage)
            p0 = [0, 0, v_range/2.65, -2]
        if not p0[2]:
            v_range = max(self.fit_voltage) - min(self.fit_voltage)
            p0[2] = v_range/2.65
        if center:
            centers = np.linspace(0.1, 4.1, 9)
            for c in centers:
                p0[3] = -c
                self.original_fit(p0 = p0)
        ffit = lambda pos, x1, x2, x3, x4: rso_response(pos, x1, x2, x3, x4,
                force_center = force_center, force_sign = force_sign)
        if force_center:
            p0[3] = force_center
        popt, pcov = curve_fit(ffit, self.position, self.fit_voltage, p0 = p0)
        sq_res_sum = sum((rso_response(self.position, *popt, ignore_lims = True) - self.fit_voltage)**2)
        chi2 = sq_res_sum / (len(self.position) - len(popt))
        if force_center or not self.squid_fit_chi2 or chi2 < self.squid_fit_chi2 or force_update:
            self.squid_popt, self.squid_pcov, self.squid_fit_chi2 = popt, pcov, chi2
            # calculate actual chi2
            sq_res_sum = sum((rso_response(self.position, *popt, ignore_lims = True) - self.voltage)**2)
            self.squid_chi2 = sq_res_sum / (len(self.position) - len(popt))

    def fit(self, center = False, force_center = False, p0 = None,
            print_new_fit = True, force_update = False, force_sign = 0):
        # check if current fit is outside parameter limits
        if self.best_popt is not None:
            if not within_lims(*self.best_popt):
                force_update = True
        if p0 is None:
            # default fitting values
            # magnitude ~ voltage range / 2.5
            v_range = max(self.voltage) - min(self.voltage)
            p0 = [0, 0, v_range/2.65, -2]
        if not p0[2]:
            v_range = max(self.voltage) - min(self.voltage)
            p0[2] = v_range/2.65
        old_chi2 = self.best_chi2
        old_popt = self.best_popt
        if print_new_fit:
            new_fit_msg = ("Fit parameters changed for temperature {0}:\n"
                    "x1: {1:3e} -> {5:3e}\n"
                    "x2: {2:3e} -> {6:3e}\n"
                    "x3: {3:3e} -> {7:3e}\n"
                    "x4: {4:3e} -> {8:3e}\n"
                    "chi2: {9:3e} -> {10:3e}")
        if center:
            # try to fit from muliple center values
            centers = np.linspace(0.1, 4.1, 9)
            for c in centers:
                p0[3] = -c
                self.fit(p0 = p0, print_new_fit = False,
                        force_update = False, force_sign = force_sign)
            if print_new_fit and old_chi2 != self.best_chi2:
                param_change = ([self.temperature] + list(old_popt) +
                        list(self.best_popt) + [old_chi2, self.best_chi2])
                print(new_fit_msg.format(*param_change))
            if old_chi2 != self.best_chi2:
                return True
        # fit SQUID voltage response
        ffit = lambda pos, x1, x2, x3, x4: rso_response(pos, x1, x2, x3, x4,
                force_center = force_center, force_sign = force_sign)
        if force_center:
            p0[3] = force_center
        popt, pcov = curve_fit(ffit, self.position, self.voltage, p0 = p0)
        # calculate fit
        sq_res_sum = sum((rso_response(self.position, *popt, ignore_lims = True) - self.voltage)**2)
        chi2 = sq_res_sum / (len(self.position) - len(popt))
        if print_new_fit and not center and old_chi2 != self.best_chi2:
            param_change = [self.temperature] + list(old_popt) + list(self.best_popt) + [old_chi2, self.best_chi2]
            print(new_fit_msg.format(*param_change))
        if force_update or not self.best_chi2 or chi2 < self.best_chi2:
            self.best_popt, self.best_pcov, self.best_chi2 = popt, pcov, chi2
            # not sure why this was being returned
            # return (popt, pcov, chi2)
            return True
        return False

    def plot_fit(self):
        try:
            ax1, ax2, ax3, ax4 = self.parent.rmf.get_axes()
        except:
            self.parent.figures(4)
            ax1, ax2, ax3, ax4 = self.parent.rmf.get_axes()
        # ax1 for raw voltage and best fit
        plt.sca(ax1)
        plt.cla()
        ax1.set_title('Best fit')
        ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
        ax1.set_ylabel('Scaled voltage (V)')
        plt.plot(self.position, self.voltage, linestyle = '', marker = 'o', label = 'Raw')
        y_fit = rso_response(self.position, *self.best_popt, ignore_lims = True)
        plt.plot(self.position, y_fit, label = 'Best Fit')
        plt.setp(ax1.get_xticklabels(), visible = False)
        # ax2 for raw voltage + SQUID fit
        plt.sca(ax2)
        plt.cla()
        ax2.set_title('SQUID fit')
        plt.plot(self.position, self.voltage, linestyle = '', marker = 'o', label = 'Raw')
        plt.plot(self.position, self.fit_voltage, linestyle = '', marker = 'v', label = 'SQUID fit')
        plt.setp(ax2.get_xticklabels(), visible = False)
        plt.setp(ax2.get_yticklabels(), visible = False)
        if self.squid_popt is not None:
            y_sqdfit = rso_response(self.position, *self.squid_popt, ignore_lims = True)
            plt.plot(self.position, y_sqdfit)
        # ax3 for residuals of raw + best fit
        plt.sca(ax3)
        plt.cla()
        ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
        ax3.set_xlabel('Position (cm)')
        ax3.set_ylabel('Scaled voltage residual (V)')
        y_resid = y_fit - self.voltage
        plt.plot(self.position, y_resid, linestyle = '--', marker = 'o')
        # ax4 for residuals of SQUID fit
        plt.sca(ax4)
        plt.cla()
        ax4.set_xlabel('Position (cm)')
        y_sqdresid = self.fit_voltage - self.voltage
        plt.plot(self.position, y_sqdresid, linestyle = '--', marker = 'o')
        plt.setp(ax4.get_yticklabels(), visible = False)
        set_sci_format(self.parent.rmf.number)
        self.parent.rmf.canvas.draw()

    def update_offset(self, index = 0, slope = None, plot_only = False):
        try:
            ax1, ax2 = self.parent.of.get_axes()
        except:
            self.parent.figures(3)
            ax1, ax2 = self.parent.of.get_axes()
        slopes = []
        d_vol = []
        d_pos = []
        for i in range(len(self.position) - 1):
            # calcluate slope between successive points
            d_pos.append(self.position[i + 1] - self.position[i])
            d_vol.append(self.voltage[i + 1] - self.voltage[i])
            slopes.append(d_vol[i]/d_pos[i])
        if not index:
            fit_figure = self.plot_fit()
            ax1, ax2 = self.parent.of.get_axes()
            # ax1 for data + fit
            plt.sca(ax1)
            plt.cla()
            ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
            ax1.set_ylabel('Scaled voltage (V)')
            plt.plot(self.position, self.voltage, linestyle = '--', marker = 'o')
            plt.setp(ax1.get_xticklabels(), visible = False)
            # ax2 for derivative plot
            plt.sca(ax2)
            plt.cla()
            ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
            ax2.set_xlabel('Position (cm)')
            ax2.set_ylabel('dV/dx (V/cm)')
            mid_pos = np.array(self.position[:-1]) + np.array(d_pos)
            plt.plot(mid_pos, slopes, linestyle = '--', marker = 'o')
            for i in range(len(mid_pos)):
                plt.text(mid_pos[i], slopes[i], str(i))
            if not plot_only:
                input_index = input("Index to correct: ")
                input_slope = input("Correct slope (blank for auto): ")
                slope = float(input_slope) if input_slope else None
                if input_index:
                    self.update_offset(index = int(input_index), slope = slope)
            set_sci_format(self.parent.of.number)
            self.parent.of.canvas.draw
        else:
            if slope is None:
                slope = 0.5 * (slopes[index - i] + slopes[index + 1])
            slope_fix = slope * d_pos[index]
            for j in range(index + 1, len(self.voltage)):
                self.voltage[j] = self.voltage[j] - d_vol[index] + slope_fix
            self.fit(center = True)
            ax1, ax2 = self.parent.of.get_axes()
            plt.sca(ax1)
            plt.plot(self.position, self.voltage, linestyle = '--', marker = 'o')
            plt.sca(ax2)
            mid_pos = np.array(self.position[:-1]) + np.array(d_pos)
            plt.plot(mid_pos, slopes, linestyle = '--', marker = 'o')
            self.parent.of.canvas.draw()

    def get_best_moment(self):
        return self.best_popt[2] * squid_factor

    def get_squid_moment(self):
        return self.squid_popt[2] * squid_factor

    def reset_offset(self):
        self.voltage = self.data[vol]

def rso_response(pos, x1, x2, x3, x4, force_center = None, force_sign = 0, ignore_lims = False):
    # rso scans start and end in middle of scan
    # account for voltage drift as function of time
    if not ignore_lims:
        if x1_min is not None:
            if x1 < x1_min:
                return 10000
        if x2_min is not None:
            if x2 < x2_min:
                return 10000
        if x3_min is not None:
            if x3 < x3_min:
                return 10000
        if x4_min is not None:
            if x4 < x4_min:
                return 10000
        if x1_max is not None:
            if x1 > x1_max:
                return 10000
        if x2_max is not None:
            if x2 > x2_max:
                return 10000
        if x3_max is not None:
            if x3 > x3_max:
                return 10000
        if x4_max is not None:
            if x4 > x4_max:
                return 10000
    if force_center:
        # force x4 within 0.02 cm of given value
        if force_center - 0.02 > x4 or force_center + 0.02 < x4:
            return 10000
    if force_sign:
        # force sign of magnetization (x3)
        if force_sign * x3 < 0:
            return 10000
    # constants
    R = 0.97
    L = 1.519
    # index array for voltage drift
    index_array = np.linspace(0, 1, len(pos))
    X = R**2 + (pos + x4)**2
    Y = R**2 + (L + (pos + x4))**2
    Z = R**2 + (-L + (pos + x4))**2
    return x1 + x2 * index_array + x3 * (2 * X**(-3/2) - Y**(-3/2) - Z**(-3/2))

def load_file(fname = None, dep = 'temperature'):
    if fname is None:
        fname = filedialog.askopenfilename()
    if fname is None:
        return
    sqd_data = Squid(fname, dep)
    # sqd_data.plot_dependence()
    return sqd_data

def dump_fit(squid, fname = None):
    if fname is None:
        fname = filedialog.asksaveasfile(mode = 'w', defaultextension = '.csv',
                filetypes = [('Comma separated values', '.csv'), ('All files', '.*')])
    if fname is None:
        return
    columnnames = ['Scan index', 'Time (s)', 'Field (Oe)', 'Temperature (K)',
            'Long Moment fitted (EMU)', 'Long Moment Original (EMU)',
            'Delta T (K)', 'Reduced chi2 fitted', 'Reduced chi2 original']
    fname.write(','.join(columnnames) + '\n')
    for i, s in enumerate(squid.scans):
        data = [i, s.time, s.field, s.temperature, s.get_best_moment(),
                s.get_squid_moment(), s.delta_temp, s.best_chi2, s.squid_chi2]
        fname.write(','.join([str(d) for d in data]) + '\n')
    fname.close()

def main(args):
    window = Window()
    time.sleep(1)
    dep = {0:'temperature', 1:'field'}[args.dependence]
    if args.FILE is None or args.FILE == []:
        sqd_data = [load_file(dep = dep)]
    else:
        sqd_data = [load_file(f, dep) for f in args.FILE]
    s = sqd_data[0]
    window.load_squid(s)
    plt.ion()
    time.sleep(1)
    plt.show()
    time.sleep(1)
    embed(display_banner = False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("FILE", metavar = 'input data files', nargs = '*', default = None,
            help = "Raw squid data files")
    parser.add_argument("--dependence", action = "store_true",
            help = "Changes magnetic moment dependence to Field if specified (default to temperature)")
    args = parser.parse_args()
    main(args)
