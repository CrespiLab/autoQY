# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:41:10 2025

@author: jorst136
Modified from LiYuanhe (liyuanhe211/Python_Lib/My_Lib_PyQt.py)
"""
# import sys
# import pathlib
import platform
import os
from PyQt5 import Qt
# from PyQt5 import uic
from PyQt5.Qt import QApplication
# import platform
import traceback

QApplication.setAttribute(Qt.Qt.AA_EnableHighDpiScaling, True)
QApplication.setAttribute(Qt.Qt.AA_UseHighDpiPixmaps)

def set_Windows_scaling_factor_env_var():
    # Sometimes, the scaling factor of PyQt is different from the Windows system scaling factor, reason unknown
    # For example, on a 4K screen sets to 250% scaling on Windows, PyQt reads a default 300% scaling,
    # causing everything to be too large, this function is to determine the ratio of the real DPI and the PyQt DPI

    if platform.system() == 'Windows':
        import ctypes
        try:
            import win32api
            MDT_EFFECTIVE_DPI = 0
            monitor = win32api.EnumDisplayMonitors()[0]
            dpiX, dpiY = ctypes.c_uint(), ctypes.c_uint()
            ctypes.windll.shcore.GetDpiForMonitor(monitor[0].handle, MDT_EFFECTIVE_DPI, ctypes.byref(dpiX), ctypes.byref(dpiY))
            DPI_ratio_for_monitor = (dpiX.value + dpiY.value) / 2 / 96
        except Exception as e:
            traceback.print_exc()
            print(e)
            DPI_ratio_for_monitor = 0

        DPI_ratio_for_device = ctypes.windll.shcore.GetScaleFactorForDevice(0) / 100
        PyQt_scaling_ratio = QApplication.primaryScreen().devicePixelRatio()
        print(f"Windows 10 High-DPI debug:", end=' ')
        Windows_DPI_ratio = DPI_ratio_for_monitor if DPI_ratio_for_monitor else DPI_ratio_for_device
        if DPI_ratio_for_monitor:
            print("Using monitor DPI.")
            ratio_of_ratio = DPI_ratio_for_monitor / PyQt_scaling_ratio
        else:
            print("Using device DPI.")
            ratio_of_ratio = DPI_ratio_for_device / PyQt_scaling_ratio

        if ratio_of_ratio > 1.05 or ratio_of_ratio < 0.95:
            use_ratio = "{:.2f}".format(ratio_of_ratio)
            print(f"{DPI_ratio_for_monitor=}, {DPI_ratio_for_device=}, {PyQt_scaling_ratio=}")
            print(f"Using GUI high-DPI ratio: {use_ratio}")
            print("----------------------------------------------------------------------------")
            os.environ["QT_SCALE_FACTOR"] = use_ratio
        else:
            print("Ratio of ratio near 1. Not scaling.")

        return Windows_DPI_ratio, PyQt_scaling_ratio


def get_matplotlib_DPI_setting(Windows_DPI_ratio):
    matplotlib_DPI_setting = 60
    if platform.system() == 'Windows':
        matplotlib_DPI_setting = 60 / Windows_DPI_ratio
    if os.path.isfile("__matplotlib_DPI_Manual_Setting.txt"):
        matplotlib_DPI_manual_setting = open("__matplotlib_DPI_Manual_Setting.txt").read()
        # if is_int(matplotlib_DPI_manual_setting):
        if int(matplotlib_DPI_manual_setting):
            matplotlib_DPI_setting = matplotlib_DPI_manual_setting
    else:
        with open("__matplotlib_DPI_Manual_Setting.txt", 'w') as matplotlib_DPI_Manual_Setting_file:
            matplotlib_DPI_Manual_Setting_file.write("")
    matplotlib_DPI_setting = int(matplotlib_DPI_setting)
    print(
        f"\nMatplotlib DPI: {matplotlib_DPI_setting}. \nSet an appropriate integer in __matplotlib_DPI_Manual_Setting.txt if the preview size doesn't match the output.\n")

    return matplotlib_DPI_setting
