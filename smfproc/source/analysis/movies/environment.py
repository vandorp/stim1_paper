"""
Properties of experimental hardware in different labs. 
"""

# Python packages
import numpy as np

# Program modules
from smfproc.source.io.datatypes import TimeSeries


class HardWare(object):
    """ Base class for hardware properties."""

    def __init__(self, mFile):
        self.movieFile = mFile
        self.get = mFile.get_timeseries
        self.devices = {
                'lasers':{},
                'shutters':{}}
    
    def add_device(self, devtype, devname, output, state):
        device = {
                'output': [output, float],
                'state': [state, bool]}
        self.devices[devtype][devname] = device

    def has_device(self, devtype, devname):
        if self.devices.has_key(devtype):
            return self.devices[devtype].has_key(devname)
        else:
            return False

    def device_is_active(self, devtype, devname):
        if self.has_device(devtype, devname):
            dev = self.get(*self.devices[devtype][devname]["output"])
            return dev.size > 0
        else:
            return False

    def get_laser_intensity(self):
        return np.array([])


class BrungerLab(HardWare):
    """ Brunger lab hardware. """

    def __init__(self, mFile):
        HardWare.__init__(self, mFile)
    


class LewisLab(HardWare):
    """ Lewis lab hardware. """

    shutterCodes = {
            'laser488':16,
            'laser532':8,
            'laser638':32}

    def __init__(self, mFile):
        HardWare.__init__(self, mFile)
        
        self.add_device(
                "lasers", 
                "laser488", 
                "Sapphire-PowerSetpoint",
                "Sapphire-State")

        self.add_device(
                "lasers", 
                "laser532", 
                "Laser Obis 532-PowerReadback",
                "Laser Obis 532-State")

        self.add_device(
                "lasers", 
                "laser638", 
                "Laser Obis 638-PowerReadback",
                "Laser Obis 638-State")

        self.add_device(
                "shutters",
                "arduino",
                "Arduino-Switch-State",
                "Arduino-Shutter-OnOff")


    def get_laser_intensity(self, name, shutter):
        dev = self.devices

        laser = self.get(*dev["lasers"][name]["output"])
        shutter = self.get(*dev["shutters"]["arduino"]["output"]) == shutter
        if laser.size>0:
            return shutter*laser
        else:
            return np.array([])
    
    def get_all_lasers(self):
        intens = []
        time = []
        descr = []
        for laser in self.shutterCodes:
            if self.device_is_active("lasers", laser):
                laserIntens = self.get_laser_intensity(
                    laser, self.shutterCodes[laser])
                intens.append(laserIntens)
                descr.append(laser)
                
                time = np.arange(len(laserIntens))*self.movieFile.header['dt']
        return TimeSeries(intens, time), descr

    def get_all_devices(self):
        lasers, descr = self.get_all_lasers()

        dt = self.movieFile.header['dt']
        nrFrames = self.movieFile.header['nrFrames']
        time = np.arange(nrFrames)*dt

        descr = ['time']+descr
        #devices = TimeSeries(np.vstack((lasers.time, lasers.data)))
        
        if len(descr)>1:
            devices = TimeSeries(np.vstack((time, lasers.data)))
        else:
            descr = descr+['laser532', 'laser638']
            devices = TimeSeries(np.vstack((time, np.ones(time.shape), np.zeros(time.shape))))

        return devices, descr

