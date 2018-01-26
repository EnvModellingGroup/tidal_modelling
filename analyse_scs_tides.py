#!/usr/bin/env python
import vtktools
import vtk
import datetime
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import argparse
from scipy.stats import linregress
import uptide
import uptide.netcdf_reader as netcdf
import csv
import sys
import os
import GFD_basisChange_tools as sp
import zipfile

#
# Caldwell, P. C., M. A. Merrifield, P. R. Thompson (2015), 
# Sea level measured by tide gauges from global oceans - the 
# Joint Archive for Sea Level holdings (NCEI Accession 0019568), 
# Version 5.5, NOAA National Centers for Environmental Information, 
# Dataset, doi:10.7289/V5V40S7W.
#
# Assumes the PVTU is on the sphere. Will not work with UTM-type projections, hence SW tides

startdate = datetime.datetime(1700, 1, 1, 0, 0)

def main():


    parser = argparse.ArgumentParser(
         prog="analyse tidal model run",
         description="""Analyse a Fluidity model "on the sphere" against tidal gauge data"""
    )
    parser.add_argument(
        '-v', 
        '--verbose', 
        action='store_true', 
        help="Verbose output: mainly progress reports.",
        default=False
    )
    parser.add_argument(
        '-n', 
        '--netcdf', 
        action='store_true', 
        help="The input file for the model is a NetCDF. Default is a (P)VTU. If NetCDF, fields should be M2Amp, M2Phase, etc, as in PVTU",
        default=False
    )
    parser.add_argument(
        '-f', 
        '--fes', 
        help="Also compare against FES data. Note the FES data should be in 'AMCG' NetCDF format. Give filename.",
        default=None
    )
    parser.add_argument(
        'tide_gauge',
        metavar='tide_gauge',
        help='Tide gauge data. Can be a CSV with sensible headers (see script comments), a directory (assumes from http://www.psmsl.org/data/hf/) or a list of files.',
        nargs='*'
    )
    parser.add_argument(
        'model_input',
        metavar='model_input',
        help='Either a (P)VTU file that contains the tidal data or a NetCDF file that contains same.'
    )
    parser.add_argument(
        'output_stub',
        metavar='output_stub',
        help="The output file stub. Five files will be generated. A cross plot of amplitudes (Model against Obs, but possibly Model-FES and FES-Obs), polar plot of phases, CSV file of error measures (D_n) and two CSV files of stations used for Fluidity and FEs comparison"
    )

    
    args = parser.parse_args()
    verbose = args.verbose
    netcdf_input = args.netcdf
    fes_input = args.fes
    model_input = args.model_input
    output_stub = args.output_stub
    tide_gauges = args.tide_gauge

    # first find the known tidal constituents in the model data
    if model_input.endswith(".pvtu") or model_input.endswith(".vtu"):
        # This assumes the are M2Amp or M2Phase etc
        vt_file=vtktools.vtu(model_input)
        pvtu_vars = vt_file.GetFieldNames()
        tidal_components = uptide.tidal.omega.keys()
        known_tidal_components = []
        tidal_vars = []
        for t in tidal_components:
            if t+"amp" and t+"phase" in pvtu_vars:
                known_tidal_components.append(t)
                tidal_vars.append(t+"amp")
                tidal_vars.append(t+"phase")

    known_tidal_components.sort()
    tidal_vars.sort()

    # read in tide gauge data, so let's do easiest to hardest...
    # how long is the input? 
    tide_gauge_data = {}
    if len(tide_gauges) == 1:
        # if dir - then this is an unzipped download from www.psmsl.org, deal with accordingly
        if os.path.isdir(tide_gauges[0]):
            tide_gauge_data = read_from_pmsl(tide_gauges[0], known_tidal_components)
        else:
            # if file check if a csv file
            try:
                with open(tide_gauges[0], 'r') as csvfile:
                    dialect = csv.Sniffer().sniff(csvfile.read(1024))
                    csvfile.seek(0)
                    # Expect file to be:
                    # Name, lat, lon, M2Amp, M2Phase, etc
                    # Header should be as above, with capitalisation etc, but order is unimportant
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        temp = dict(row) # copy
                        temp.pop['Name'] # remove name
                        tide_gauge_data[row['Name']] = temp
            except csv.Error:
                # it's not, so no idea what the heck has been thrown at us
                print "Sorry, I could not decipher your tide gauge data. Exiting."
                sys.exit(1)
    else:
        # list of fles - read each in turn and dump in hash. 
        # NOT YET IMPLEMENTED
        print "Sorry, a list of files is not yet implemented. Code Python? Can you add this please!"
        sys.exit()

    # now grab data from FES if requested
    if fes_input is not None:
        fes_data = []
        lats = []
        lons = []
        for t in sorted(tide_gauge_data): # loop in alphabetical order of names
            lats.append(tide_gauge_data[t]['lat'])
            lons.append(tide_gauge_data[t]['lon'])

        nci = netcdf.NetCDFInterpolator(fes_input, ('longitude', 'latitude'), ('longitude', 'latitude'))
        fes_data = []
        for v in tidal_vars:
            nci.set_field(v)
            data = []
            for ll,lt in zip(lons,lats):
                data.append(nci.get_val((ll, lt),allow_extrapolation=False)/100.0) # in cm, so convert to m
            fes_data.append([v, data])


    # tide gauge data comes back as two nested dicts
    #{location_name: {M2Amp: x, M2Phase:, y, etc, etc} 

    # deal with fluidity - tells us what constits to grab later
    if model_input.endswith(".pvtu") or model_input.endswith(".vtu"):
        # if pvtu - extract points at the locations.
        # grab arrays of all lats and lons
        lats = []
        lons = []
        for t in sorted(tide_gauge_data): # loop in alphabetical order of names
            lats.append(tide_gauge_data[t]['lat'])
            lons.append(tide_gauge_data[t]['lon'])

        # grab the values of each variable at each location (not sure what the function passess back yet - need to test)
        model_data = probe_spherical_vtu(model_input, tidal_vars, 6.37101e+06, 0, lats, lons, False)
    else:
        # try opening as a NetCDF
        nci = uptide.NetCDFInterpolator(model_input, ('lon', 'lat'), ('longitude', 'latitude'))
        
        # loop through the lat/lon of the tide gauges and extract values.
        pass
        # implement this!


    # perform error analysis
    errors = []
    av_err = []
    i = 0
    average_amp = []
    for t in known_tidal_components:
        obs_amps = []
        obs_phases = []
	lons = []
        lats = []
        names = []
        for l in sorted(tide_gauge_data):
            obs_amps.append(tide_gauge_data[l][t+"amp"])
            obs_phases.append(tide_gauge_data[l][t+"phase"])
	    lats.append(tide_gauge_data[l]['lon'])
            lons.append(tide_gauge_data[l]['lat'])
            names.append(l)

        model_amps = fes_data[i][1]
        model_phases = fes_data[i+1][1]
        obs_amps = np.array(obs_amps)
        obs_phases = np.array(obs_phases)
 	lats = np.array(lats)
        lons = np.array(lons)
        names = np.array(names)
        
        index_to_remove = []
        for j in range(0,len(model_amps)):
            #if np.isnan(model_amps[j]) or model_amps[j] > 10000 or model_amps[j] == '--': #null value
            if type(model_amps[j]) is np.ma.core.MaskedConstant:
                index_to_remove.append(j)

        model_amps = np.delete(model_amps, index_to_remove)
        model_phases = np.delete(model_phases, index_to_remove)
        obs_amps = np.delete(obs_amps, index_to_remove)
        obs_phases = np.delete(obs_phases, index_to_remove)
	lats = np.delete(lats, index_to_remove)
        lons = np.delete(lons, index_to_remove)
        names = np.delete(names, index_to_remove)
        average_amp.append(np.average(obs_amps))

        errors.append(uptide.analysis.error_analysis(model_amps, model_phases, obs_amps, obs_phases)[0])
        av_err.append(uptide.analysis.error_analysis(model_amps, model_phases, obs_amps, obs_phases)[1])        
        i+=2 # 'cos we have amp and phase
      
    average_amp = np.array(average_amp)
    print "Error to FES"
    print "components:", known_tidal_components
    print "error:", errors
    print "Gauge av amp:", average_amp
    print "Relative err:", 1. - (errors / average_amp)
    print "No. stations valid:", len(obs_amps)
    
    with open(output_stub+"_fes_stations.csv", 'w') as f:
    	writer = csv.writer(f)
        for lt,ln,name in zip(lats,lons,names):
    	    writer.writerow([name,ln,lt])

    errors = []
    av_err = []
    i = 0
    average_amp = []
    for t in known_tidal_components:
        obs_amps = []
        obs_phases = []
	lons = []
        lats = []
        names = []
        for l in sorted(tide_gauge_data):
            obs_amps.append(tide_gauge_data[l][t+"amp"])
            obs_phases.append(tide_gauge_data[l][t+"phase"])
	    lats.append(tide_gauge_data[l]['lon'])
            lons.append(tide_gauge_data[l]['lat'])
            names.append(l)

        model_amps = model_data[i][1]
        model_phases = model_data[i+1][1]
        obs_amps = np.array(obs_amps)
        obs_phases = np.array(obs_phases)
 	lats = np.array(lats)
        lons = np.array(lons)
        names = np.array(names)
        
        index_to_remove = []
        for j in range(0,len(model_amps)):
            if np.isnan(model_amps[j]) or model_amps[j] > 10000: #null value
                index_to_remove.append(j)

        model_amps = np.delete(model_amps, index_to_remove)
        model_phases = np.delete(model_phases, index_to_remove)
        obs_amps = np.delete(obs_amps, index_to_remove)
        obs_phases = np.delete(obs_phases, index_to_remove)
	lats = np.delete(lats, index_to_remove)
        lons = np.delete(lons, index_to_remove)
        names = np.delete(names, index_to_remove)
        average_amp.append(np.average(obs_amps))

        errors.append(uptide.analysis.error_analysis(model_amps, model_phases, obs_amps, obs_phases)[0])
        av_err.append(uptide.analysis.error_analysis(model_amps, model_phases, obs_amps, obs_phases)[1])

        i+=2 # 'cos we have amp and phase
      
    average_amp = np.array(average_amp)
    print "Error to Model"
    print "components:", known_tidal_components
    print "error:", errors
    print "Gauge av amp:", average_amp
    print "Relative err:", 1. - (errors / average_amp)
    print "No. stations valid:", len(obs_amps)
    with open(output_stub+"_fluidity_stations.csv", 'w') as f:
    	writer = csv.writer(f)
        for lt,ln,name in zip(lats,lons,names):
    	    writer.writerow([name,ln,lt])


    # fluidity to tide gauge 
    fig = plt.figure(figsize=(15,15),dpi=180)
    i = 0
    p = 1
    for t in known_tidal_components:
        obs_amps = []
        obs_phases = []
        for l in sorted(tide_gauge_data):
            obs_amps.append(tide_gauge_data[l][t+"amp"])
            obs_phases.append(tide_gauge_data[l][t+"phase"])

        model_amps = model_data[i][1]
        obs_amps = np.array(obs_amps)
        
        index_to_remove = []
        for j in range(0,len(model_amps)):
            if np.isnan(model_amps[j]) or model_amps[j] > 10000: #null value
                index_to_remove.append(j)

        model_amps = np.delete(model_amps, index_to_remove)
        obs_amps = np.delete(obs_amps, index_to_remove)

        ax = fig.add_subplot(2,2,p)

        gradient, intercept, r_value, p_value, std_err = linregress(model_amps,obs_amps)
        ax.plot(model_amps,obs_amps,'bx')
        yLim = ax.get_ylim()
        xLim = ax.get_xlim()
        lim = max((xLim, yLim))
        ax.plot(lim, lim, 'k-')
        ax.set_title(known_tidal_components[p-1]+" correlation")
        plt.axis('equal')
        ax.set_xlabel("Model (m)")
        ax.set_ylabel("Tide Gauges (m)")
        i += 2 #  counter for model data

    plt.savefig(output_stub+"_fluidity_obs.pdf", dpi=180)
     
    # fluidity to fes
    fig = plt.figure(figsize=(15,15),dpi=180)
    i = 0
    p = 1
    for t in known_tidal_components:

        model_amps = model_data[i][1]
        obs_amps = fes_data[i][1]
        
        index_to_remove = []
        for j in range(0,len(model_amps)):
            if ((np.isnan(model_amps[j]) or model_amps[j] > 10000) or 
                (np.isnan(obs_amps[j]) or obs_amps[j] > 10000)): #null value
                    index_to_remove.append(j)

        model_amps = np.delete(model_amps, index_to_remove)
        obs_amps = np.delete(obs_amps, index_to_remove)

        ax = fig.add_subplot(2,2,p)

        gradient, intercept, r_value, p_value, std_err = linregress(model_amps,obs_amps)
        ax.plot(model_amps,obs_amps,'bx')
        yLim = ax.get_ylim()
        xLim = ax.get_xlim()
        lim = max((xLim, yLim))
        ax.plot(lim, lim, 'k-')
        ax.set_title(known_tidal_components[p-1]+" correlation")
        plt.axis('equal')
        ax.set_xlabel("Model (m)")
        ax.set_ylabel("FES 2012 (m)")
        i += 2 #  counter for model data

    plt.savefig(output_stub+"_fluidity_fes.pdf", dpi=180)
   

def rmse(predictions, targets):
    predictions = np.array(predictions)
    targets = np.array(targets)
    return np.sqrt(((predictions - targets) ** 2).mean())


def read_from_pmsl(directory, constituents):

    import glob

    # the directory structure is this:
    # root--+
    #       |
    #       +-atlantic
    #       |
    #       +-indian
    #       |
    #       +-pacific-+
    #                 |
    #                 +-hourly-+
    #                          |
    #                          + h001a.zip
    #                          |
    #                          + h699a.zip
    #
    # The zip files have the following naming convention
    # hSSSv.zip where SSS is the station number and v is the series ID
    # Unzipping them gives:
    #   - h641a75.dat
    #   - h641a76.dat
    #   - etc
    # where 641 is the station number, a is the series ID and 75 is the year
    # The first letter (h) is 1900-1999, i is 2000-2999 and g is 1800-1899
    # So we want to traverse the directory structure, getting all zip files
    # within the zip file, find the last year (they only put up complete years)
    #
    # once found, the file has the following structure:
    # 302A Balboa             Panama              2016 08579N 079344W 0000 3 00000R MM
    # station number        1-3     3     exactly 3 digits
    # station version         4     1     letter from A to Z 
    # station name         6-23    18     Abbreviated if necessary  
    # region/country      25-43    19     Abbreviated if necessary     
    # year                45-48     4
    # latitude            50-55     6     degrees, minutes, tenths
    #                                      (implied decimal), and hemisphere
    # longitude           57-63     7     degrees, minutes, tenths
    #                                      (implied decimal), and hemisphere
    # GMT offset          65-68     4     time data are related to in terms
    #                                      of difference from GMT in hours
    #  	                                   and tenths (implied decimal) with
    #					   East longitudes positive*
    # decimation method      70     1     Coded as 1 : filtered
    #                                              2 : simple average of all samples
    #						   3 : spot readings
    #						   4 : other
    # reference offset    72-76     5     constant offset to be added to each data value 
    #                                      for data to be relative to tide staff zero or primary
    #					   datum in same units as data
    # reference code         77     1     R = data referenced to datum
    #                                     X = data not referenced to datum
    # units               79-80     2     always millimeters, MM
    #
    # The file then has the data arrange like:
    # 302A Balb  201601011 5641 5956 5867 5346 4637 3737 2905 2288 2272 2715 3555 4450
    # 302A Balb  201601012 5244 5628 5595 5188 4522 3768 3016 2469 2341 2560 3269 4184
    #
    # So each day has two lines, with the date in the third block (minus the digit), then
    # the time is 1 2 3 4, etc to give a full 24 hours
    tide_gauge_data = {}

    oceans = ['atlantic', 'indian', 'pacific']
    for o in oceans:
        zipfiles = glob.glob(os.path.join(directory,o,"hourly","*.zip"))
        for z in zipfiles:
            # unzip and list files
            with zipfile.ZipFile(z, 'r') as myzip:
                files = myzip.namelist()
                # find the file we want...
                files.sort()
                # grab last one!
                filename = files[-1]
                # files are created with bytes as the chuncks, so we use the
                # fact that strings are lists and partition by character.
                with myzip.open(filename) as f:
                    i = 0
                    sl_data = []
                    sl_time = []
                    nam = ''
                    for line in f:
                        if i == 0:
                            # header
                            name = line[5:24].rstrip()
                            name.replace(" ","_")
                            lat = line[49:56].rstrip()
                            lon = line[56:64].rstrip()
                            lat = dms_to_dd(lat[0:2], lat[2:4], lat[4], lat[5])
                            lon = dms_to_dd(lon[0:3], lon[3:5], lon[5], lon[6])
                            tide_gauge_data[name] = {'lat': lat, 'lon': lon}
                            i += 1
                        else:
                            # data. Split by characters, then remove all whitespace
                            date = line[11:19]
                            date = datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8]))
                            if line[19] == '1':
                                time_addition = 0
                            else:
                                time_addition = 12
                            data = line[20:]
                            n = 5
                            # seperate out into chucks of the right size (- 4 digits)
                            data = [data[i:i+n] for i in range(0, len(data), n)]
                            data = data[:12] # removes trailing crap
                            hour = 0
                            for d in data:
                                time = datetime.time(hour+time_addition)
                                sl_time.append(datetime.datetime.combine(date, time))
                                sl_data.append(float(d) / 1000.) #mm -> m
                                hour += 1
                    # now analyse to get the components we want
                    amps, phases = analysis_tide(sl_time, sl_data, constituents)
                    temp = tide_gauge_data[name]
                    for c, a, p in zip(constituents, amps, phases): #in the same order as amp and phases
                        temp[c+"amp"] = a
                        temp[c+"phase"] = p
                        temp['year'] = date.year
                    tide_gauge_data[name] = temp

    return tide_gauge_data


def dms_to_dd(degrees, minutes, decimal, direction):
    dd = float(degrees) + float(minutes+'.'+decimal)/60.;
    if direction == 'W' or direction == 'S':
        dd *= -1
    return dd    


def analysis_tide(times, levels, constituents):

        # times needs to be seconds since start
        times_s = []
        for t in times:
            times_s.append((t - times[0]).total_seconds())
        # remove NaN data
        i = 0
        remove = []
        for s in levels:
            if s == None:
                remove.append(i)
            i += 1
        times_s  = np.delete(times_s, remove)
        ssh      = np.delete(levels, remove)
        # now we can push this through uptide
        tide = uptide.Tides(constituents)
        tide.set_initial_time(times[0])
        amp,phase = uptide.analysis.harmonic_analysis(tide, ssh, times_s)

        return amp, phase


def probe_spherical_vtu(input_file, variable_metadata, sphereSurfaceRadius, sliceDepth, lats, lons, verbose):
    # largely written by Alexandros Avdis (https://github.com/AlexandrosAvdis) and copied from older scripts

    import copy
    # Express output grid into polar-stereographic and prepare into
    # a format used by vtktools probing routine
    grid_polarStereographicCoords = []
    grid_lonLatRadCoords = []
    for latitude,longitude in zip(lats,lons):
        grid_lonLatRadCoords.append([longitude, latitude, sphereSurfaceRadius-sliceDepth])
        polarStereographicCoords = sp.lonlatradius_2_polarStereographic(
                    [longitude, latitude, sphereSurfaceRadius-sliceDepth])
        polarStereographicCoords[2] = np.double(0)   #Force transformed surface to
                                          # be on z=0 for probing to work.
        grid_polarStereographicCoords.append([polarStereographicCoords[0],
                                              polarStereographicCoords[1],
                                              polarStereographicCoords[2]]
                                      )
    grid_polarStereographicCoords = vtktools.arr(grid_polarStereographicCoords)
    grid_lonLatRadCoords = vtktools.arr(grid_lonLatRadCoords)

    #Open (p)vtu and check all requested variables actually exist.
    if (verbose):
        sys.stdout.write('Probing '+input_file+' ... \n')
        sys.stdout.flush()
    vt_file=vtktools.vtu(input_file)
    vtk_npoints = vt_file.ugrid.GetNumberOfPoints()
    for variable_name in variable_metadata:
        if variable_name not in vt_file.GetFieldNames():
            raise Exception('Error: variable '+variable_name+' not in file.')
    if (verbose):
        sys.stdout.write('    Input grid has '+str(vtk_npoints)+' points. \n')
        sys.stdout.write('    Input grid has '+str(vt_file.ugrid.GetNumberOfCells())+' cells. \n')
        sys.stdout.flush()
    #Remove unwanted variables to make probing faster
    if (verbose):
        sys.stdout.write('    Removing not-probed variables ... \n')
        sys.stdout.flush()
    for var in vt_file.GetFieldNames():
        if var not in variable_metadata:
            vt_file.RemoveField(var)
    # Transform the point coordinates into polar-stereographic.
    # However also keep a copy of the cartesian coordinates, they
    # are used later on.
    if (verbose):
        sys.stdout.write('    Transforming coordinates from (p)vtu ... \n')
        sys.stdout.flush()
    VT_CartesianCoords_List=[]
    for i in range (vtk_npoints):
      (x,y,z) = vt_file.ugrid.GetPoint(i)
      VT_CartesianCoords_List.append([x,y,z])
      new_x, new_y, new_z = sp.cartesian_2_polarStereographic([x,y,z]) #MUST ADD A SURFACE RADIUS ARGUMENT
      new_z = np.double(0.0) #Force transformed surface to be on z=0 for probing to work.
      vt_file.ugrid.GetPoints().SetPoint(i, new_x, new_y, new_z)
    # Now probe the data. If probing at a non-zero depth we must probe
    # the full data, which can be very slow for large data-sets. If
    # we are probing at the surface (depth=0) the process is made
    # much faster by extracting the surface of the domain and then
    # probing that surface.
    if sliceDepth!=0:
        if (verbose):
            sys.stdout.write('    Probing the data ... ')
            sys.stdout.flush()
        probed_variableData = vt_file.ProbeData(grid_polarStereographicCoords, variable_name)
    else:
        # Initialise probe
        points = vtk.vtkPoints()
        points.SetDataTypeToDouble()
        ilen, jlen = grid_polarStereographicCoords.shape
        for i in range(ilen):
            points.InsertNextPoint(grid_polarStereographicCoords[i][0], grid_polarStereographicCoords[i][1], grid_polarStereographicCoords[i][2])
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        probe = vtk.vtkProbeFilter()
        try:
            probe.SetInputData(polydata)
            probe.SetSourceData(vt_file.ugrid)
        except AttributeError:
            probe.SetInput(polydata)
            probe.SetSource(vt_file.ugrid)
        if (verbose):
            sys.stdout.write('    Probing the data ... \n')
            sys.stdout.flush()
        probe.Update()
        # Generate a list invalidNodes, containing the indices of points
        # outside the domain.
        valid_ids = probe.GetValidPoints()
        valid_loc = 0
        invalidNodes = []
        for i in range(ilen):
          if valid_ids.GetTuple1(valid_loc) == i:
            valid_loc += 1
          else:
            invalidNodes.append(i)
        if (verbose):
            sys.stdout.write('        Found '+str(len(invalidNodes))+' points outside the domain.\n')
            sys.stdout.flush()
        # Get final updated values
        pointdata=probe.GetOutput().GetPointData()
        probed_variables = []
        for variable in variable_metadata:
            vtk_variable_name = variable
            if (verbose):
                sys.stdout.write('      Processing variable '+vtk_variable_name+' ...\n')
                sys.stdout.flush()
            vtkdata=pointdata.GetArray(vtk_variable_name)
            nComponents = vtkdata.GetNumberOfComponents()
            nTuples = vtkdata.GetNumberOfTuples()
            probed_variableData = np.array([vtkdata.GetValue(i) for i in range(nTuples * nComponents)])
            #Assign Nan's to "invalid" points
            if (verbose):
                sys.stdout.write("        Assigning NaN's to points outside the domain ... \n")
                sys.stdout.flush()
            for invalidNode in invalidNodes:
                probed_variableData[invalidNode] = np.nan
                #Reshape data and store to final object prior to returning
            probed_variables.append([vtk_variable_name, probed_variableData])
    if (verbose):
        sys.stdout.write('Done probing.\n')
        sys.stdout.flush()

    return probed_variables


if __name__ == "__main__":
    main()




