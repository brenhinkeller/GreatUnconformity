# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:54:57 2016

@author: jhusson
"""

################ IMPORT MODULES AND DEFINE FUNCTIONS
import urllib2
import json
import numpy as np
import matplotlib.pyplot as plt
import pickle
from show_color import show_color,cnames
from addGTS import addGTS
from scipy import ndimage


#%% LOAD FROM MACROSTRAT API

def json2struct(myjson,mynames,myformats):
    mylist = np.zeros(len(myjson),dtype={'names':mynames,'formats':myformats})
    for i,item in enumerate(myjson):
        for name in mynames:
            if not item[name]:
                mylist[name][i]=0
            else:
                mylist[name][i]=item[name]            
    return mylist


#QUERY THE API FOR SEDIMENTARY UNITS
api =  'https://macrostrat.org/api/v2/'
u_query = 'units?lith_class=sedimentary' +\
                '&lith_type=metasedimentary'+\
                '&project_id=1' #+\
#                '&project_id=5' +\
#                '&project_id=6' +\
#                '&project_id=7' +\
#                '&environ_class=marine'

response = urllib2.urlopen(api+u_query)
sed = json.load(response) 
sed=sed['success']['data']

#CONVERT TO STRUCTURED ARRAY
mynames = ['unit_id','project_id','col_id','col_area','min_thick','max_thick','strat_name_id','unit_name','b_age','t_age','pbdb_occurrences']
myformats=['i4','i4','i4','f4','f4','f4','i4','|S100','f4','f4','f4']

sed = json2struct(sed,mynames,myformats)

cols=np.unique(sed['col_id'])


#QUERY THE API FOR MARINE SEDIMENTARY UNITS
api =  'https://macrostrat.org/api/v2/'
u_query = 'units?lith_class=sedimentary' +\
                '&lith_type=metasedimentary'+\
                '&project_id=1'+\
                '&environ_class=marine'

#                '&project_id=5' +\
#                '&project_id=6' +\
#                '&project_id=7' +\

response = urllib2.urlopen(api+u_query)
marinesed = json.load(response) 
marinesed=marinesed['success']['data']

#CONVERT TO STRUCTURED ARRAY
mynames = ['unit_id','project_id','col_id','col_area','min_thick','max_thick','strat_name_id','unit_name','b_age','t_age','pbdb_occurrences']
myformats=['i4','i4','i4','f4','f4','f4','i4','|S100','f4','f4','f4']

marinesed = json2struct(marinesed,mynames,myformats)



##%% SAVE RESPONSE
#f = open('macrostratNAm.pckl', 'wb')
#pickle.dump([sed,cols], f)
#f.close()

##%% LOAD FROM FILE
#f = open('macrostrat.pckl', 'rb')
#[sed,cols] = pickle.load(f)
#f.close()

#%% PARSE INTO TIME LIST

time_list=np.zeros(4000,dtype={'names': ['age', 'coverage', 'volume', 'volumeinterp'],'formats':['f4','i4','f4','f4']})
time_list['age']=np.arange(4000)+0.5

# Calculate average thickness in kilometers and duration in Myr
average_thick=np.mean((sed['min_thick']+sed['max_thick'])/2000);
average_duration=np.mean(sed['b_age']-sed['t_age']);
#average_time=np.mean(sed['b_age']-sed['t_age']);
#average_rate = average_thick/average_time;

nearby_window = 2000;

#COLLECT TIMESCALE AND TIMESERIES DATA
for i,age in enumerate(time_list):
    #SAMPLING METHOD FOR COVERAGE
    if age['age']<66:
        # Exclude recent alluvium by looking only at marine sediments where sed type is well-known
        idx = np.where( np.logical_and( marinesed['t_age']<age['age'], marinesed['b_age']>=age['age']) )
        tmp_units = marinesed[idx];
        idx = np.where( np.logical_and( marinesed['t_age']<age['age']+nearby_window, marinesed['b_age']>=age['age']-nearby_window) )
        nearby_units = marinesed[idx];
    else:
        idx = np.where( np.logical_and( sed['t_age']<age['age'], sed['b_age']>=age['age']) )
        tmp_units = sed[idx];
        idx = np.where( np.logical_and( sed['t_age']<age['age']+nearby_window, sed['b_age']>=age['age']-nearby_window) )
        nearby_units = sed[idx];
    

    time_list['coverage'][i]=len(set(tmp_units['col_id']));
    
    # Calculate average thickness in kilometers and duration in Myr
    tmp_thick=(tmp_units['min_thick']+tmp_units['max_thick'])/2000;
    tmp_duration=(tmp_units['b_age']-tmp_units['t_age']);
    time_list['volume'][i]=sum(tmp_thick*tmp_units['col_area']/tmp_duration);
    
    average_thick=np.mean((nearby_units['min_thick']+nearby_units['max_thick'])/2000);
    average_duration=np.mean(nearby_units['b_age']-nearby_units['t_age']);

    tmp_thick[tmp_thick==0]=average_thick/average_duration*tmp_duration[tmp_thick==0];
    time_list['volumeinterp'][i]=sum(tmp_thick*tmp_units['col_area']/tmp_duration);

    
#%%
    
######## FIGURE 1 - SEDIMENTARY UNITS NORMALIZED
fig1=plt.figure(1,figsize=(7.2,5.2625))
fig1.clf()
fig1.patch.set_facecolor('white')

#coverage = time_list['coverage']/float(len(cols));
coverage = ndimage.filters.gaussian_filter1d(time_list['coverage']/float(len(cols)),10,mode='reflect');

#TIME SERIES - EARTH HISTORY
ax1 = fig1.add_subplot(111)
ax1.plot(time_list['age'],
         coverage,
         '-',lw=2) # ,color=cnames['salmon']

#Put figure window on top of all other windows
#fig1.canvas.manager.window.raise_()

#FORMAT THE PLOT
ax1.axis([0, 4500, 0, 0.4])
ax1.axes.set_yticks([0,0.1,0.2,0.3,0.4])
ax1.grid(False)
ax1.invert_xaxis()  

#ADD GEOLOGICAL TIMESCALE
#addGTS(ax1,'international periods',0.05,0,text=False)
#addGTS(ax1,'international epochs',0.035,0.05,text=False)
#addGTS(ax1,'international eons',0.035,0.05,text=False)

#ADD X AND Y AXIS LABELS
ax1.axes.set_xlabel('Age (Ma)')
ax1.axes.set_ylabel('Covered fraction')

#%%
#ADD GLACIATIONS
glaciations = [720,635,580,2180];
for i in glaciations:
    ax1.plot([i,i],[0,1])


## ADD AVERAGES
boundaries = [720,2180,4000];
for i in range(1,np.size(boundaries)):
    segment = time_list['coverage'][(time_list['age']>boundaries[i-1]) & (time_list['age']<boundaries[i])]/float(len(cols))
    ax1.plot([boundaries[i-1],boundaries[i]],[segment.mean(),segment.mean()])




#%%
######## FIGURE 2 - SEDIMENTARY VOLUME 
fig2=plt.figure(2,figsize=(7.2,5.2625))
fig2.clf()
fig2.patch.set_facecolor('white')


from scipy import ndimage

#volume = ndimage.filters.gaussian_filter1d(time_list['volume'],10,mode='reflect'); # Not interpolated
volume = ndimage.filters.gaussian_filter1d(time_list['volumeinterp'],10,mode='reflect'); # Interpolated

#TIME SERIES - EARTH HISTORY
ax2 = fig2.add_subplot(111)
ax2.plot(time_list['age'],
         volume/1e6,
         '-',lw=2) # ,color=cnames['salmon']

#FORMAT THE PLOT
#ax2.axis([0, 4500, 0, 1.6])
#ax2.axes.set_yticks([0,0.4,0.8,1.2,1.6])
ax2.axis([0, 4500, 0, 0.27])
ax2.grid(False)
ax2.invert_xaxis()  

#ADD X AND Y AXIS LABELS
ax2.axes.set_xlabel('Age (Ma)')
ax2.axes.set_ylabel('Volume flux (km3/yr)')

#Put figure window on top of all other windows
#fig2.canvas.manager.window.raise_()

#ADD GLACIATIONS
glaciations = [720,635,580,2180];
for i in glaciations:
    ax2.plot([i,i],[0,1.6])

#gMyr3e21 = 3e21/1000/2700/1e9/1e6;
#ax2.plot([0,4500],[gMyr3e21,gMyr3e21])

#%%
######## FIGURE 3 - SEDIMENTARY VOLUME for last 1600 Myr, Scaled to global 
fig3=plt.figure(3,figsize=(7.2,5.2625))
fig3.clf()
fig3.patch.set_facecolor('white')


from scipy import ndimage

#volume_interpolated = time_list['volumeinterp'];
volume_interpolated = ndimage.filters.gaussian_filter1d(time_list['volumeinterp'],10,mode='reflect');

#TIME SERIES - EARTH HISTORY
ax3 = fig3.add_subplot(111)
ax3.plot(time_list['age'],
         volume_interpolated*6.07/1e6,
         '-',lw=2) # color=cnames['salmon']

#FORMAT THE PLOT
ax3.axis([0, 1600, 0, 1.6])
ax3.axes.set_yticks([0,0.4,0.8,1.2,1.6])
ax3.grid(False)
ax3.invert_xaxis()  

#ADD X AND Y AXIS LABELS
ax3.axes.set_xlabel('Age (Ma)')
ax3.axes.set_ylabel('Volume flux  (km3/yr)')

#Put figure window on top of all other windows
fig3.canvas.manager.window.raise_()

#ADD GLACIATIONS
glaciations = [720,635,580,541];
for i in glaciations:
    ax3.plot([i,i],[0,1.6])
    
gMyr3e21 = 3e21/1000/2700/1e9/1e6;
ax3.plot([0,4500],[gMyr3e21,gMyr3e21])


#%%
######## FIGURE 3 - SEDIMENTARY VOLUME for last 2500 Myr, Scaled to global 
fig4=plt.figure(4,figsize=(7.2,5.2625))
fig4.clf()
fig4.patch.set_facecolor('white')


from scipy import ndimage

volume_interpolated = ndimage.filters.gaussian_filter1d(time_list['volumeinterp'],10,mode='reflect');

#TIME SERIES - EARTH HISTORY
ax4 = fig4.add_subplot(111)
ax4.plot(time_list['age'],
         volume_interpolated*6.07/1e6,
         '-',lw=2) # color=cnames['salmon']

#FORMAT THE PLOT
ax4.axis([0, 2500, 0, 1.6])
ax4.axes.set_yticks([0,0.4,0.8,1.2,1.6])
ax4.grid(False)
ax4.invert_xaxis()  

#ADD X AND Y AXIS LABELS
ax4.axes.set_xlabel('Age (Ma)')
ax4.axes.set_ylabel('Volume flux  (km3/yr)')

#Put figure window on top of all other windows
fig4.canvas.manager.window.raise_()

#ADD GLACIATIONS
glaciations = [720,635,580,541];
for i in glaciations:
    ax4.plot([i,i],[0,1.6])
    
gMyr3e21 = 3e21/1000/2700/1e9/1e6;
ax4.plot([0,4500],[gMyr3e21,gMyr3e21])