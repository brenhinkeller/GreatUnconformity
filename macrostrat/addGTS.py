# -*- coding: utf-8 -*-
"""
Created on Sun May 29 13:22:47 2016

@author: jhusson
"""

import requests
import matplotlib.patches as patches

#ADD GEOLOGICAL TIMESCALES TO FIGURE   
def addGTS(ax,name,height,ypos,text):
    q = 'https://macrostrat.org/api/defs/intervals?timescale=' + name
    success = requests.get(q)
    timescale = success.json()
    timescale=timescale['success']['data']

    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    
    for age in timescale:
        rect = patches.Rectangle(((xlim[0]-age['b_age'])/(xlim[0]-xlim[1]), ypos), (age['b_age']-age['t_age'])/(xlim[0]-xlim[1]), height, transform=ax.transAxes,
                              facecolor=age['color'],zorder=100)
        ax.add_patch(rect)
        
        if age['abbrev'] and age['int_id']!=421 and age['b_age']<= xlim[0] and text and age['t_age']>= xlim[1]:
            ax.text((age['b_age']-age['t_age'])/2+age['t_age'], 
                    (ypos+height/1.5)*(ylim[1]-ylim[0])+ylim[0],
                    age['abbrev'],zorder=101,ha='center',va='center',fontsize=10)
