# -*- coding: utf-8 -*-
"""
Created on Sun May 29 13:17:48 2016

@author: jhusson
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#LIST OF COLORS
cnames={}
for name, hex in matplotlib.colors.cnames.iteritems():
   cnames[name]=hex.encode('ascii', 'ignore')
   
#GRAPHICALLY SHOW AVAILABLE PLOTTING COLORS  
def show_color(cnames):
    #MAKE A FIGURE AND AXIS
    fig=plt.figure(figsize=(24,6))
    fig.patch.set_facecolor('white')
    ax=fig.add_subplot(111)
    
    #EXPAND THE FIGURE
    pos1 = ax.get_position()
    pos2 = [pos1.x0-0.12, pos1.y0 + 0.035,  pos1.width+0.2, pos1.height]
    ax.set_position(pos2)
    
    #LOOP THROUGH COLOR DICTIONARY
    for i,c in enumerate(cnames):
        #MAKE A COLORED RECTANGLE
        rect = patches.Rectangle((float(i)/len(cnames),0.05), 1.0/len(cnames), 0.85,
                              facecolor=cnames[c],edgecolor=cnames[c])
        ax.add_patch(rect)
        
        #ADD LABELS (TOP AND BOTTOM)
        if i % 2==0:
            ax.text(float(i)/len(cnames), 
                    0.91,
                    c,rotation=45,ha='left',va='bottom',fontsize=10,color=cnames['black'])
        else:   
            ax.text(float(i)/len(cnames), 
                   0.049,
                    c,rotation=-45,ha='left',va='top',fontsize=10,color=cnames['slategrey'])
    
    #ADD A BLACK BAR TO ASSOCIATE TOP LABELS WITH RECTANGLE TOPS                    
    rect = patches.Rectangle((0,0.89), 1.0, 0.01,
                              facecolor=cnames['black'],edgecolor=cnames['black'])
    ax.add_patch(rect)
    
    #ADD A GREY BAR TO ASSOCIATE TOP LABELS WITH RECTANGLE TOPS                    
    rect = patches.Rectangle((0,0.05), 1.0, 0.01,
                              facecolor=cnames['slategrey'],edgecolor=cnames['slategrey'])
    ax.add_patch(rect)

    #FORMAT FIGURE
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.axis('off')