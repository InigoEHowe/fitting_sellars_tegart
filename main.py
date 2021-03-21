import numpy as np
from bokeh.plotting import figure
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource,Range1d, CheckboxGroup, PreText
from bokeh.layouts import layout
from scipy.optimize import least_squares
from decimal import Decimal

def Constitutive_law(Q,A,n,SigmaP,SigmaR,Temp,Rate):    
    # Constants 
    R = 8.31
    Temp = Temp + 273 # convert to Kelvin
    z = Rate*np.exp(Q/(R*Temp))
    yield_stress = 1e6*(SigmaR*np.arcsinh((z/A)**(1.0/n)) + SigmaP)
    return yield_stress

def Constitutive_law_fit(inputs,exp_data,activations,original_values):
    
    # Sigma R
    if activations[0] == 1:
        SigmaR = inputs[0]
    else:
        SigmaR = original_values[0]
        
    # Sigma P
    if activations[1] == 1:
        SigmaP = inputs[1]
    else:
        SigmaP = original_values[1]
        
    # Q
    if activations[2] == 1:
        Q = inputs[2]
    else:
        Q = original_values[2]
    
    # A
    if activations[3] == 1:
        A = inputs[3]
    else:
        A = original_values[3]
        
    # n
    if activations[4] == 1:
        n = inputs[4]
    else:
        n = original_values[4]
    
    # Constants 
    R = 8.31
    
    strains = np.arange(0.05,0.7,0.05)
    temps = [490,520,560]
    rates = [0.1,1,10]
    
    yield_PRO3 = np.zeros(len(strains)*len(temps)*len(rates))
    count = 0
    for Temp in temps:
        Temp = Temp + 273 # convert to Kelvin
        for Rate in rates:
            for Strain in strains:
                z = Rate*np.exp(Q/(R*Temp))
                yield_PRO3[count] = (SigmaR*np.arcsinh((z/A)**(1.0/n)) + SigmaP)
                count = count + 1
         
    diff = abs(yield_PRO3 - exp_data[:,0])
    return diff

# Variables
Q = 1.56e5
A = 7.00e8
n = 5.0
SigmaP = 3.19
SigmaR = 17.39

# Load in the experimental data
corrected_exp_data = np.load('data/corrected_exp_data.npy')

#### Produce the plot
rates = [0.1,1,10] 
temp_range = np.arange(400,610,10)
fig=figure() # create a figure of log(sigma) vs T
curdoc().theme = 'light_minimal' # figure theme
colours = ['Red','Green','Blue']
count = 0

## Plotting experimental data
for rate in rates:
    # Experiment
    exp_data = corrected_exp_data[abs(corrected_exp_data[:,3]-rate)<10**(-6)]
    fig.circle(exp_data[:,2],np.log10(exp_data[:,0]),color=colours[count],legend_label=str(rate)+'s-1')
    count = count + 1

## Plotting Sellars-Tegart for the 3 rates
SellarsTegart_01 = np.zeros(len(temp_range))
SellarsTegart_1 = np.zeros(len(temp_range))
SellarsTegart_10 = np.zeros(len(temp_range))
for i in range(len(temp_range)):
    SellarsTegart_01[i] = np.log10(Constitutive_law(Q,A,n,SigmaP,SigmaR,temp_range[i],0.1)*10**(-6)) 
    SellarsTegart_1[i]  = np.log10(Constitutive_law(Q,A,n,SigmaP,SigmaR,temp_range[i],1  )*10**(-6)) 
    SellarsTegart_10[i] = np.log10(Constitutive_law(Q,A,n,SigmaP,SigmaR,temp_range[i],10 )*10**(-6)) 
    
source_01=ColumnDataSource(dict(x=temp_range, y=SellarsTegart_01))
source_1=ColumnDataSource(dict(x=temp_range, y=SellarsTegart_1))
source_10=ColumnDataSource(dict(x=temp_range, y=SellarsTegart_10))

# Baseline
fig.line(temp_range, SellarsTegart_01,color='Red',line_width=2,alpha=0.2)
fig.line(temp_range, SellarsTegart_1,color='Green',line_width=2,alpha=0.2)
fig.line(temp_range, SellarsTegart_10,color='Blue',line_width=2,alpha=0.2)

# Updated
fig.line('x','y',source=source_01,color='Red',line_width=2)
fig.line('x','y',source=source_1,color='Green',line_width=2)
fig.line('x','y',source=source_10,color='Blue',line_width=2)
    
## Set plot aesthetics
fig.xaxis[0].axis_label = 'T / Degrees C'
fig.yaxis[0].axis_label = 'log(Sigma / MPa)'
fig.x_range=Range1d(400, 600)
fig.y_range=Range1d(1, 2)
fig.xaxis.axis_label_text_font_size = "15pt"
fig.xaxis.major_label_text_font_size = "15pt"
fig.yaxis.axis_label_text_font_size = "15pt"
fig.yaxis.major_label_text_font_size = "15pt"
fig.title.text_font_size = '10pt'

# Set plot title
ST_values = [SigmaR,SigmaP,Q,A,n]
values_string = ('SigmaR = ' + "{0:.2E}".format(ST_values[0]) + 
                     'MPa, SigmaP = ' + "{0:.2E}".format(ST_values[1]) + 
                     'MPa, Q = ' + "{0:.2E}".format(ST_values[2]) + 
                     'MPa, A = '+ "{0:.2E}".format(ST_values[3]) + 
                     'MPa, n = ' + "{0:.2E}".format(ST_values[4]))
fig.title.text = values_string

#create a sliders for the variables
def callback(attrname, old, new):
    
    # Set up inputs
    inputs      = [SigmaR,SigmaP,Q,A,n]
    activations = np.zeros(5)
    checked_ind = checkbox_group.active
    
    # make a list of the checkboxs that are active
    for ind in checked_ind:
        activations[ind] = 1
    
    # Least squares fit on the experimental data
    fitted_values = least_squares(Constitutive_law_fit, inputs,
                                  bounds = [(0,0,0.5e5,1.00e7,3),(50,50,3e5,1.00e9,10)],
                                  args=([corrected_exp_data,activations,inputs]))
    
    # Sigma R
    if activations[0] == 1:
        SigmaR_new = fitted_values.x[0]
    else:
        SigmaR_new = SigmaR
        
    # Sigma P
    if activations[1] == 1:
        SigmaP_new = fitted_values.x[1]
    else:
        SigmaP_new = SigmaP
        
    # Q
    if activations[2] == 1:
        Q_new = fitted_values.x[2]
    else:
        Q_new = Q
    
    # A
    if activations[3] == 1:
        A_new = fitted_values.x[3]
    else:
        A_new = A
        
    # n
    if activations[4] == 1:
        n_new = fitted_values.x[4]
    else:
        n_new = n
    
    SellarsTegart_01 = np.zeros(len(temp_range))
    SellarsTegart_1  = np.zeros(len(temp_range))
    SellarsTegart_10 = np.zeros(len(temp_range))
    for i in range(len(temp_range)):
        SellarsTegart_01[i] = np.log10(Constitutive_law(Q_new,A_new,n_new
                                                        ,SigmaP_new,SigmaR_new,temp_range[i],0.1)*10**(-6))
        SellarsTegart_1[i]  = np.log10(Constitutive_law(Q_new,A_new,n_new
                                                        ,SigmaP_new,SigmaR_new,temp_range[i],1  )*10**(-6))
        SellarsTegart_10[i] = np.log10(Constitutive_law(Q_new,A_new,n_new
                                                        ,SigmaP_new,SigmaR_new,temp_range[i],10)*10**(-6))
        
    source_01.data = dict(x=temp_range, y=SellarsTegart_01)
    source_1.data  = dict(x=temp_range, y=SellarsTegart_1) 
    source_10.data = dict(x=temp_range, y=SellarsTegart_10) 
    
    # update text box
    ST_values = [SigmaR_new,SigmaP_new,Q_new,A_new,n_new]
    
    values_string = ('SigmaR = ' + "{0:.2E}".format(ST_values[0]) + 
                     'MPa, SigmaP = ' + "{0:.2E}".format(ST_values[1]) + 
                     'MPa, Q = ' + "{0:.2E}".format(ST_values[2]) + 
                     'MPa, A = '+ "{0:.2E}".format(ST_values[3]) + 
                     'MPa, n = ' + "{0:.2E}".format(ST_values[4]))
    fig.title.text = values_string
    
LABELS = ["SigmaR", "SigmaP", "Q", "A", "n"]
checkbox_group = CheckboxGroup(labels=LABELS, active=[0, 0, 0, 0, 0])
checkbox_group.on_change('active', callback)

plot_layout = layout([[fig,checkbox_group]])
curdoc().add_root(plot_layout)#serve it via "bokeh serve main.py --show --allow-websocket-origin=*"