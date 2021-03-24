import numpy as np
from bokeh.plotting import figure
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource,Range1d, CheckboxGroup, RangeSlider
from bokeh.layouts import layout, column
from scipy.optimize import least_squares

def Constitutive_law(Q,A,n,SigmaP,SigmaR,Temp,Rate):    
    # Constants 
    R = 8.31
    Temp = Temp + 273 # convert to Kelvin
    z = Rate*np.exp(Q/(R*Temp))
    yield_stress = 1e6*(SigmaR*np.arcsinh((z/A)**(1.0/n)) + SigmaP)
    return yield_stress

def activation_Constitutive_law_fit(i,activations,inputs,original_values):
    if activations[i] == 1:
        return inputs[i]
    else:
        return original_values[i]

def Constitutive_law_fit(inputs,exp_data,activations,original_values):
    
    # activate the inputs which have been checked
    SigmaR = activation_Constitutive_law_fit(0,activations,inputs,original_values)
    SigmaP = activation_Constitutive_law_fit(1,activations,inputs,original_values)
    Q = activation_Constitutive_law_fit(2,activations,inputs,original_values)
    A = activation_Constitutive_law_fit(3,activations,inputs,original_values)
    n = activation_Constitutive_law_fit(4,activations,inputs,original_values)
    
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
values_string = ('SigmaR = ' + "{0:.3}".format(ST_values[0]) + 
                 ', SigmaP = ' + "{0:.3}".format(ST_values[1]) + 
                     ', Q = ' + "{0:.2E}".format(ST_values[2]) + 
                     ', A = '+ "{0:.2E}".format(ST_values[3]) + 
                     ', n = ' + "{0:.3}".format(ST_values[4]))
fig.title.text = values_string

# deal with the case that the input value is not in the range of the sliders
def input_outside_slider(value,slider_min,slider_max):
    if slider_max < value:
        return slider_max
    elif slider_min > value:
        return slider_min
    else:
        return value
    
def activation_value(i,activations,fitted_values,original_values):
    if activations[i] == 1:
        return fitted_values.x[i]
    else:
        return original_values[i]
    

#create a sliders for the variables
def callback(attrname, old, new):
    
    # Make a list of the checkboxes that are active
    activations     = np.zeros(5)
    checked_ind     = checkbox_group.active
    for ind in checked_ind:
        activations[ind] = 1
        
    # Deal with the case where the slider values do not include the original value
    SigmaR_input = input_outside_slider(SigmaR,SigmaR_range_slider.value[0],SigmaR_range_slider.value[1])
    SigmaP_input = input_outside_slider(SigmaP,SigmaP_range_slider.value[0],SigmaP_range_slider.value[1])
    Q_input = input_outside_slider(Q,Q_range_slider.value[0],Q_range_slider.value[1])
    A_input = input_outside_slider(A,A_range_slider.value[0],A_range_slider.value[1])
    n_input = input_outside_slider(n,n_range_slider.value[0],n_range_slider.value[1])
    
    # Store the original values 
    original_values = [SigmaR,SigmaP,Q,A,n]
    
    # Set the inputs for optimisation
    inputs          = [SigmaR_input,SigmaP_input,Q_input,A_input,n_input]
    
    # Least squares fit on the experimental data
    fitted_values = least_squares(Constitutive_law_fit, inputs,
                                  bounds = [(SigmaR_range_slider.value[0],
                                             SigmaP_range_slider.value[0],
                                             Q_range_slider.value[0],
                                             A_range_slider.value[0],
                                             n_range_slider.value[0]),
                                            (SigmaR_range_slider.value[1],
                                             SigmaP_range_slider.value[1],
                                             Q_range_slider.value[1],
                                             A_range_slider.value[1],
                                             n_range_slider.value[1])],
                                  args=([corrected_exp_data,activations,original_values]))
    
    # Set the new values of the ST equation
    # Sigma R
    SigmaR_new = activation_value(0,activations,fitted_values,original_values)
    SigmaP_new = activation_value(1,activations,fitted_values,original_values)
    Q_new = activation_value(2,activations,fitted_values,original_values)
    A_new = activation_value(3,activations,fitted_values,original_values)
    n_new = activation_value(4,activations,fitted_values,original_values)
    
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
    
    # update the title
    values_string = ('SigmaR = ' + "{0:.3}".format(ST_values[0]) + 
                     ', SigmaP = ' + "{0:.3}".format(ST_values[1]) + 
                     ', Q = ' + "{0:.2E}".format(ST_values[2]) + 
                     ', A = '+ "{0:.2E}".format(ST_values[3]) + 
                     ', n = ' + "{0:.3}".format(ST_values[4]))
    fig.title.text = values_string

# Set checkboxes
LABELS = ["SigmaR", "SigmaP", "Q", "A", "n"]
checkbox_group = CheckboxGroup(labels=LABELS)
checkbox_group.on_change('active', callback)

# insert sliders for the max and min values
SigmaR_range_slider = RangeSlider(start=0, end=100, value=(0,50), step=.1, title="SigmaR")
SigmaP_range_slider = RangeSlider(start=0, end=100, value=(0,50), step=.1, title="SigmaP")
Q_range_slider = RangeSlider(start=0.5e5, end=3e5, value=(0.5e5,3e5), step=1e3, title="Q")
A_range_slider = RangeSlider(start=1.00e7, end=1.00e9, value=(1.00e7,1.00e9), step=1.00e7, title="A")
n_range_slider = RangeSlider(start=0, end=10, value=(3,6), step=.1, title="n")

# When slider values change call the callback function
SigmaR_range_slider.on_change('value', callback)
SigmaP_range_slider.on_change('value', callback)
Q_range_slider.on_change('value', callback)
A_range_slider.on_change('value', callback)
n_range_slider.on_change('value', callback)

# Set the layout of the sliders
sliders = column(SigmaR_range_slider,SigmaP_range_slider,Q_range_slider,A_range_slider,n_range_slider)

# Set layout of the figure
plot_layout = layout([[fig,[[checkbox_group],[sliders]]]])
curdoc().add_root(plot_layout)#serve it via "bokeh serve main.py --show --allow-websocket-origin=localhost:5006"