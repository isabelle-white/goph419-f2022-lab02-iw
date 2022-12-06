def main():
    
    """
    Created on Thu Dec  1 10:35:01 2022
    
    In this script I load the air and water vs temperature data stored in text files.
    This script then plots the spline function found using 100 evenly spaced temperature values within the range of the 
    origional data at orders 1, 2, and 3 for both air density and water density data. 
    
    The origional data points are plotted on the same graphs as red dots. 
    
    Parameters: 
    x_air: temperature data from the air_density_vs_temp_eng_toolbox txt file 
    y_air: air density data from the air_density_vs_temp_eng_toolbox txt file
    min_air_temp: minimum air temperature in the file 
    max_air_temp:maximum air temperature in the file 
    f_air1 = 1st order spline function 
    y_air1 = air densities using the spline function, order =1 
    f_air2 = 2nd order spline function 
    y_air2 = air densities using the spline function, order =2
    f_air3 = 3rd order spline function 
    y_air3 = air densities using the spline function, order =3 
    
    x_water:temperature data from the water_density_vs_temp_usgs txt file 
    y_water: water denisty data from the  water_density_vs_temp_usgs txt file
    min_water_temp: minimum water temperature in the file 
    max_water_temp:maximum water temperature in the file 
    f_water1 = 1st order spline function 
    y_water1 = water densities using the spline function, order =1 
    f_water2 = 2nd order spline function 
    y_water2 = water densities using the spline function, order =2
    f_water3 = 3rd order spline function 
    y_water3 = water densities using the spline function, order =3
    
    Inputs:
    air_density_vs_temp_eng_toolbox.txt
    water_density_vs_temp_usgs.txt
    
    Outputs: 
    3x2 subplot of order 1, 2, and 3 spline functions using the loaded data 
    
    
    @author: izzywhite
    """
    
    
    import numpy as np
    import matplotlib.pyplot as plt
    from spline_function_new import spline_function
    
    #set the subplot shape as 2 x 3
    fig, ax = plt.subplots(2, 3)
    fig.tight_layout(pad=3)
    
    #loading air density data 
    dens_data_air = np.loadtxt('air_density_vs_temp_eng_toolbox.txt')
    x_air = dens_data_air[:, 0]
    y_air = dens_data_air[:, 1]
    
    #getting 100 evenly spaced temperature values
    min_air_temp = np.min(x_air)
    max_air_temp = np.max(x_air)
    x_air_temp = np.linspace(min_air_temp, max_air_temp, 100)
    
    #get y_air estimates using 1st order spline function
    f_air1 = spline_function(x_air, y_air, order=1)
    y_air1 = f_air1(x_air_temp)
    
    #plot graph
    ax[0][0].plot(x_air_temp, y_air1)
    ax[0][0].plot(x_air, y_air, 'ro', markersize='2')
    ax[0][0].set_xlabel('Temperature')
    ax[0][0].set_ylabel('Air density')
    ax[0][0].set_title('Order=1')
    
    
    #get y_air estimates when order = 2 for the spline function
    f_air2 = spline_function(x_air, y_air, order=2)
    y_air2 = f_air2(x_air_temp)
    
    #plot the graph
    ax[0][1].plot(x_air_temp, y_air2)
    ax[0][1].plot(x_air, y_air, 'ro', markersize='2')
    ax[0][1].set_xlabel('Temperature')
    ax[0][1].set_ylabel('Air density')
    ax[0][1].set_title('Order=2')
    
    #get y_air estimates when order=3 for the spline function 
    f_air3 = spline_function(x_air, y_air, order=3)
    y_air3 = f_air3(x_air_temp)
    
    #plot the graph 
    ax[0][2].plot(x_air_temp, y_air3)
    ax[0][2].plot(x_air, y_air, 'ro', markersize = '2')
    ax[0][2].set_xlabel('Temperature')
    ax[0][2].set_ylabel('Air density')
    ax[0][2].set_title('Order=3')
    
    
    #loading water density data 
    dens_data_water = np.loadtxt('water_density_vs_temp_usgs.txt')
    x_water = dens_data_water[:, 0]
    y_water = dens_data_water[:, 1]
    
    #getting 100 evenly spaced temperature values
    min_water_temp = np.min(x_water)
    max_water_temp = np.max(x_water)
    x_water_temp = np.linspace(min_water_temp, max_water_temp, 100)
    
    #get y_water estimates using 1st order spline function
    f_water1 = spline_function(x_water, y_water, order=1)
    y_water1 = f_water1(x_water_temp)
    
    #plot graph
    ax[1][0].plot(x_water_temp, y_water1)
    ax[1][0].plot(x_water, y_water, 'ro', markersize='2')
    ax[1][0].set_xlabel('Temperature')
    ax[1][0].set_ylabel('Water density')
    ax[1][0].set_title('Order=1')
    
    
    #get y_water estimates using 2nd order spline function
    f_water2 = spline_function(x_water, y_water, order=2)
    y_water2 = f_water2(x_water_temp)
    
    #plot graph
    ax[1][1].plot(x_water_temp, y_water2)
    ax[1][1].plot(x_water, y_water, 'ro', markersize='2')
    ax[1][1].set_xlabel('Temperature')
    ax[1][1].set_ylabel('Water density')
    ax[1][1].set_title('Order=2')
    
    #get y_water estimates using 1st order spline function
    f_water3 = spline_function(x_water, y_water, order=3)
    y_water3 = f_water3(x_water_temp)
    
    #plot graph
    ax[1][2].plot(x_water_temp, y_water3)
    ax[1][2].plot(x_water, y_water, 'ro', markersize='2')
    ax[1][2].set_xlabel('Temperature')
    ax[1][2].set_ylabel('Water density')
    ax[1][2].set_title('Order=3')
    
    
    plt.savefig("Density of air and water vs temperature")
    
    if __name__ == "__main__":
        main()
