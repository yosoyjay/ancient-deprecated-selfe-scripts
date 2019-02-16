from numpy import *
from gr import *
import model_util
from model_reader import mreader
from geometry import *
import pdb
from shape_procs import *
from pylab import *
import matplotlib

# meccs is the final processing class for the python Model-based Estuary Characterization and Classification
# System.  It accepts data from the mreader (for raw data) and shape_procs (shape processing) class, which
# extract raw model data (mreader), and process it over given shapes
# meccs is responsible for performing time and space averaging, and calculating final quantities such as salt
# intrusion, Hansen & Rattray numbers, etc

class meccs:
  class _struct:
    pass

  def __init__(self):
    pass

  def get_salinity_intrusion(self, shp, salt_data, dt, time_interval, threshold):
    #shp in the form of x, y pts, assumes its a transect: determines mouth vs. head by lowest time-averaged maximum salt
    #salt data is a list of a 2D array in max_salt (through depth): [pt data, timesteps]
    tstep_int = ceil(time_interval/dt)  #keeping it simple and slightly innacurate
    nsteps = int(salt_data.shape[1]/tstep_int)

    ret = zeros(nsteps, float)

    #make sure the mouth is at the start of the salt_data array
    if mean(salt_data[0,:]) < mean(salt_data[-1,:]):
      tmp = salt_data
      tmpshp = shp
      for i in range(salt_data.shape[0]):
        salt_data[-i-1,:] = tmp[i,:]
        shp[-i-1,:] = tmpshp[i,:]

    seg_dist = distance(shp[0:-1,0], shp[0:-1,1], shp[1:,0], shp[1:,1])
    cum_dist = append(array(0, float), cumsum(seg_dist))

    for i in range(nsteps):
      salty = where(mean(salt_data[:,i*tstep_int:(i+1)*tstep_int], 1) >= threshold)
      if len(salty) != 0:
        ret[i] = cum_dist[max(salty[0])]
      else:
        ret[i] = 0

    return ret      

  def get_Hansen_Rattray(self, surf_salt, bott_salt, da_salt, dsurf_vel_data, mag_vel_data, dt, time_interval):
    tstep_int = ceil(time_interval/dt)  #keeping it simple and slightly innacurate
    nsteps = int(surf_salt.shape[0]/tstep_int)

    ret = zeros((2, nsteps), float)  #x and y axes are returned

    ds_das = (bott_salt - surf_salt)/da_salt
    su_um = -dsurf_vel_data/mag_vel_data #assumes (-) values are out of the estuary

    for i in range(nsteps):
      ret[1,i] = mean(ds_das[i*tstep_int:(i+1)*tstep_int])
      #      ret[0,i] = mean(su_um[i*tstep_int:(i+1)*tstep_int])/mean(mag_vel_data[
      ret[0,i] = mean(-dsurf_vel_data[i*tstep_int:(i+1)*tstep_int])/mean(mag_vel_data[i*tstep_int:(i+1)*tstep_int])
    return ret

  def get_time_avg_data(self, data, dt, time_interval):
    tstep_int = floor(float(time_interval)/dt)  #keeping it simple and slightly innacurate
    nsteps = int(data.shape[0]/tstep_int)

    ret = zeros(nsteps, float)  
    for i in range(nsteps):
      ret[i] = mean(data[i*tstep_int:(i+1)*tstep_int])

    return ret

  def get_amplitude(self, data, dt, time_interval):
    tstep_int = floor(time_interval/dt)  #keeping it simple and slightly innacurate
    nsteps = int(data.shape[0]/tstep_int)

    ret = zeros(nsteps, float)  
    for i in range(nsteps):
      ret[i] = max(data[i*tstep_int:(i+1)*tstep_int]) - min(data[i*tstep_int:(i+1)*tstep_int])

    return ret

  def get_Prandle():
    pass

  def save_xy_graph(self, xdata, ydata, xlabel, ylabel, title):
    pass

  ################# plotting functions ################
#  def hnr_plot(self, su_u, ds_s, name, legend_str):
#    colors=['b', 'r', 'g', 'y', 'c']
#    markers=['o', '*', '+', 'x', 's', 'd', 'h']
#
#    pdb.set_trace()
#
#    figure1=self.hnr_figure_base();
#
#    for i in range(su_u.shape[1]):
#      matplotlib.scatter(su_u[:,i], ds_s[:,i], 'MarkerFaceColor', colors[mod(i,5)], \
#              'MarkerEdgeColor', colors[mod(i,5)], 'Marker', markers[mod(i,len(markers))])
#
#    title(name)
#    if su_u.shape<2:
#      legend('off')
#    else:
#      legend (legend_str)
#
#    xlabel('2 Layer Circulation (U_s/U_f)')
#    ylabel('Stratification (dS/S_0)')
#    legend(legend_str)
#
#    astr = name+".png"
#    savefig(astr)
#
#
#  def hnr_figure_base(self):
#    xlm = [1, 100000]
#    ylm = [0.001, 10]
#
#    pdb.set_trace()
#
#    fig = matplotlib.Figure()
##    fig = figure('Color',[1, 1, 1])
#    set(fig, 'Position', [20, 50, 640, 480])
#    hold
#    xlim('manual')
#    set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', xlm, 'YLim', ylm, 'XLimMode', 'manual')
#
#    return fig

#  def hnr_figure_annotations(self, do_hnr_original):
#    # assumes current figure and axis
#
#    colr = [0.95, 0.95, 0.95]
#    fsize = 10
#    fweight = 'bold'
#
#    xlm = [1, 100000]
#    ylm = [0.001, 10]
#
#    plot(xlm, [10.^-1, 10.^-1], 'Color', 'k', 'LineWidth', 2)
#
#    #1 (1b)
#    xcoords = [1.085, 10.^0.27]
#    ycoords = [10.^0.1, 10.^-1]
#    x = [xcoords(1), xcoords(1), xcoords(2), xcoords(2)]
#    y = [ycoords(1), ycoords(2), ycoords(2), ycoords(1)]
#    patch(x, y, colr, 'FaceAlpha', 0.3, 'LineStyle', 'none')
#    text(10.^0.05, 10.^-0.4, '1b', 'FontSize', fsize, 'FontWeight', fweight)
#
#    #2 (1a)
#    ycoords = [10.^-1, 10.^-2.33]
#    x = [xcoords(1), xcoords(1), xcoords(2), xcoords(2)]
#    y = [ycoords(1), ycoords(2), ycoords(2), ycoords(1)]
#    patch(x, y, colr, 'FaceAlpha', 0.3, 'LineStyle', 'none')
#    text(10.^0.05, 10.^-1.6, '1a', 'FontSize', fsize, 'FontWeight', fweight)
#
#    #3 (2b)
#    xcoords = [10.^0.35, 10.^0.8, 10.^1.7];
#    ycoords = [10.^0.1, 10.^-0.1, 10.^-1];
#    x = [xcoords(1), xcoords(2), xcoords(3), xcoords(1)];
#    y = [ycoords(1), ycoords(2), ycoords(3), ycoords(3)];
#    patch(x, y, colr, 'FaceAlpha', 0.3, 'LineStyle', 'none');
#    text(10.^0.65, 10.^-0.6, '2b', 'FontSize', fsize, 'FontWeight', fweight);
#    
#    %4 (2a)
#    slope = (log10(ycoords(2)) - log10(ycoords(3)))./(log10(xcoords(2)) - log10(xcoords(3)));
#    ycoords = [10.^-1, 10.^-2.33];
#    x = [xcoords(3), 10.^(log10(xcoords(3))+(log10(ycoords(2))-log10(ycoords(1)))./slope), xcoords(1), xcoords(1)];
#    y = [ycoords(1), ycoords(2), ycoords(2), ycoords(1)];
#    patch(x, y, colr, 'FaceAlpha', 0.3, 'LineStyle', 'none');
#    text(10.^0.95, 10.^-1.9, '2a', 'FontSize', fsize, 'FontWeight', fweight);
#    
#    xcoord = 10.^(1 + (-1-(-0.1))./slope);
#    %5 (3b)
#    vpow = [0.22, 1.6];
#    p = 3.1;%0.9;
#    xpow = (2:-0.01:1);
#    ypow = -((4.*p).*(xpow-vpow(1))).^0.14 + vpow(2);
#    x = [10.^1, xcoord, xlm(2), xlm(2), 10.^xpow];
#    y = [10.^-0.15, 10.^-1, 10.^-1, 10.^0.05, 10.^ypow];
#    patch(x, y, colr, 'FaceAlpha', 0.3, 'LineStyle', 'none');
#    text(10.^3, 10.^-0.4, '3b', 'FontSize', fsize, 'FontWeight', fweight);
#    
#    % 6 (3a)
#    x = [xcoord, 10.^(log10(xcoord) + (log10(ycoords(2))-(-1))./slope), xlm(2), xlm(2)];
#    y = [10.^-1, ycoords(2), ycoords(2), 10.^-1];
#    patch(x, y, colr, 'FaceAlpha', 0.3, 'LineStyle', 'none');
#    text(10.^4, 10.^-1.6, '3a', 'FontSize', fsize, 'FontWeight', fweight);
#    
#    % 7 (4)
#    xpow = (0.95:-0.01:0.22);
#    ypow = -((4.*p).*(xpow-vpow(1))).^0.14 + vpow(2);
#    x = [1.085, 1.085, xcoords(1), xcoords(2), 10.^0.95, 10.^xpow];
#    y = [9.99, 10.^0.5, 10.^0.15, 10.^-0.05, 10.^-0.08, 10.^ypow];
#    patch(x, y, colr, 'FaceAlpha', 0.3, 'LineStyle', 'none');
#    text(10.^0.09, 10.^0.8, '4', 'FontSize', fsize, 'FontWeight', fweight);
#    
#    % black
#    xpow = (0.22:0.01:2);
#    ypow = -((4.*p).*(xpow-vpow(1))).^0.14 + vpow(2);
#    x = [10.^xpow, xlm(2), xlm(2)];
#    y = [10.^ypow, 10.^0.05, 10];
#    patch(x, y, 'k', 'FaceAlpha', 0.3, 'LineStyle', 'none');
#    
#    % vpow = [0.22, 1.6];
#    % p = 3.1;%0.9;
#    % xpow = (0.22:0.01:5);
#    % ypow = -((4.*p).*(xpow-vpow(1))).^0.14 + vpow(2);
#    % plot(10.^xpow, 10.^ypow, 'r', 'LineWidth', 2);
#
#if (exist('do_hnr_original')~=0)   %matlab doesn't seem to do logical short circuit...
#    if (do_hnr_original)
#        label_offset = 0.04
#        
#        xpow = [0.25, 0.70];
#        ypow = [0.05, -0.30];
#        plot(10.^xpow, 10.^ypow, ':k', 'Marker', 'o', 'MarkerFaceColor', 'k');
#        text(10.^(xpow(1)+label_offset), 10.^(ypow(1)+label_offset), 'C_h', 'FontSize', fsize, 'FontWeight', fweight);
#        text(10.^(xpow(2)+label_offset), 10.^(ypow(2)+label_offset), 'C_l', 'FontSize', fsize, 'FontWeight', fweight);
#        
#        xpow = [0.21, 0.83];
#        ypow = [0.3, 0.05];
#        plot(10.^xpow, 10.^ypow, ':k', 'Marker', 'o', 'MarkerFaceColor', 'k');
#        text(10.^(xpow(1)+label_offset), 10.^(ypow(1)+label_offset), 'M_h', 'FontSize', fsize, 'FontWeight', fweight);
#        text(10.^(xpow(2)+label_offset), 10.^(ypow(2)+label_offset), 'M_l', 'FontSize', fsize, 'FontWeight', fweight);
#        
#        xpow = [1.12, 1.4];
#        ypow = [-0.6, -0.75];
#        plot(10.^xpow, 10.^ypow, 'k', 'Marker', 'o', 'MarkerFaceColor', 'k');
#        text(10.^(xpow(1)+label_offset), 10.^(ypow(1)+label_offset), 'J_1_7', 'FontSize', fsize, 'FontWeight', fweight);
#        text(10.^(xpow(2)+label_offset), 10.^(ypow(2)+label_offset), 'J_1_1', 'FontSize', fsize, 'FontWeight', fweight);
#        
#        xpow = [1.8];
#        ypow = [-1.3];
#        plot(10.^xpow, 10.^ypow, 'k', 'Marker', 'o', 'MarkerFaceColor', 'k');
#        text(10.^(xpow(1)+label_offset), 10.^(ypow(1)+label_offset), 'NM', 'FontSize', fsize, 'FontWeight', fweight);
#        
#        xpow = [2.05, 2.55];
#        ypow = [-0.9, -1.4];
#        plot(10.^xpow, 10.^ypow, ':k');
#        text(10.^((xpow(1)+xpow(2))./2), 10.^((ypow(1)+ypow(2))./2), 'JF', 'FontSize', fsize, 'FontWeight', fweight);
#        
#        xpow = [3.8, 4.55];
#        ypow = [-0.25, -0.93];
#        plot(10.^xpow, 10.^ypow, ':k');
#        text(10.^((xpow(1)+xpow(2))./2), 10.^((ypow(1)+ypow(2))./2), 'S_h', 'FontSize', fsize, 'FontWeight', fweight);
#
#        xpow = [4.5, 5.1];
#        ypow = [-0.8, -1.35];
#        plot(10.^xpow, 10.^ypow, ':k');
#        text(10.^((xpow(1)+xpow(2))./2), 10.^((ypow(1)+ypow(2))./2-0.15), 'S_h', 'FontSize', fsize, 'FontWeight', fweight);
#    end
#end
