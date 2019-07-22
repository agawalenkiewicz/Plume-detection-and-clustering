import numpy as np
import netCDF4 as nc
import os, sys, re
import matplotlib.pyplot as plt
import random
from scipy import integrate


"""
# DO NOT USE THIS FUNCTION!
# TCWV IS CALCULATED BY RTTOV ALREADY AND IN THE JACOBIANS FILE
def get_tcwv_from_prior(q, hpa, t):
	h = np.genfromtxt('.txt') #height of model levels
    R = 287.058
	p = hpa * 100.
    rho_air = p/(R*t) #density of dry air
    Q = rho_air*q  #kg/kg -> kg/m3
    tcwv = -integrate.simps(Q, h)
    return tcwv
"""
# function for reading in CCI SST and 
def read_in_apriori_state(prior_and_k_file):
	"""
	read in the SST, wv profile and TCWV from the example_k files
	"""
	### for sst ###############################
	f = open(prior_and_k_file)
	for i, line in enumerate(f):
		if i == 79:
			SST = np.float(line[45:57].strip())
	print 'SST from prior' , SST
	#### for water vapour #####################
	f = open(prior_and_k_file)
	p = []
	t = []
	q = []
	for i, line in enumerate(f):
		if 93 < i < 154:
			#lvl, p_elem, t_elem, q_elem = np.genfromtxt(line, delimiter='   ')
			p_elem = line[7:16].strip()
			t_elem = np.float(line[18:30].strip())
			q_elem = np.float(line[33:45].strip())
			p.append(p_elem)
			t.append(t_elem)
			q.append(q_elem)
	print 'humidity profile' , q
	##### for TCWV ###########################
	f = open(prior_and_k_file)
	for i, line in enumerate(f):
		if i == 154:
			TCWV = np.float(line[10:22].strip())
	print 'TCWV' , TCWV
	###for x_a ###############################
	x_a = np.array((SST , TCWV))
	print 'x_a matrix' , x_a
	return x_a, TCWV, q, SST
	
def sa_covariance(tcwv, SST):
	SST=8.
	S_a = np.diag((SST **2, (0.3*tcwv *(0.1+(7.5 - tcwv)/15.)) **2))
	print 's_a matrix' , S_a
	return S_a

def seps_covariance(sat):
	if sat == 'Landsat':
		s_eps = np.zeros((2,2))
		bt10_at280k = 0.053 #TIRS Noise-Equivalent-Change-in-Temperature (NEΔT) - Landsat's User's Handboook
		bt11_at280k = 0.059 #TIRS Noise-Equivalent-Change-in-Temperature (NEΔT) - Landsat's User's Handboook
		s_eps = np.diag((bt10_at280k**2, bt11_at280k ** 2))
	if sat == 'ASTER':
		"""
		For the TIR subsystem it is convenient to establish the subsystem noise
		requirement in terms of a noise equivalent delta temperature (NEDT). The
		subsystem requirement is that the NEDT be less than 0.3 K for a 300 K target.
		The accuracy requirement on the TIR subsystem is given for each of several
		brightness temperature ranges as follows: 200-240 K, 3 K; 240-270 K, 2 K; 270-
		340 K, 1 K; and 340-370 K, 2 K
		"""
		bt_for_all_channels = 0.3
		s_eps = np.zeros((5,5))
		s_eps = np.diag((bt_for_all_channels**2 , bt_for_all_channels**2, bt_for_all_channels**2, bt_for_all_channels**2, bt_for_all_channels**2))
	print 's_eps matrix' , s_eps
	return s_eps


def read_in_observations(bt_netcdf_file):
	"""
	Read in each .nc file for variable BT10 and BT11

	"""
	nc_file = nc.Dataset(bt_netcdf_file)
	BT10_obs = np.array(nc_file.variables['BT_band10'])
	BT11_obs = np.array(nc_file.variables['BT_band11'])
	# need to mask land and clouds on the bt arrays
	# y is a matrix of bt10 and bt11
	# O.E. happens for every pixel
	y = np.array((BT10_obs, BT11_obs))
	print 'Observed BT by Landsat for the whole scene' , y
	return y

# function for reading in simulated BT ####this reads in the whole array - need to read in just one day - change !!!
def read_in_simulations(sat, prior_and_k_file):
	f = open(prior_and_k_file)
	for i, line in enumerate(f):
		if i == 160:
			if sat == 'Landsat':
				BT10_sim = np.float(line[3:10].strip())
				BT11_sim = np.float(line[11:18].strip())
				
				Kx_a = np.array((BT10_sim, BT11_sim))
				print 'Kx_a matrix of simulated BT output by RTTOV fwd model' , Kx_a
			
			elif sat == 'ASTER':
				BT10_sim = np.float(line[3:10].strip())
				BT11_sim = np.float(line[11:18].strip())
				BT12_sim = np.float(line[19:26].strip())
				BT13_sim = np.float(line[27:34].strip())
				BT14_sim = np.float(line[35:42].strip())

				Kx_a = np.array((BT10_sim, BT11_sim, BT12_sim, BT13_sim, BT14_sim))
				print 'Kx_a matrix of simulated BT output by RTTOV fwd model' , Kx_a
	return Kx_a

def get_jacobian(sat, sst_jac_file, wv_jac_file, tcwv, q):
	"""
	sat - type of satellite
	sst_jac_file - jacobians for surface SST for different channels (on cluster: from jacobians_skin folder)
	wv_jac_file - jacobians for wv (on cluster: from jacobians_airtemp_wv folder)
	q - humididty profile from observations; needed for multiplication for the wv_jac_file to get the wv_tanlin
	tcwv - one value of TCWV, needed to calculate wv jacobians for each channel as wv_tanlin / tcwv
	"""
	wv_level_jacobian10 = []
	wv_level_jacobian11 = []

	if sat == 'Landsat':
		K = np.zeros((2,2))
		a = open(sst_jac_file)
		for j, line in enumerate(a):
			if j == 8:
				sst_jacobian10 = np.float(line[5:18].strip())
			elif j == 9:
				sst_jacobian11 = np.float(line[5:18].strip())
		
		print 'jacobians for each sst' , sst_jacobian10, sst_jacobian11
		
		wv_levels = np.array(q)
		b = open(wv_jac_file)
		for k, line in enumerate(b):
			if 250 < k < 311 :
				wv_jacobian10_elem = np.float(line[33:45].strip())
				wv_level_jacobian10.append(wv_jacobian10_elem)
			if 313 < k < 374:
				wv_jacobian11_elem = np.float(line[33:45].strip())
				wv_level_jacobian11.append(wv_jacobian11_elem)
		
		print 'jacobians for each wv level' , wv_level_jacobian10, wv_level_jacobian11
		
		wv_tanlin_10 = np.sum(np.array(wv_level_jacobian10) * wv_levels)
		wv_tanlin_11 = np.sum(np.array(wv_level_jacobian11) * wv_levels)
		
		print ' wv tanlin' , wv_tanlin_10, wv_tanlin_11
		
		wv_jacobian10 = wv_tanlin_10 / tcwv
		wv_jacobian11 = wv_tanlin_11 / tcwv

		print ' WV jacobians ' , wv_jacobian10, wv_jacobian11
		
		K = np.array(((sst_jacobian10, wv_jacobian10),(sst_jacobian11, wv_jacobian11)))
		
	if sat == 'ASTER':
		K = np.zeros((5,2))
		a = open(sst_jac_file)
		for j, line in enumerate(a):
			if j == 8:
				sst_jacobian10 = np.float(line[5:20].strip())
			elif j == 9:
				sst_jacobian11 = np.float(line[5:20].strip())
			elif j == 10:
				sst_jacobian12 = np.float(line[5:20].strip())
			elif j == 11:
				sst_jacobian13 = np.float(line[5:20].strip())
			elif j ==12:
				sst_jacobian14 = np.float(line[5:20].strip())
			
		print 'jacobians for each sst' , sst_jacobian10, sst_jacobian11, sst_jacobian12, sst_jacobian13, sst_jacobian14
		
		wv_levels = np.array(q)
		b = open(wv_jac_file)
		for k, line in enumerate(b):
			if 250 < k < 311 :
				wv_jacobian10_elem = np.float(line[33:49].strip())
				wv_level_jacobian10.append(wv_jacobian10_elem)
			if 313 < k < 374:
				wv_jacobian11_elem = np.float(line[33:49].strip())
				wv_level_jacobian11.append(wv_jacobian11_elem)
			if 376 < k < 431:
				wv_jacobian12_elem = np.float(line[33:49].strip())
				wv_level_jacobian12.append(wv_jacobian12_elem)
			if 433 < k < 494:
				wv_jacobian13_elem = np.float(line[33:49].strip())
				wv_level_jacobian13.append(wv_jacobian13_elem)
			if 496 < k < 557:
				wv_jacobian14_elem = np.float(line[33:49].strip())
				wv_level_jacobian14.append(wv_jacobian14_elem)

		
		print 'jacobians for each wv level' , wv_level_jacobian10, wv_level_jacobian11, wv_level_jacobian12, wv_level_jacobian13, wv_level_jacobian14
		
		wv_tanlin_10 = np.sum(np.array(wv_level_jacobian10) * wv_levels)
		wv_tanlin_11 = np.sum(np.array(wv_level_jacobian11) * wv_levels)
		wv_tanlin_12 = np.sum(np.array(wv_level_jacobian12) * wv_levels)
		wv_tanlin_13 = np.sum(np.array(wv_level_jacobian13) * wv_levels)
		wv_tanlin_14 = np.sum(np.array(wv_level_jacobian14) * wv_levels)
		
		print ' wv tanlin' , wv_tanlin_10, wv_tanlin_11, wv_tanlin_12, wv_tanlin_13, wv_tanlin_14
		
		wv_jacobian10 = wv_tanlin_10 / tcwv
		wv_jacobian11 = wv_tanlin_11 / tcwv
		wv_jacobian12 = wv_tanlin_12 / tcwv
		wv_jacobian13 = wv_tanlin_13 / tcwv
		wv_jacobian14 = wv_tanlin_14 / tcwv

		print ' WV jacobians ' , wv_jacobian10, wv_jacobian11, wv_jacobian12, wv_jacobian13, wv_jacobian14
		
		K = np.array(((sst_jacobian10, wv_jacobian10),(sst_jacobian11, wv_jacobian11),(sst_jacobian12, wv_jacobian12),(sst_jacobian13, wv_jacobian13),(sst_jacobian14, wv_jacobian14)))
		
		print 'K matrix ' , K
	return K

# Optimal Estimation (as in Rodgers, 2000)
def optimal_estimation(y, Kx_a, K, S_a, S_eps, x_a):
	"""
	A function reproduced from Rodgers (2000)
	'Inverse Methods for Atmospheric Sounding'
	Equation 4.3 on page 67
	x_a - water column vapour and SST
	S_a - 2x2 diagonal matrix with the uncertainties (squared) of the TCWV and SST
	S_eps - 2x2 diagonal matrix combining the measurement uncertainty and RTTOV model error (for RTTOV assume 0.15)
	K -  the matrix of partial derivatives for the reduced state vector
	"""
	eps = y - Kx_a 	#observed array - simulated single value for bt10 and single value for bt11
	#eps = np.transpose(eps)
	KS_a = np.ma.dot(K, S_a)
	KS_aK_T = np.ma.dot(KS_a, np.transpose(K))
	KS_aK_T_inv = np.linalg.inv(KS_aK_T + S_eps) #inverse
	S_aK_T = np.ma.dot(S_a, np.transpose(K))
	x_hat = x_a + np.ma.dot(np.ma.dot(S_aK_T, KS_aK_T_inv), eps)
	return x_hat


################################################
def oe_main(satellite_string, netcdf_file, prior_and_k, skin_k_file):
	##### make sure that file1, file2, sst_jacobian_file and wv_jacobian_file
	##### all have the same date !!!
	"""
	##### can be found in the filename what is the date
	prior_and_k = sys.argv[1] #sst, bt_sim, tcwv, q, q_jac
	netcdf_file = sys.argv[2] #.nc file with bt_obs 2d array
	sst_k_file = sys.argv[3] # skin sst_jac
	"""
	x_a, tcwv, q_list, sst_prior = read_in_apriori_state(prior_and_k)
	s_a = sa_covariance(tcwv, 8.0)
	s_eps = seps_covariance(satellite_string) #'Landsat' or 'ASTER'
	y = read_in_observations(netcdf_file)
	kx_a = read_in_simulations(satellite_string, prior_and_k)
	k = get_jacobian(satellite_string, skin_k_file, prior_and_k, tcwv, q_list)
	# loop through the pairs of observations (bt10 and bt11) for each pixel in the scene
	print 'bt observed' , y

	new_sst = np.zeros((401,401))
	for i in np.arange(0,401):
		for j in np.arange(0,401):
			"""
			bt10_pair = y[0,i,j]
			bt11_pair = y[1,i,j]
			y_pair = np.array((bt10_pair, bt11_pair))
			"""
			y_pair = y[:,i,j] #for i in np.arange(0,401) for j in np.arange(0,401)
			#y_pair - observations, 
			#s_eps - observation covariance matrix, 
			#s_a - prior covariance matrix, 
			# k - jacobians for sst and wv fpr both channels
			# kx_a - simulated bt for both channels
			# x_a - prior sst and tcwv
			new_sst[i,j], tcwv = optimal_estimation(y_pair, kx_a, k, s_a, s_eps, x_a)
			#print y_pair, new_sst
	print new_sst
	return new_sst, sst_prior
	
#################################################

#oe_main()
	
	
