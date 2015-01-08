## Should provide all the necessary methods for building a light curve
##based on http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/data/pyLikelihood/light_curve.py and https://github.com/kialio/LATAnalysisScripts
#Author: Matthew Brady
#Last modified: 8/8/2014


import os
import time
import multiprocessing as mp
import ConfigParser 

from gt_apps import *
from BinnedAnalysis import *
from UnbinnedAnalysis import *
import pyLikelihood as pyLike
import IntegralUpperLimit as IUL
import pyfits as fits

import numpy as np


def bin_shell(infobinandself): #This method is a shell for the run_bin() method, which is required for this arrangment of multiprocessing 
	infobin, lc = infobinandself
	lc.run_bin(infobin)

class Light_Curve(object):
	
	events_list = 'events.txt'
	debug = False	

	def __init__(self, config_path):
		self.config_path = config_path
		config = ConfigParser.RawConfigParser()
		config.read(config_path)

		self.base_name = config.get('Object', 'Base Name')
		self.datastart = config.getint('Object', 'Data Start')
		self.datastop = config.getint('Object', 'Data Stop')
		self.ra = config.getfloat('Object', 'ra')
		self.dec = config.getfloat('Object', 'dec')
		self.main_source = config.get('Object', 'Object 2FGL Name')
		self.source = config.getboolean('Object', 'Use Sourceclass Photons')
		self.binned = config.getboolean('Object', 'Use Binned Analysis')
		self.tbin = config.getint('Object', 'Seconds per Bin')
		self.ul_confidence = config.getfloat('Object', 'Upperlimit Confidence')
		self.data_path = config.get('Object', 'Data Folder')#

		self.fgl_path = config.get('General', 'FGL Catalog Path')#
		self.catalog_file = config.get('General', 'FGL Catalog File')
		self.makexml_path = config.get('General', 'make2FGLxml Path')#

		self.evclass = config.getint('General', 'Event Class')
		self.roiradius = config.getfloat('General', 'ROI Radius')
		self.emin = config.getfloat('General', 'Minimum Energy')
		self.emax = config.getfloat('General', 'Maximum Energy')
		self.zmax = config.getfloat('General', 'Maximum Zenith Angle')

		self.roicut = config.get('General', 'roicut')
		self.filter = config.get('General', 'filter')

		self.dcostheta = config.getfloat('General', 'dcostheta')
		self.ltcube_binsz = config.getfloat('General', 'ltcube bin size')

		self.ccube_binsz = config.getfloat('General', 'CCUBE binsize')
		self.num_ebins = config.get('General', 'num_ebins')
		self.ccube_xpix = config.getint('General', 'CCUBE xpix')
		self.ccube_ypix = config.getint('General', 'CCUBE ypix')

		self.irfs = config.get('General', 'Irfs')
		self.iso_base = config.get('General', 'Isotropic Base')
		self.iso_path = config.get('General', 'Path to Isotropic Base')
		self.gal_diff_base = config.get('General', 'Galactic Diffuse Base')
		self.gal_diff_path = config.get('General', 'Path to Galactic Diffuse Base')

		self.expmap_xpix = config.getint('General', 'Expmap xpix')
		self.expmap_ypix = config.getint('General', 'Expmap ypix')
		self.expmap_binsz = config.getfloat('General', 'Expmap bin size')

		self.uexpmap_radius = config.getfloat('General', 'Unbinned Expmap Radius')
		self.uexpmap_long = config.getint('General', 'Unbinned Expmap Longitude Points')
		self.uexpmap_lat = config.getint('General', 'Unbinned Expmap Latitude Points')
		self.uexpmap_num_ebins = config.getint('General', 'Unbinned Expmap Number of Energies')

		self.optimizer = config.get('General', 'Optimizer')

		print config.items('General')
		print config.items('Object')



		

		if (self.source):
			self.irfs = 'P7REP_SOURCE_V15'
			self.iso_base = 'iso_source_v05'
			self.evclass = 2
			print 'source photons'


	def select_data(self, tstart, tstop, base_name, clobber=False):
		print '....Running gtselect....'		
		filter['evclass'] = self.evclass
		filter['ra'] = self.ra
		filter['dec'] = self.dec
		#filter['rad'] = self.roiradius #gtexpmap (unbinned) doesn't like this
		filter['tmin'] = tstart
		filter['tmax'] = tstop
		filter['zmax'] = self.zmax
		filter['infile'] = '@' + self.events_list
		filter['outfile'] = base_name + '_filtered.fits'
		if(not os.path.exists(base_name + '_filtered.fits') or clobber):
			filter.run()
		else:
			print base_name + '_filtered.fits already exists'
		print '....Running gtmktime....'
		maketime['scfile'] = self.base_name+ '_SC.fits'
		maketime['roicut'] = self.roicut
		maketime['filter'] = self.filter
		maketime['evfile'] = base_name + '_filtered.fits'
		maketime['outfile'] = base_name + '_gti.fits'
		
		if(not os.path.exists(base_name + '_gti.fits') or clobber):
			maketime.run()
		else:
			print base_name + '_filtered.fits already exists'

	def livetime_cube(self, base_name, clobber=False):
		print "....Running gtltcube...."
		expCube['evfile'] = base_name + '_gti.fits'
		expCube['scfile'] = self.base_name + '_SC.fits'
		expCube['outfile'] = base_name + '_ltcube.fits'
		expCube['dcostheta'] = self.dcostheta
		expCube['binsz'] = self.ltcube_binsz
		if (not os.path.exists(base_name + '_ltcube.fits') or clobber):
			expCube.run()
		else:
			print base_name +'_ltcube.fits already exists'

	def xml_model(self, base_name, clobber=False):
		if(not os.path.exists(base_name + '_model.xml') or clobber):
			import make2FGLxml
			model = make2FGLxml.srcList('gll_psc_v08.fit', base_name + '_gti.fits',base_name + '_model.xml')
			model.makeModel(self.gal_diff_base + '.fit',self.gal_diff_base, self.iso_base + '.txt',self.iso_base)
		else:
			print base_name + '_model.xml already exists'


	def binned_ccube(self, base_name, clobber=False):
		print '....Running gtbin:CCUBE....'				
		command = 'gtbin evfile=' + base_name + '_gti.fits'\
			+ ' outfile=' + base_name + '_CCUBE.fits'\
			+ ' algorithm=CCUBE ebinalg=LOG emin=' + str(self.emin)\
			+ ' emax=' + str(self.emax)\
			+ ' enumbins=' + str(self.num_ebins)\
			+ ' nxpix=' + str(self.ccube_xpix)\
			+ ' nypix=' + str(self.ccube_ypix)\
			+ ' binsz=' + str(self.ccube_binsz)\
			+ ' axisrot=0. proj=AIT coordsys=CEL xref=' + str(self.ra)\
			+ ' yref=' + str(self.dec)\
			+ ' scfile=NONE'
		if(not os.path.exists(base_name + '_CCUBE.fits') or clobber):
			os.system(command)
		else:
			print base_name + '_CCUBE.fits already exists'

	def binned_expmap(self, base_name, clobber=False):
		print '....Running gtexpcube2....'				
		command = 'gtexpcube2 infile=' + base_name + '_ltcube.fits'\
			+ ' cmap=none outfile=' + base_name + '_bin_expmap.fits'\
			+ ' irfs=' + self.irfs\
			+ ' nxpix=' + str(self.expmap_xpix)\
			+ ' nypix=' + str(self.expmap_ypix)\
			+ ' binsz=' + str(self.expmap_binsz)\
			+ ' axisrot=0. proj=AIT coordsys=CEL xref=' + str(self.ra)\
			+ ' yref=' + str(self.dec)\
			+ ' emin=' + str(self.emin)\
			+ ' emax=' + str(self.emax)\
			+ ' enumbins=' + str(self.num_ebins)
		if(not os.path.exists(base_name + '_bin_expmap.fits') or clobber):
			os.system(command)
		else:
			print base_name + '_bin_expmap.fits already exists'

	def unbinned_expmap(self, base_name, clobber=False):
		print '....Running gtexpmap...'
		expMap['evfile'] = base_name + '_gti.fits'
		expMap['scfile'] = self.base_name + '_SC.fits'
		expMap['expcube'] = base_name + '_ltcube.fits'
		expMap['outfile'] = base_name + '_expmap.fits'
		expMap['irfs'] = self.irfs
		expMap['srcrad'] = self.uexpmap_radius
		expMap['nlong'] = self.uexpmap_long
		expMap['nlat'] = self.uexpmap_lat
		expMap['nenergies'] = self.uexpmap_num_ebins
		if(not os.path.exists(base_name + '_expmap.fits') or clobber):
			expMap.run()
		else:
			print base_name + '_expmap.fits already exists'

	
	def binned_srcMaps(self, base_name, clobber=False):
		print '....Running gtsrcmaps....'
		command = 'gtsrcmaps scfile=' + base_name + '_SC.fits'\
			+ ' expcube=' + base_name + '_ltcube.fits'\
			+ ' cmap=' + base_name + '_CCUBE.fits'\
			+ ' srcmdl=' + base_name + '_model.xml'\
			+ ' bexpmap=' + base_name + '_bin_expmap.fits'\
			+ ' outfile=' + base_name + '_srcMaps.fits'\
			+ ' irfs=' + self.irfs
		if(not os.path.exists(base_name + '_srcMaps.fits') or clobber):
			os.system(command)
		else:
			print base_name + '_srcMaps.fits already exists'

	def run_like(self, base_name):
		print '....Running gtlike....'
		if self.binned:
			obs = BinnedObs(srcMaps = base_name + '_srcMaps.fits', expCube = base_name + '_ltcube.fits', binnedExpMap = base_name + '_bin_expmap.fits', irfs = self.irfs)
			like = BinnedAnalysis(obs, base_name + '_model.xml', optimizer=self.optimizer)
		else:
			obs = UnbinnedObs(base_name + '_gti.fits', self.base_name + '_SC.fits', expMap = base_name + '_expmap.fits', expCube = base_name + '_ltcube.fits', irfs = self.irfs)
			like = UnbinnedAnalysis(obs, base_name + '_model.xml', optimizer=self.optimizer)
		
		if ('alpha' in like[self.main_source].funcs['Spectrum'].params):
			like.freeze(like.par_index(self.main_source, 'beta'))
			self.spectral_index_name = 'alpha'
		else:
			self.spectral_index_name = 'Index'
		
		like.tol = .1
		
		#like = self.fix_most(like)
		like = self.fix_nonsource(like)
	
		
		#like.freeze(like.par_index(self.main_source, 'alpha'))
		#like_return = pyLike.Minuit(like.logLike)
		like.fit(verbosity=1, covar=True)
		#print like_return.getQuality()
		like.logLike.writeXml(base_name + '_out.xml')
		#like.writeCountsSpectra()
		
		
		print like.Ts(self.main_source)
		return like

	def run_upperlimit(self, like):
		like.freeze(like.par_index(self.main_source, self.spectral_index_name))
		return IUL.calc_chi2(like, self.main_source, cl=self.ul_confidence, freeze_all=True, verbosity=1)
		#return IUL.calc_int(like, self.main_source, cl=self.ul_confidence, freeze_all=True, verbosity=1)
		
	def remove_low_Ts(self, like, Ts_min, fix=False):
		#basically a copy of method by same name in LatAnalysisScripts
		for source in like.sourceNames():
			ts = like.Ts(source)
			if (ts < Ts_min):
				print 'Removing ' + source + ' with Ts=' + str(ts)
				like.deleteSource(source)
		return like

		
	def fix_nonsource(self, like):
		saved_source = like.deleteSource(self.main_source)
		for i in range(len(like.params())):
			like.freeze(i)
		like.addSource(saved_source)
		return like
	
	def fix_most(self, like): ##fix all but source and diffuse/extragalactic 
		saved_sources = []
		saved_sources.append(like.deleteSource(self.main_source))
		saved_sources.append(like.deleteSource(self.iso_base))
		saved_sources.append(like.deleteSource(self.gal_diff_base))
		for i in range(len(like.params())):
			like.freeze(i)
		for x in saved_sources:
			like.addSource(x)
		return like

	def lightcurve(self, processes=1, binlist=[], binfile=''):
		systime = time.time()

		if (binfile==''):
			timebins = range(self.datastart, self.datastop, self.tbin)
			starttimes = timebins[:len(timebins)-1]
			stoptimes = timebins[1:]
		
		else:
			custom_bins = np.loadtxt(binfile)
			starttimes = custom_bins[:,0]
			stoptimes = custom_bins[:,1]

		curve = []
		curvepath = self.base_name + '_curve.txt'
		#curvefile = open(curvepath, 'w')
		#curvefile.write('#'+self.base_name + ' flux and flux error,ts, npred (alpha and beta fixed)\n')
		#curvefile.close()


		if(len(binlist)==0):
			binlist = range(len(starttimes))
	
		pool = mp.Pool(processes)
		args = [[[starttimes[i], stoptimes[i], i], self] for i in binlist]

		pool.map(bin_shell, args)

		
		totaltime = time.time()-systime
		print 'Total time (seconds): ' + str(totaltime)
		print 'Total time (hours): ' + str(totaltime/3600)
		print 'Seconds per bin: ' + str(totaltime/len(starttimes))




	def create_adaptive_input(self):
		if (not os.path.exists(self.adaptive_config)):
			print 'Creating Adaptive Bin Config ' + self.adaptive_config 
			config_file = open(self.adaptive_config, 'w')
			config_file.write('ft2file = ' + self.base_name + '_SC.fits\n'\
					+ 'ft1file = ' + self.base_name + '_full_gti.fits\n'\
					+ 'MacroC = TS_estimate_P7.C\n'
					+ 'rspfunc=' +  self.adaptive_rspfunc + '\n'\
					+ 'critval =' + str(self.critval) + '\n'\
					+ 'crit = ' + self.crit + '\n'\
					+ 'reverse = ' + self.reverse + '\n'\
					+ 'expfile = ' + self.base_name + '_adaptive_exp.txt\n'\
					+ 'anafile= ' + self.adaptive_out_file + '\n'\
					+ 'SRC_RA = ' + str(self.ra) + '\n'\
					+ 'SRC_DEC = ' + str(self.dec) + '\n'\
					+ 'Index = ' + str(self.spectral_index) + '\n'\
					+ 'Flux = ' + str(self.seed_flux) + '\n'\
					+ 'Elowerl = ' + str(self.emin) + '\n'\
					+ 'normeg=1\n'\
					+ 'normgal=1\n'\
					+ 'indexgal=0')
				##NOTE: THIS MAY END UP CALCULATING TOO MUCH EXPOSURE
				##GTI for whole interval, MUST BE SORTED BY DATE
		else:
			print self.adaptive_config + ' already exists.'
			
	def find_adaptive_bins(self, seed_flux = 'DEFAULT', spectral_index = 'DEFAULT'):
		#load adaptive config section
		config = ConfigParser.RawConfigParser()
		config.read(self.config_path)
		
		self.adaptive_config = config.get('Adaptive', 'Adaptive Config File Name')
		self.adaptive_rspfunc = config.get('Adaptive', 'rspfunc')
		self.critval = config.getfloat('Adaptive', 'critval')
		self.crit = config.get('Adaptive', 'crit')
		self.reverse = config.get('Adaptive', 'reverse')
		self.adaptive_out_file = config.get('Adaptive', 'Out File Name')
		self.adaptive_tstart = config.get('Adaptive', 'Start Time')
		self.adaptive_tstop = config.get('Adaptive', 'Stop Time')
		
		print config.items('Adaptive')
		#find default seed flux/spectral index values if necessary
		if(seed_flux == 'DEFAULT' or spectral_index == 'DEFAULT'):
			catalog_table = fits.open(self.catalog_file)[1].data
			names = catalog_table.field('Source_Name')
			full_catalog_name = '2FGL ' + self.main_source[5:]
			i = np.where(names == full_catalog_name)[0][0] #sets i to row of table holding our object
			self.seed_flux = catalog_table.field('Flux_Peak')[i] * .30 #seed flux is not that important of a number, and usually this is a reasonable guess
			self.spectral_index = catalog_table.field('Spectral_Index')[i]
			
		#if the seed flux/index is provided, use that (is there a more pythonic way to do this?)
		if(seed_flux != 'DEFAULT'):
			self.seed_flux = seed_flux
		if(spectral_index != 'DEFAULT'):
			self.spectral_index = spectral_index
				
		print 'Index: ' + str(self.spectral_index)
		print 'Seed_flux: ' + str(self.seed_flux)
		
		#create adaptive input file
		self.create_adaptive_input()
		
		#select data and sort
		self.select_data(self.adaptive_tstart, self.adaptive_tstop, base_name = self.base_name + '_full')
		os.system('fsort ' + self.base_name + '_full_gti.fits TIME method=\"heap\"')
		
		#find exposure
		os.system('python comp_exposure_phi.py ' + self.adaptive_config)
		
		#find bins
		os.system('python time_estimate_B.py ' + self.adaptive_config) 
				



	
	def run_bin(self, bininfolist):
		starttime, stoptime, bin_number = bininfolist
		base_name = self.base_name + str(bin_number+1)
		
		folder = './'+ self.base_name + '_bin' + str(bin_number+1)		
		try:
			os.mkdir(folder)
		except OSError:
			if (not os.path.isdir(folder)):
				raise
		base_name = folder + '/' + base_name

		print starttime, stoptime
		try:
			self.select_data(starttime, stoptime, base_name)
		except:
			if(self.debug):
				raise
			results = open(base_name + '_results.txt', 'w')
			results.write('#' + str(starttime) + '     Failure (data selection)')
			results.close()
			print '#####  FAILURE (data selection)  ####'
			return
			
		self.livetime_cube(base_name = base_name)##this is a really non-object oriented way to do it...
		self.xml_model(base_name = base_name)

		if(self.binned):
			self.binned_ccube(base_name = base_name)	#should probably create a seperate object to run this object repeatedly with everything intialized correctly
			self.binned_expmap(base_name = base_name)
			self.binned_srcMaps(base_name = base_name)
		else:
			self.unbinned_expmap(base_name = base_name)
			

		try:
			like = self.run_like(base_name = base_name)
			flux = like.flux(self.main_source)
			fluxerr = like.fluxError(self.main_source)
			
			
			index = like.model[self.main_source].funcs['Spectrum'].getParam(self.spectral_index_name).value()
			indexerr = like.model[self.main_source].funcs['Spectrum'].getParam(self.spectral_index_name).error()
			ts = like.Ts(self.main_source)
			npred = like.NpredValue(self.main_source)

			upperlimit = self.run_upperlimit(like)

			results = open(base_name + '_results.txt', 'w')
			results.write(str(starttime) + '     '\
					+ str(stoptime)+ '     '\
					+ str(Light_Curve.met_to_mjd(starttime)) + '     '\
					+ str(Light_Curve.met_to_mjd(stoptime)) + '     '\
					+ str(flux) + '     ' \
					+ str(fluxerr)	+ '     '\
					+ str(upperlimit[0]) + '     '\
					+ str(index) + '     '\
					+ str(indexerr) +  '     '\
					+ str(ts) +  '     '\
					+ str(npred))
			results.close()
		except:
			if(self.debug):
				raise
			results = open(base_name + '_results.txt', 'w')
			results.write('#' + str(starttime) + '     Failure (likelihood)')
			results.close()
			print '#####  FAILURE (likelihood) ####'
			return

		return

	def read_bins(self):#go through consecutive folders starting with bin1 and pull results 
		i = 1
		results = []
		path = './' + self.base_name + '_bin' + str(i) + '/' + self.base_name + str(i) + '_results.txt'
		while(os.path.exists(path)):
			resultbin = open(path, 'r')
			results.append(resultbin.readline())
			resultbin.close()
			i+=1
			path = './' + self.base_name + '_bin' + str(i) + '/' + self.base_name + str(i) + '_results.txt'
		outfile = open('results.txt', 'w')
		outfile.write('#' + self.base_name + ' metstart, flux, flux error, upperlimit, spectral index, index error, ts, npred (all fixed but norm and spectral index)\n')
		for j in range(len(results)):
			outfile.write(results[j] + '\n')
		outfile.close()
	
	def read_bin_list(self, binlist):
		results = []
		for i in binlist:
			path = './' + self.base_name + '_bin' + str(i) + '/' + self.base_name + str(i) + '_results.txt'
			if(os.path.exists(path)):
				resultbin = open(path, 'r')
				results.append(resultbin.readline())
				resultbin.close()
			else:
				results.append('#bin ' + str(i) + ' not found')
		
		outfile = open('results.txt', 'w')
		outfile.write('#' + self.base_name + ' metstart, metstop, MJDstart, MJDstop, flux, flux error, upperlimit, spectral index, index error, ts, npred (all fixed but norm and index)\n')
		for j in range(len(results)):
			outfile.write(results[j] + '\n')
		outfile.close()		

	@staticmethod
	def text_to_fits(filepath = 'results.txt'):  ##Note: this is a temporary workaround
		data = np.loadtxt(filepath)
		columns = []
		#mjd = Light_Curve.met_to_mjd(data[:,0])
		#data = np.c_[data[:,0], mjd, data[:,1:]]
		names = ['METstart', 'METstop', 'MJDstart', 'MJDstop', 'Flux', 'Flux_err', 'UL', 'Index', 'Index_err', 'TS', 'Npred']
		units = ['seconds', 'seconds', 'days', 'days', 'PH/cm^2/s', 'PH/cm^2/s', 'PH/cm^2/s', '', '', '', 'Counts']
		for i in range(len(data[0])):
			col = fits.Column(name = names[i], format='E', array=data[:,i], unit = units[i])
			columns.append(col)
		
		bintable = fits.new_table(columns)
		bintable.writeto(filepath[:-3]+'fits') #filepath[:-3] is the name without the .txt 

	@staticmethod
	def met_to_mjd(time_met):
		return time_met/86400.0 + 51910 + 0.000742870370

	def setup_links(self):	
		if(not os.path.exists(self.iso_base + '.txt')):
			os.system('ln -s ' + self.iso_path + ' ./')
		#if(not os.path.exists('iso_source_v05.txt')):
		#	os.system('ln -s $FERMI_DIR/refdata/fermi/galdiffuse/iso_source_v05.txt ./')
		if(not os.path.exists(self.gal_diff_base + '.fit')):
			os.system('ln -s ' + self.gal_diff_path + ' ./')
		if(not os.path.exists(self.catalog_file)):
			os.system('ln -s ' + self.fgl_path + ' ./')
		if(not os.path.exists('make2FGLxml.py')):
			os.system('ln -s ' + self.makexml_path + ' ./')


		photon_files = glob.glob(self.data_path + '*_PH??.fits')
		for path in photon_files:
			if(not os.path.exists(path[len(self.data_path):])): 
				os.system('ln -s ' + path + ' ./')

		if(not os.path.exists('events.txt')):
			events = open('events.txt', 'w')
			for path in photon_files:
				events.write(path[len(self.data_path):] + '\n') 
			events.close()



		scfile = glob.glob(self.data_path + '*_SC00.fits')
		if(not os.path.exists(self.base_name + '_SC.fits')):
			os.system('ln -s ' + scfile[0] + ' ./' + self.base_name + '_SC.fits')


		
			