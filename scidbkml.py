'''
ScidbKML
~~~~~~~~

Quick library to convert Scidb arrays to Kml
    - requires ==> python-gdal, matplotlib, numpy, simplekml
    - requires ==> scidbload, scidbproject 
        {www.github.com/nicksteiner/scidbload, www.github.com/nicksteiner/scidbproject}
'''


import os
import sys
import shutil
import zipfile
import datetime
import subprocess as sp

#import yaml
import ogr
from osgeo import gdal, osr
import numpy as np
import ipdb
from matplotlib import pyplot, mpl, cm
import matplotlib.pyplot as plt
import simplekml

from scidbload import scidbload
from scidbproject import AlaskaAlbers1km


'''
Plotting Arguments
~~~~~~~~~~~~~~~~~~
'''

SDB_FMT = { 'default': {
                'grid_name':'alaska_albers', 
                'clim': (225, 325), 
                'cmap': 'spectral', 
                'scale': 0.01,
                'null': -9999.0, 
                'label':''
                },
            'modnrt_lst': {
                'grid_name':'alaska_albers', 
                'clim': (225, 325), 
                'cmap': 'spectral', 
                'scale': 1,
                'null': -9999.0, 
                'label':''
                },
            'testkml_x': {
                'grid_name':'alaska_albers', 
                'clim': (0, 255), 
                'cmap': 'spectral', 
                'scale': 1,
                'null': -9999.0, 
                'label':''
                },
            'test_regrid_ch10v_asc': {
                'grid_name':'alaska_albers', 
                'clim': (180, 270),
                'cmap': 'spectral', 
                'scale': .01,
                'null': -9999.0, 
                'label':''
                },
            'test_loading_ch36v_asc': {
                'grid_name':'alaska_albers', 
                'clim': (180, 310),
                'cmap': 'spectral', 
                'scale': .01,
                'null': -9999.0, 
                'label':''
                }
        }


'''
OBJFactories
~~~~~~~~~~~~
'''

def write_kml(kml_path, sdb_array_obj=None, attribute_list=None, query=False, plt_arg=None):
    if not sdb_array_obj:
        if query:
            sdb_array_obj = scidbload.Query(sdb_array)
        else:
            sdb_array_obj = scidbload.ScidbArray(sdb_array)
    
    if not attribute_list:
        attribute_list = sdb_array_obj.attributes
    
    if not all([att in sdb_array_obj.attsLoaded for att in attribute_list]):
        sdb_array_obj.load(attribute_list)

    gtiff_list = [get_geotiff(sdb_array_obj, att, plt_arg=plt_arg) for att in attribute_list]
    kml = EasyKml(kml_path, gtiff_list, sdb_array_obj.attributes)
    kml.write_kmz()

def get_geotiff(scidb_array, attribute, plt_arg=None):
    npy_array = scidb_array.sa[attribute]
    if not plt_arg:
        plt_arg = '{}_{}'.format(scidb_array.description.Name, attribute)
    #(self, array_name, attribute_name, proj_obj, npy_array, file_path='.', pltargs=None):
    geotiff = Geotiff(scidb_array.description.Name,
            attribute,  AlaskaAlbers1km, npy_array,
            pltargs=SDB_FMT[plt_arg])
    return geotiff

'''
OO-Stuff
~~~~~~~~
'''

class EasyKml(object):
    def __init__(self, kml_path, gtiff_list, trim=None):
        # needs geotiff and colorbar for each folder
        # create paths
        try:
            os.mkdir(kml_path)
        except: pass
        self.kml_path = kml_path
        self.img_path = os.path.join(kml_path, 'img')
        try:
            os.mkdir(self.img_path)
        except: pass
        kml_root, kml_fold = os.path.split(kml_path)
        self.kml_fold = kml_path
        self.kml_file = os.path.join(kml_path, '{0}.kml'.format(kml_fold))
        # initialize gtiffs
        self.gtiff_list = gtiff_list
        self.kml = simplekml.Kml()
        self.trim = trim

    def write_kmz(self):
        # write tiffs into the kml format
        for gtiff in self.gtiff_list:
            gtiff.file_path = self.img_path
            gtiff.write_geotiff()
            gtiff.reproject_latlon()
            gtiff.write_colorbar()
            self.write_layer(gtiff)

        self.kml.save(self.kml_file)
        # compress and add to kmz
        shutil.make_archive(self.kml_fold, 'zip', self.kml_path)
        shutil.copy(self.kml_fold + '.zip', self.kml_fold  + '.kmz')
        os.remove(self.kml_fold + '.zip')

    def write_layer(self, gtiff):
        # make folder
        fol = self.kml.newfolder(name=gtiff.name)
        # write tiff
        #gtiff = gtiff.reproject_latlon()
        # add the groundcover
        grnd = fol.newgroundoverlay(name=gtiff.name)
        grnd.icon.href = os.path.join('img', gtiff.file_name)
        grnd.gxlatlonquad.coords = gtiff.bbox
        # add the colorbar
        key_path = os.path.join(self.kml_path, gtiff.key_file)
        #shutil.copy(key_path, os.path.join(self.img_path, key_file))
        screen = fol.newscreenoverlay(name = 'key_{0}'.format(gtiff.name))
        screen.icon.href = os.path.join('img', os.path.basename(gtiff.key_file))
        screen.overlayxy = simplekml.OverlayXY(x=0,y=0,xunits=simplekml.Units.fraction,
                                       yunits=simplekml.Units.fraction)
        screen.screenxy = simplekml.ScreenXY(x=5,y=25,xunits=simplekml.Units.pixel,
                                     yunits=simplekml.Units.pixel)
        screen.size.x = 500
        screen.size.y = 100
        screen.size.xunits = simplekml.Units.pixel
        screen.size.yunits = simplekml.Units.pixel

class Geotiff(object):

    driver = gdal.GetDriverByName('GTiff')

    def __init__(self, array_name, attribute_name, proj_obj, npy_array, file_path='.', pltargs=None):
        self.proj_obj = proj_obj()
        self.srs = osr.SpatialReference()
        self.srs.ImportFromProj4(self.proj_obj.proj_str)
        self.wkt = self.srs.ExportToWkt()
        self.geotransform = self.proj_obj.geotransform
        self.geo_srs = self.srs.CloneGeogCS()

        self.trans = osr.CoordinateTransformation(self.srs, self.geo_srs)
        self.itrans = osr.CoordinateTransformation(self.geo_srs, self.srs)

        self.file_name = '{}_{}.tif'.format(array_name, attribute_name)

        self.npy_array = npy_array

        self.file_path = file_path

        self.name = '{}_{}'.format(array_name, attribute_name)
        
        self.array_name = array_name
        self.attribute_name = attribute_name

        self.cmap = cm.ScalarMappable(cmap=pltargs['cmap'])
        self.cmap.set_clim(pltargs['clim'])
        if not pltargs:
            self.pltargs = SDB_FMT['default']
        else:
            self.pltargs = pltargs
        self.carray = npy_array
        self.nodata = pltargs['null']
        
    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, value):
        self._file_name = value
    

    @property
    def file_path(self):
        return self._file_path
    @file_path.setter
    def file_path(self, file_path):
        self._file_path = os.path.realpath(file_path)
        self._geotiff = self.init_geotiff()


    @property
    def geotiff(self):
        if not hasattr(self, '_geotiff'):
            self._geotiff = self.init_geotiff()
        return self._geotiff

    @geotiff.setter
    def geotiff(self, value):
        self._geotiff = value

    def init_geotiff(self):
        # write geotiff
        target_path = os.path.join(self.file_path, self.file_name)
        dst_ds = self.driver.Create(target_path,self.npy_array.shape[1],self.npy_array.shape[0], 3, gdal.GDT_Byte)
        dst_ds.SetGeoTransform(self.geotransform)
        dst_ds.SetProjection(self.wkt)
        return dst_ds

    @property
    def carray(self):
        return self._carray
    @carray.setter
    def carray(self, npy_array):
        self._carray = 255 * self.cmap.to_rgba(npy_array * self.pltargs['scale'])

    @property
    def nodata(self):
        return self._nodata
    @nodata.setter
    def nodata(self, null):
        self._nodata =  255 * np.array(self.cmap.to_rgba(null))

    def write_geotiff(self):
        for band in range(3):
            bnd_ = self.geotiff.GetRasterBand(band + 1)
            dat_8bit = self.carray[:,:,band]
            if self.pltargs['null'] == 'null':
                dat_8bit[np.isnull(dat_8bit)] = 0
            else:
                dat_8bit[dat_8bit == self.nodata[band]] = 0
            bnd_.WriteArray(dat_8bit.astype('uint8'))
            bnd_.SetNoDataValue(0)
        self.geotiff.FlushCache()
        self.write_colorbar()

    def write_colorbar(self):
        #plot colorbar
        fig = plt.figure(figsize=(8,1.5))
        ax1 = plt.subplot(111)
        vmin, vmax = self.pltargs['clim']
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(ax1, norm=norm, cmap=self.pltargs['cmap'], orientation='horizontal')
        cb1.set_label(self.pltargs['label'])
        # PATHS
        file_name = 'key_{0}_{1}.png'.format(self.array_name, self.attribute_name)
        self.key_file = os.path.join(self.file_path, file_name)
        plt.tight_layout()
        plt.savefig(self.key_file)
        plt.close()

    def reproject_latlon(self):
        
        target_file = self.file_name.rstrip('.tif') + '_ll.tif'
        target_path = os.path.join(self.file_path, target_file)
        try:
            os.remove(target_path)
            print('target removed')
        except: pass
        
        #always writing at ~1.5km

        cmd_ = ['gdalwarp',
                '-t_srs', 'EPSG:4326',
                os.path.join(self.file_path, self.file_name), target_path,
                 '-tr', str(.015), str(.015)]
        #NOTE: WILL ONLY FIX ALASKA CASE !!!!!
        trim = self.check_trim()
        cmd_.extend(trim)
        print(' '.join(cmd_))
        sp.Popen(cmd_).communicate()
        #self.latlon_path = target_path
        # set geotiff to new file
        self.geotiff = gdal.Open(target_path)
        self.file_name = target_file
        self.reset_projection()
        
    def reset_projection(self):
        self.geotransform = self.geotiff.GetGeoTransform()
        self.srs = osr.SpatialReference(self.geotiff.GetProjectionRef())
        self.geo_srs = self.srs.CloneGeogCS()
        self.trans = osr.CoordinateTransformation(self.srs, self.geo_srs)
        self.itrans = osr.CoordinateTransformation(self.geo_srs, self.srs)

    def check_trim(self):
        bbox = np.array(self.bbox)
        lats, lons = bbox[:,1], bbox[:,0]
        signs = sum(np.sign(lons))
        # defaults to western hemisphere if straddling
        if abs(signs) != 4:
            if signs <= 0:
                trim = [-180, self.proj_obj.latitude.min(), max(lons[lons<0]), self.proj_obj.latitude.max()]
            if signs > 0:
                trim = [min(lons[lons>0]), self.proj_obj.latitude.min(), 180, self.proj_obj.latitude.max()]
            trim = ['-te'] + [str(tr) for tr in trim]
        else:
            trim = []
        return trim


    @property
    def bbox(self):
        return self.gen_bbox()

    def gen_bbox(self):
        bbox_cells = (
                (0, self.geotiff.RasterYSize),
                (self.geotiff.RasterXSize, self.geotiff.RasterYSize),
                (self.geotiff.RasterXSize, 0),
                (0., 0.)
                    )
        geo_pts = []
        for x, y in bbox_cells:
            x2 = self.geotransform[0] + self.geotransform[1] * x + self.geotransform[2] * y
            y2 = self.geotransform[3] + self.geotransform[4] * x + self.geotransform[5] * y
            geo_pt = self.trans.TransformPoint(x2, y2)[:2]
            geo_pts.append(geo_pt)
            print x, y, '->', x2, y2, '->', geo_pt
        return geo_pts