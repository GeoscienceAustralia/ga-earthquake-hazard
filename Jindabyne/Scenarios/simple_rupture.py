import shapefile
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth

def checkfloat(floatstr):
    try:
        return float(floatstr)
    except:
        from numpy import nan
        return nan
        

def checkint(intstr):
    try:
        return int(intstr)
    except:
        from numpy import nan
        return nan


def get_field_index(sf, field):
    fields = sf.fields[1:]
    findex = -1
    for i, f in enumerate(fields):
        if f[0] == field:
            findex = i      
    return findex


def get_field_data(sf, field, datatype):
    '''
    datatype = float, str
    '''      
    # get index
    index = get_field_index(sf, field)
    
    # get records
    recs = sf.records()
    
    # now loop thru recs and get data
    data = []
    for rec in recs:
        if datatype == 'str':
            data.append(rec[index])
        elif datatype == 'float':
            data.append(checkfloat(rec[index]))
        elif datatype == 'int':
            data.append(checkint(rec[index]))
            
    return np.array(data)


def get_line_parallels(pts, rngkm):
    '''
    pts are an N x [lon, lat] matrix, i.e.:
                   [[lon1, lat1],
                    [lon2, lat2]]
    '''    
    # set outputs
    posazpts = []
    negazpts = []
    
    for j, pt in enumerate(pts):
        # if 1st point
        if j == 0:
            rngm, az, baz = gps2dist_azimuth(pts[j][1], pts[j][0], \
                                            pts[j+1][1], pts[j+1][0])
            
        # if last point
        elif j == len(pts)-1:
            rngm, az, baz = gps2dist_azimuth(pts[j-1][1], pts[j-1][0], \
                                            pts[j][1], pts[j][0])
                                           
        # use points either side (assumes evenly spaced)
        else:
            rngm, az, baz = gps2dist_azimuth(pts[j-1][1], pts[j-1][0], \
                                            pts[j+1][1], pts[j+1][0])
           
        # get azimuth for new points
        azpos = az + 90.
        azneg = az - 90.
        # get points
        posazpts.append(reckon(pts[j][1], pts[j][0], rngkm, azpos))
        negazpts.append(reckon(pts[j][1], pts[j][0], rngkm, azneg))
    
    '''    
    # for testing only
    x=[]
    y=[]
    xx=[]
    yy=[]
    xxx=[]
    yyy=[]
    for j, pt in enumerate(pts):
        x.append(pt[0])
        y.append(pt[1])
        xx.append(posazpts[j][0])
        yy.append(posazpts[j][1])
        xxx.append(negazpts[j][0])
        yyy.append(negazpts[j][1])
        
    plt.plot(x,y,'b-')
    plt.plot(xx,yy,'r-')
    plt.plot(xxx,yyy,'g-')
    plt.show()
    '''
    return posazpts, negazpts


def reckon(lat1d, lon1d, rngkm, brngd):
    # calculate lat lon from range (km) and bearing (degrees)
    from math import radians, asin, sin, cos, atan2, degrees   
    
    R = 6378.1 #Radius of the Earth
    brng = radians(brngd)
    
    lat1r = radians(lat1d) #Current lat point converted to radians
    lon1r = radians(lon1d) #Current long point converted to radians
    
    lat2r = asin(sin(lat1r)*cos(rngkm/R) + cos(lat1r)*sin(rngkm/R)*cos(brng))

    lon2r = lon1r + atan2(sin(brng)*sin(rngkm/R)*cos(lat1r), \
            cos(rngkm/R)-sin(lat1r)*sin(lat2r))

    lat2d = degrees(lat2r)
    lon2d = degrees(lon2r)
    
    return [lon2d, lat2d]

def mag2rupwid_WC94(mw, ftype):
    if ftype == 'all':
        a = -1.01
        b = 0.32
    elif ftype == 'ss':
        a = -0.76
        b = 0.27
    elif ftype == 'rs':
        a = -1.61
        b = 0.41
    elif ftype == 'ns':
        a = -1.14
        b = 0.35
    return 10**(a + b * mw) # in km




shpfile = "Morwell/Scenarios/Rosedale_monocline/Rosedale_monocline.shp"
print('Reading source shapefile...')
sf = shapefile.Reader(shpfile)

ids = get_field_data(sf, 'ID', 'str')
mw = get_field_data(sf, 'MW', 'float')
dip = get_field_data(sf, 'DIP', 'float')

shapes = sf.shapes()

for i, shape in enumerate(shapes):
	
	#rupwid = mag2wid_L10(mw[i], 'scr')
	rupwid = mag2rupwid_WC94(mw[i], 'rs')
	print("Rupture width %s" % rupwid)
    #surface projection of top of rupture
	rngkm = rupwid * np.cos(np.radians(dip[i]))
	print("Surface rupture width:  %s" % rngkm)
	rupdep = rupwid * np.sin(np.radians(dip[i]))
	print("Rupture depth:  %s" % rupdep)
	
	
	posazpts, negazpts = get_line_parallels(shape.points, rngkm)
	
	ftxt = '\t'.join((str('%0.3f' % shape.points[0][0]),
                      str('%0.3f' % shape.points[0][1]), '0.0')) + '\n'
	ftxt += '\t'.join((str('%0.3f' % shape.points[1][0]),
                       str('%0.3f' % shape.points[1][1]), '0.0')) + '\n\n'
	
	ftxt += '\t'.join((str('%0.3f' % posazpts[1][0]),
                       str('%0.3f' % posazpts[1][1]), 
                       str('%0.1f' % rupdep))) + '\n'
	ftxt += '\t'.join((str('%0.3f' % posazpts[1][0]),
                       str('%0.3f' % posazpts[0][1]), 
                       str('%0.1f' % rupdep))) + '\n\n'
	
	ftxt += '\t'.join((str('%0.3f' % negazpts[1][0]),
                       str('%0.3f' % negazpts[1][1]), 
                       str('%0.1f' % rupdep))) + '\n'   
	ftxt += '\t'.join((str('%0.3f' % negazpts[0][0]),
                       str('%0.3f' % negazpts[0][1]), 
                       str('%0.1f' % rupdep))) + '\n\n'
	
	ftxt += '\t'.join((str('%0.3f' % shape.points[0][0]),
                       str('%0.3f' % shape.points[0][1]), '0.0'))                 
	print(ftxt)	
	print('%s_fault.txt' % str(ids[i]))
	f = open('%s_fault.txt' % str(ids[i]), 'w')
	print(ftxt)
	f.write(ftxt)

	#print(ftxt, f)
	f.close()

