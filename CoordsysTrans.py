import math
from math import sin,cos,atan2,sqrt,fabs
from math import pi as PI


# define ellipsoid
a = 6378245.0  # 长半轴
f = 1/298.3
b = a*(1-f)
#ee = 0.00669342162296594323  扁率
ee = 1-(b*b)/(a*a)


# check if the point in CHINA
def out_of_china(lng,lat):
    if lng < 72.004 or lng > 137.8347:
        return True
    if lat < 0.8293 or lat > 55.8271:
        return True
    return False


def transformLat(lng,lat):
    ret = -100.0 + 2*lng +3.0*lat + 0.2*lat*lat + 0.1*lng*lat + 0.2*sqrt(fabs(lng))
    ret = ret + (20.0* sin(6.0*lng*PI) + 20.0*sin(2.0*lng*PI))*2.0/3.0
    ret = ret + (20.0* sin(lat*PI) + 40.0*sin(lat/3.0*PI))*2.0/3.0
    ret = ret + (160.0*sin(lat/12.0*PI) + 320.0*sin(lat*PI/30.0))*2.0/3.0
    return ret


def transformLng(lng,lat):
    ret = 300.0 + lng + 2.0*lat + 0.1*lng*lng + 0.1*sqrt(fabs(lng))
    ret = ret + (20.0*sin(6.0*lng*PI)+ 20.0*sin(2.0*lng*PI))*2.0/3.0
    ret = ret + (20.0*sin(lng*PI) + 40.0*sin(lng/3.0*PI))*2.0/3.0
    ret = ret + (150.0*sin(lng/12.0*PI) + 300.0*sin(lng*PI/30.0))*2.0/3.0
    return ret


def wgs2gcj(wgslng,wgslat):
    if (out_of_china(wgslng,wgslat)):
        return [wgslng,wgslat]
    dlat = transformLat(wgslng-105.0,wgslat-35.0)
    dlng = transformLng(wgslng-105.0,wgslat-35.0)
    radlat = wgslat/180.0*PI
    magic = math.sin(radlat)
    magic = 1 - ee*magic*magic
    sqrtMagic = sqrt(magic)
    dlat = (dlat*180.0)/((a*(1-ee))/(magic*sqrtMagic)*PI)
    dlng = (dlng*180.0)/(a/sqrtMagic*cos(radlat)*PI)
    gcjLat = wgslat + dlat
    gcjLng = wgslng + dlng
    return (gcjLng,gcjLat)


def gcj2wgs(gcjLng,gcjLat):
    g0 = (gcjLng,gcjLat)
    w0 = g0
    g1 = wgs2gcj(w0[0],w0[1])
    # w1 = w0-(g1-g0)
    w1 = tuple([x[0]-(x[1]-x[2]) for x in zip(w0,g1,g0)])
    # delta = w1-w0
    delta = tuple([x[0]-x[1] for x in zip(w1,w0)])
    while (abs(delta[0])>= 1e-6 or abs(delta[1]) >= 1e-6):
        w0 = w1
        g1 = wgs2gcj(w0[0],w0[1])
        w1 = tuple([x[0]-(x[1]-x[2]) for x in zip(w0,g1,g0)])
        delta = tuple([x[0]-x[1] for x in zip(w1,w0)])
    return w1


def gcj2bd(gcjlng,gcjlat):
    z = sqrt(gcjlng*gcjlng + gcjlat*gcjlat) + 0.00002*sin(gcjlat*PI*3000.0/180.0)
    theta = atan2(gcjlat,gcjlng) + 0.000003*cos(gcjlng*PI*3000.0/180.0)
    bdlng = z*cos(theta) + 0.0065
    bdlat = z*sin(theta) + 0.006
    return (bdlng,bdlat)


def bd2gcj(bdLng,bdLat):
    x = bdLng - 0.0065
    y = bdLat - 0.006
    z = sqrt(x*x + y*y) - 0.00002*sin(y*PI*3000.0/180)
    theta = atan2(y,x) - 0.000003*cos(x*PI*3000.0/180)
    gcjlng = z*cos(theta)
    gcjlat = z*sin(theta)
    return (gcjlng,gcjlat)


def wgs2bd(wgslng,wgslat):
    gcj = wgs2gcj(wgslng,wgslat)
    bd = gcj2bd(gcj[0],gcj[1])
    return bd

def bd2wgs(bdlng,bdlat):
    gcj = bd2gcj(bdlng,bdlat)
    wgs = gcj2wgs(gcj[0],gcj[1])
    return wgs


if __name__ == '__main__':
    # wgs2gcj
    # coord = (112, 40)
    # trans = WGS2GCJ()
    print(wgs2gcj(112, 40))
    print(gcj2wgs(112.00678230985764, 40.00112245823686))

    # gcj2wgs