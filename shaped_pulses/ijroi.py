# Copyright: Luis Pedro Coelho <luis@luispedro.org>, 2012
#            Tim D. Smith <git@tim-smith.us>, 2015
# License: MIT

import numpy as np
from matplotlib import pyplot as plt

from skimage.draw import (line, polygon, circle,
                          circle_perimeter,
                          ellipse, ellipse_perimeter,
                          bezier_curve)

def read_roi(fileobj):
    '''
    points = read_roi(fileobj)

    Read ImageJ's ROI format. Points are returned in a nx2 array. Each row
    is in [row, column] -- that is, (y,x) -- order.
    '''
    # This is based on:
    # http://rsbweb.nih.gov/ij/developer/source/ij/io/RoiDecoder.java.html
    # http://rsbweb.nih.gov/ij/developer/source/ij/io/RoiEncoder.java.html

    SPLINE_FIT = 1
    DOUBLE_HEADED = 2
    OUTLINE = 4
    OVERLAY_LABELS = 8
    OVERLAY_NAMES = 16
    OVERLAY_BACKGROUNDS = 32
    OVERLAY_BOLD = 64
    SUB_PIXEL_RESOLUTION = 128
    DRAW_OFFSET = 256

    class RoiType:
        POLYGON = 0
        RECT = 1
        OVAL = 2
        LINE = 3
        FREELINE = 4
        POLYLINE = 5
        NOROI = 6
        FREEHAND = 7
        TRACED = 8
        ANGLE = 9
        POINT = 10
    
    RoiTypeStrings = {RoiType.POLYGON:'polygon',
	                            RoiType.RECT:'rect',
				    RoiType.OVAL:'oval',
				    RoiType.LINE:'line',
				    RoiType.FREELINE:'freeline',
				    RoiType.POLYLINE:'polyline',
				    RoiType.NOROI:'noroi',
				    RoiType.FREEHAND:'freehand',
				    RoiType.TRACED:'traced',
				    RoiType.ANGLE:'angle',
				    RoiType.POINT:'point'
				    }
    def get8():
        s = fileobj.read(1)
        if not s:
            raise IOError('readroi: Unexpected EOF')
        return ord(s)

    def get16():
        b0 = get8()
        b1 = get8()
        return (b0 << 8) | b1

    def get32():
        s0 = get16()
        s1 = get16()
        return (s0 << 16) | s1

    def getfloat():
        v = np.int32(get32())
        return v.view(np.float32)

    magic = fileobj.read(4)
    if magic != b'Iout':
        raise ValueError('Magic number not found')
    version = get16()

    # It seems that the roi type field occupies 2 Bytes, but only one is used
    roi_type = get8()
    # Discard second Byte:
    get8()
    
#    print()
#    print( "roi_type",roi_type,RoiTypeStrings[roi_type])
#    print()

    if roi_type not in [RoiType.FREEHAND, RoiType.POLYGON, RoiType.RECT, RoiType.POINT, RoiType.OVAL]:
        raise NotImplementedError('roireader: ROI type %s not supported' % roi_type)

    top = get16()
    left = get16()
    bottom = get16()
    right = get16()
    n_coordinates = get16()
#    print("bounding box: ", top, left, bottom, right )
#    print("n_coords", n_coordinates )
    
    x1 = getfloat()
    y1 = getfloat()
    x2 = getfloat()
    y2 = getfloat()
#    print ('x1,y1,x2,y2: ', x1,y1,x2,y2 )
    
    stroke_width = get16()
    
#    print( "stroke_width: ", stroke_width )
    
    shape_roi_size = get32(); #print( "shape_roi_size: ", shape_roi_size )
    
    stroke_color = get32(); #print( "stroke_color: ", stroke_color )
    
    fill_color = get32(); #print( "fill_color: ", fill_color )
    subtype = get16(); #print( 'subtype: ', subtype )
    if subtype != 0:
        raise NotImplementedError('roireader: ROI subtype %s not supported (!= 0)' % subtype)
    options = get16(); #print("options: ", options )
    arrow_style = get8()
    arrow_head_size = get8()
    rect_arc_size = get16()
    position = get32()
    header2offset = get32()
    
    if roi_type == RoiType.OVAL:
	    height = bottom-top
	    width  = right-left
	    return RoiTypeStrings[roi_type], np.array( [[left+width/2.,top+height/2.],[width/2.,height/2.]], dtype=np.float32)

    if roi_type == RoiType.RECT:
        if options & SUB_PIXEL_RESOLUTION:
            return RoiTypeStrings[roi_type], np.array(
                [[y1, x1], [y1, x1+x2], [y1+y2, x1+x2], [y1+y2, x1]],
                dtype=np.float32)
        else:
            return RoiTypeStrings[roi_type], np.array(
                [[top, left], [top, right], [bottom, right], [bottom, left]],
                dtype=np.int16)

    if options & SUB_PIXEL_RESOLUTION:
        getc = getfloat
        points = np.empty((n_coordinates, 2), dtype=np.float32)
        fileobj.seek(4*n_coordinates, 1)
    else:
        getc = get16
        points = np.empty((n_coordinates, 2), dtype=np.int16)

    points[:, 1] = [getc() for i in range(n_coordinates)]
    points[:, 0] = [getc() for i in range(n_coordinates)]

    if options & SUB_PIXEL_RESOLUTION == 0:
        points[:, 1] += left
        points[:, 0] += top

    return RoiTypeStrings[roi_type], points


def read_roi_zip(fname):
    import zipfile
    with zipfile.ZipFile(fname) as zf:
        return [(n, read_roi(zf.open(n))) for n in zf.namelist()]
	
	
	
if __name__ == "__main__":
    fn = r"E:\documents\work\projects\Newcastle\ipython\read_roi-master\RoiSet_square_circle_oval_rectangle_polygon.zip"
#	fn = r"E:\documents\work\projects\Newcastle\ipython\read_roi-master\RoiSet_circles.zip"
#	fn = r"E:\documents\work\projects\Newcastle\ipython\read_roi-master\RoiSet_slice0.zip"

    rois = read_roi_zip(fn)
	
    for roi in rois:
        print(roi[0])
        print(roi[1][0])
        print(roi[1][1])
        
    numRows = 400
    numCols = 400
    img = np.zeros((numRows,numCols))
    
    
    
#    plt.imshow(img.reshape(numRows,numCols))
#    plt.show()
   
           
    #image_slice = t2m_i
    #img = np.zeros((3,192,192))
    #img[0]=t2m_i
    fig = plt.figure(figsize=(6,6))
    #plt.imshow(image_slice,zorder=1)
    

    
    for roi in rois:
        (roiName,(roiShape,coords))=roi
        print("roiName: ", roiName,"roiShape: ",roiShape)
        
        if roiShape == 'polygon':        
            coords=coords.transpose()
            xxx = coords[0]
            yyy = coords[1]
            (rr,cc) = polygon(xxx,yyy)
            img[rr, cc] = 100
            
        elif roiShape == 'rect':        
            coords=coords.transpose()
            xxx = coords[0]
            yyy = coords[1]
            (rr,cc) = polygon(xxx,yyy)
            img[rr, cc] = 200
                        
        elif roiShape == 'oval':
            xc = coords[0][0]
            yc = coords[0][1]
            xr = coords[1][0]
            yr = coords[1][1]
            print(coords)
            print(xc,yc)
            (rr,cc) = ellipse(yc,xc,yr,xr)
            img[rr, cc] = 300
            
    plt.imshow(img)  
	
	
