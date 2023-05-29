import scipy.io as sio
import scipy.ndimage as ndimage

class Initializer:
    
    
    def __init__(self):

        
        # BOX
        Lx = 0.27
        Ly = 0.27
                
        self.boxSize    = [Lx, Ly, 1] # in d0 - ion inertia length
        self.boxSizePxl = [460, 460, 1]

        self.dampingBoundaryWidth = [[15,15], [15,15]]
        
        self.dx = self.boxSize[0]/self.boxSizePxl[0]
        self.dy = self.boxSize[1]/self.boxSizePxl[1]
        
        mat_contents2 = sio.loadmat('./example/shot_020firstEBT_sym.mat')['sym']
        [lx,ly] = mat_contents2.shape
        print(lx,ly)
        lx2 = int(0.5*lx)-20
        ly2 = int(0.5*ly)
        mat_contents2 = mat_contents2[lx2-230:lx2+230,ly2-230:ly2+230]

        self.fluence = mat_contents2
        self.fluenceSmoothed = ndimage.gaussian_filter(self.fluence, sigma=(90,90), order=0)

        self.mpiCores  = [1,2,1]
        
        # time
        self.ts = 1e-7
        self.maxtsnum = 100001
        self.outputStride = 500
        
        # output
        self.outputDir = "output_las/"
        self.fileTemplate = "proto_"
        

    
    #   spatial: left lower corner is (0,0,0)
    #   spatial: total box length in normalized units
    def getXright(self):
        return self.boxSize[0]
    
    def getYright(self):
        return self.boxSize[1]
    
    def getZright(self):
        return self.boxSize[2]
    
    # if number of pixels is less than 2, there is no direvative
    #   total box length in pixels Lx
    def getXresolution(self):
        return self.boxSizePxl[0]
    
    #   total box length in pixels Ly
    def getYresolution(self):
        return self.boxSizePxl[1]
    
    #   total box length in pixels Lz
    def getZresolution(self):
        return self.boxSizePxl[2]
    
    def getXmpiDomainNum(self):
        return self.mpiCores[0]
    
    def getYmpiDomainNum(self):
        return self.mpiCores[1]
    
    def getZmpiDomainNum(self):
        return self.mpiCores[2]
    

    def getDampingBoundaryWidthXleft(self):
        return self.dampingBoundaryWidth[0][0]

    def getDampingBoundaryWidthXright(self):
        return self.dampingBoundaryWidth[0][1]

    def getDampingBoundaryWidthYleft(self):
        return self.dampingBoundaryWidth[1][0]

    def getDampingBoundaryWidthYright(self):
        return self.dampingBoundaryWidth[1][1]


    #   time
    def getTimestep(self):
        return self.ts
    
    def getMaxTimestepsNum(self):
        return self.maxtsnum

    
    #   output:
    #   output: folder must be created manually
    def getOutputDir(self):
        return self.outputDir
    
    def getOutputFilenameTemplate(self):
        return self.fileTemplate
    
    def getOutputTimestep(self):
        return self.outputStride
    
    
    def getProtonFluence(self, x, y, z):
        if (x >= self.boxSizePxl[0] ):
            x = self.boxSizePxl[0]-1
        if (y >= self.boxSizePxl[1] ):
            y = self.boxSizePxl[1]-1
        return self.fluence[x,y]

    def getBackgroundProtonFluence(self, x, y, z):
        if (x >= self.boxSizePxl[0] ):
            x = self.boxSizePxl[0]-1
        if (y >= self.boxSizePxl[1] ):
            y = self.boxSizePxl[1]-1
        return self.fluenceSmoothed[int(x),int(y)]
  







