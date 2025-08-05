clearAnnotations()
def entry = getProjectEntry()
def imgName = entry.getImageName().tokenize('.')[0]
print imgName

def gson = GsonTools.getInstance(true)

def json = new File("./Results/BREAST/clustering/bbox_IHC_reg/${imgName}.geojson").text
// Read the annotations
def type = new com.google.gson.reflect.TypeToken<List<qupath.lib.objects.PathObject>>() {}.getType()
def deserializedAnnotations = gson.fromJson(json, type)
addObjects(deserializedAnnotations)

def server = getCurrentServer()
def roi = getAnnotationObjects()[0].getROI()
def requestROI = RegionRequest.createInstance(server.getPath(), 1, roi)
writeImageRegion(server, requestROI, "./BREAST/IHC_reg/${imgName}.ome.tif")
