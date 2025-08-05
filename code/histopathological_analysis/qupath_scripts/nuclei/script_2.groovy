def entry = getProjectEntry()
def imgName = entry.getImageName().tokenize('.')[0]
print imgName

def gson = GsonTools.getInstance(true)

def json = new File("./Results/BREAST/nuclei_IHC/nuclei_seg/${imgName}.geojson").text
// Read the annotations
def type = new com.google.gson.reflect.TypeToken<List<qupath.lib.objects.PathObject>>() {}.getType()
def deserializedAnnotations = gson.fromJson(json, type)
addObjects(deserializedAnnotations)
fireHierarchyUpdate()
resolveHierarchy()