def entry = getProjectEntry()
def imgName = entry.getImageName().tokenize('.')[0]
print imgName

importObjectsFromFile("annotations_directory/${imgName}.geojson")
for (annotation in getAnnotationObjects()) 
{annotation.setName("Roi_pathologist")}

