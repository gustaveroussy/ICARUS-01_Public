def entry = getProjectEntry()
print entry.getImageName()
def imgName = entry.getImageName().tokenize('.')[0]
def path = buildFilePath(PROJECT_BASE_DIR, 'detection_measurements')
mkdirs(path)
path = buildFilePath(path, imgName + '.csv')
//saveAnnotationMeasurements(path)
saveDetectionMeasurements(path)
print 'Results exported to ' + path
