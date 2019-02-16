##
## SinglepartToMultipart.py
##
## Sm Jones
## June 2007
##

import arcgisscripting, os, sys, string, math, fpformat 

try:
     
    # Create the geoprocessing object
    gp = arcgisscripting.create ()
    gp.overwriteoutput = 1

    gp.addmessage ("")
    gp.addmessage ("********************************************")
    gp.addmessage ("*                                          *")
    gp.addmessage ("* Convert Singlepart to Multipart features *")
    gp.addmessage ("*                                          *")
    gp.addmessage ("*       SinglepartToMultipart.py           *")
    gp.addmessage ("*                                          *")    
    gp.addmessage ("********************************************")
    gp.addmessage ("*Sm Jones                                  *")
    gp.addmessage ("*June 2007                                 *")
    gp.addmessage ("********************************************")
    gp.addmessage ("")

    # Fetch the XYZ File
    strFile = string.replace(gp.getparameterastext(0),"\\", "/")
    #Set the workspace
    listFileparts = string.split(strFile,"/")
    gp.addmessage ("Crap") 
    for f in listFileparts:
        fn = f
    strFolder = string.replace(strFile, fn, "")
    gp.workspace = strFolder
    newfn = string.replace(fn,".xyz",".shp")

    #Create multipart feature class and addfield
    gp.createfeatureclass_management (gp.workspace, string.replace(newfn,".shp",""),  "MULTIPOINT",  "", "SAME_AS_TEMPLATE", "ENABLED")

    #Declare the insert cursor
    gp.addmessage("Declare a feature cursor ...")
    rows = gp.insertcursor(strFolder + newfn)

    gp.addmessage (strFolder + newfn)    

    #Open the file for proceesingnewfn
    gp.addmessage("Open the xyz file for processing ...")
    xyzfile = open(strFile, "r")
    xyz = xyzfile.readline() 
    xyz = xyzfile.readline()  

    # List values and make feature
    gp.addmessage ("Processing Multi Points ...")
    n = 0
    m = 0
    for xyz in xyzfile:      

        #Make a feature
        if n == 0:
            row = rows.newrow()
            mp = gp.createobject("Array")      

        #Process the points
        listXYZ = string.split(xyz,",")
        test = 0
        for pt in listXYZ:
            test = test + 1
##            if test == 3:
##                gp.addmessage(str(len(listXYZ[2])))
        if n < 3500 and test == 3 and len(listXYZ[2]) > 1:
            pt = gp.createobject("Point")
            pt.x = string.atof(listXYZ[0])
            pt.y = string.atof(listXYZ[1])
            pt.z = string.atof(listXYZ[2])
            mp.add(pt)

        #Insert the shape and reset thre counter
        if n == 3499:
            gp.addmessage("Insert the shape ...")
            row.shape = mp
            rows.insertrow(row)
            n = -1

##        gp.addmessage(str(n))
        n = n + 1       
            
      
    #Increment the counter
    m = m + 1            

    #Insert remianing shape of less than 3500 points         
    row.shape = mp
    rows.insertrow(row)


    xyzfile.close()  
    gp.addmessage (str(m) + " features")
    gp.refreshcatalog (gp.workspace)

    gp.addmessage ("")
    gp.addmessage ("***************************")
    gp.addmessage ("* Susan is seriously cool *")
    gp.addmessage ("***************************")
    gp.addmessage ("")

    #Delete the geoprocessor
    del gp

except:

    # Error encountered
    gp.addmessage ("")
    gp.addmessage ("Error encountered ...")
    gp.addmessage ("")
    del gp
