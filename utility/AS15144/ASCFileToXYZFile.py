#######################################################
#                 ASC File ToXYZ File                 #
# Convert an Ascii file (.asc) into an xyz (.xyz file #
#                                                     #
#                      SM Jones                       #
#             Eagle Technology Group Ltd              #
#                      June 2007                      #
#                   smj@eagle.co.nz                   #
#######################################################
#Import modules and declare the geoprocessor object
##import win32com.client, sys, os, string, math
##gp = win32com.client.Dispatch("esriGeoprocessing.GpDispatch.1")
import arcgisscripting, sys, os, string, fpformat
gp = arcgisscripting.create()
gp.overwriteoutput = 1

#Print message banner
gp.addmessage ("")
gp.addmessage ("#######################################################")
gp.addmessage ("#                 ASC File ToXYZ File                 #")
gp.addmessage ("# Convert an Ascii file (.asc) into an xyz (.xyz file #")
gp.addmessage ("#                                                     #")
gp.addmessage ("#                      SM Jones                       #")
gp.addmessage ("#             Eagle Technology Group Ltd              #")
gp.addmessage ("#                      June 2007                      #")
gp.addmessage ("#                   smj@eagle.co.nz                   #")
gp.addmessage ("#######################################################")
gp.addmessage ("")

########################################################################
#                    USER MUST MODIFY PARAMETERS                       #
########################################################################
# strWorkspace - location of shapefile                                 #
# strFileNames - source LAS files                                      #
# strShapeFile - name of shapefile to be created                       # 
# strProduct - ArcView, ArcEditor or ArcInfo. Not everyone has ArcInfo #
########################################################################

try:

    # Collect the input paramers
    strAsciiFile = string.replace(gp.getparameterastext(0),"\\","/")

    gp.addmessage (strAsciiFile)    

    # Get the outputfile
    gp.addmessage ("")    
    strFileNames = string.split (strAsciiFile, "/")
    for fn in strFileNames:
        strOldFile = fn
    strFolder = string.replace(strAsciiFile, strOldFile,"")
    strSplitfn = string.split(strOldFile, ".")
    strXYZFile = strFolder + strSplitfn[0] + ".xyz"

    gp.addmessage (strXYZFile)    

    # Open the ascii file
    asciifile = open(strAsciiFile, "r")
    #Get the header file
    #rownum
    coordinates = asciifile.readline()
    rowsplit = string.split(coordinates," ")
    for fn in rowsplit:
        rownum = fn
    gp.addmessage ("Number of rows " + str(rownum))

    #colnum
    coordinates = asciifile.readline()
    rowsplit = string.split(coordinates," ")
    for fn in rowsplit:
        colnum = fn
    gp.addmessage ("Number of columns " + str(colnum))

    #topX
    coordinates = asciifile.readline()
    rowsplit = string.split(coordinates," ")
    for fn in rowsplit:
        topX = fn
##    gp.addmessage ("Top X " + str(topX))

    #topY
    coordinates = asciifile.readline()
    rowsplit = string.split(coordinates," ")
    for fn in rowsplit:
        topY = fn
##    gp.addmessage ("Top Y " + str(topY))    

    #Cellsize
    coordinates = asciifile.readline()
    rowsplit = string.split(coordinates," ")
    for fn in rowsplit:
        cellsize = fn
    gp.addmessage ("Cellsize " + str(cellsize))
    
    #NoData
    coordinates = asciifile.readline()
    rowsplit = string.split(coordinates," ")
    for fn in rowsplit:
        noData = fn
##    gp.addmessage ("No data " + str(noData))

    
    #Open the xyz file for writing
    xyzfile = open(strXYZFile, "w")
    strLine = "X,Y,Z\n"
    xyzfile.write (strLine)
    xyzfile.close

    #Convert to float
    topX = string.atof(topX)
    topY = string.atof(topY)
    cellsize = string.atof(cellsize)
    

    #Go through the Ascii file
    allpoints = 0
    rowno = 0
    colno = 0
    gp.addmessage ("Processing " + strXYZFile + " ...")
##    coordinates = asciifile.readline()
    for coordinates in asciifile:
      
        splitCoordinates = string.split(coordinates, " ")

        for z in splitCoordinates:
            
            #Cycle through the matrix            
            if rowno < string.atoi(rownum):
                
                #Get the x
                x = cellsize * rowno
                x = topX + x
                    
                #Get the y
                y = cellsize * colno
                y = topY + y
                
                #Write the string for the output
                strLine = str(x) + "," + str(y) + "," + str(z) + "\n"
                xyzfile.write (strLine)
                allpoints = allpoints + 1
                rowno = rowno + 1
                
            else:

                rowno = 0
                colno = colno + 1

            #Commit records in batches of 1000
            if allpoints % 1000 == 0:
                xyzfile.close ()
                gp.addmessage (str(allpoints) + " features processed ...")
                xyzfile = open(strXYZFile, "a")              

    gp.addmessage (str(allpoints) + " features processed.\n")

    #Close the fike
    asciifile.close()
    xyzfile.close ()

    
    #Refresh the catalog and print the status
    gp.addmessage ("")
    gp.addmessage ("*****************************************************")
    gp.addmessage ("*       All XYZ Points Successfully Converted*      *")
    gp.addmessage ("*                                                   *")
    gp.addmessage ("*              Susan is seriously cool              *")
    gp.addmessage ("*****************************************************")
    gp.addmessage ("")

    #Release the geoprocessor
    del gp

except:

    #Print error message
    gp.addmessage ("")
    gp.addmessage ("************************")
    gp.addmessage ("* Error encountered... *")
    gp.addmessage ("************************")    
    gp.addmessage ("")
    asciifile.close
    del gp
