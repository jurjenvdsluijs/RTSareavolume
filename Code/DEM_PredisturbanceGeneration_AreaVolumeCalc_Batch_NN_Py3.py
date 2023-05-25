#-------------------------------------------------------------------------------
# Name:         DEM_PredisturbanceGeneration_AreaVolumeCalc_Batch_NN_Py3.py
#
# Purpose:      Takes existing high-resolution DEMs and slump shapefiles (polygons) to model predisturbance elevations for the goal of obtaining DEMs of Difference and slump area/volume estimates.
#               If delineations reflect real disturbances the DOD metrics relate to slump structural estimates, such as volume and average depth-of-thaw.
#               If delineations reflect undisturbed locations the DOD metrics relate to the accuracy and precision to which undisturbed topography can be recreated using neighbouring elevation measurements and Natural Neighbour re-interpolation.
#               Script will work on multiple shapefile input files.
#
# Output        Database file (DBF) created through ESRI ArcGIS Pro "Zonal Statistics" tool, with pixel counts and various DOD metrics (min, max, range, mean, standard deviation, median, 90th percentile distribution, and SUM.
#               The metric "SUM" represents slump volume if the inputs reflect real disturbances and should be a negative value when correct. SUM represents volumetric uncertainty if using synthetic data voids.
#               Outputs can be found in 06_FinalZonalStats and 07_FinalRMSEStats folders. The latter folder is irrelevant if working with real disturbances. If the input vector data concerns multiple shapefiles a DBF file is produced for each shapefile, along with a merged BDF file combining all features.
#               Please refer to the following publication for more information: Van der Sluijs et al., Allometric scaling of retrogressive thaw slumps, Cryosphere Discussions (in review): https://tc.copernicus.org/preprints/tc-2022-149/
#
# Author:       Jurjen van der Sluijs, Unmanned Aircraft Systems Coordinator, NWT Centre for Geomatics, Yellowknife, Northwest Territories, Canada (jurjen_vandersluijs@gov.nt.ca)
#
# Created:      17/07/2019 (Original). Last metadata edits (02/09/2023)
#
# Licence:      CC Attribution-NonCommercial-ShareAlike (CC BY-NC-SA)
#               Remix, tweak, build upon this work is allowed non-commercially and license new creations under identical terms.
#
# Citation:     Developed for the following preprint: Van der Sluijs et al., Allometric scaling of retrogressive thaw slumps, Cryosphere Discussions (in review): https://tc.copernicus.org/preprints/tc-2022-149/
#               Use associated final peer-reviewed article where possible.
#
# Dependencies:
#               • Developed and tested using ESRI ArcGIS Pro v2.7 to v2.9 with ArcPy, Spatial Analyst or 3D Analyst
#               • Developed and tested using Python 3.7.11 [MSC v.1927 64 bit (AMD64)] on win32
#
# Considerations:
#               • Shapefile must have a “UniqueID” attribute.
#               • Shapefile is allowed to have overlapping polygons, which normally provide erroneous ArcGIS Zonal Statistics as Table results as polygons are normally rasterized before raster statistics are generated.
#                 In the provided Python code features are processed iteratively in batch mode, thus this limitation has been overcome.
#               • If the spatial resolution of the input DEM is not 1 m, a multiplication is required in order to derive volumes in cubic metres. For example, for the 2 m resolution ArcticDEM this factor is 2 x 2 = 4.
#                 For a 3 m DEM such as MVAP this factor is 3 x 3 = 9. This functionality is not included in the code and needs to be applied afterwards.
#
# More information:
#                   • Preprint: Van der Sluijs et al., Allometric scaling of retrogressive thaw slumps, Cryosphere Discussions (in review): https://tc.copernicus.org/preprints/tc-2022-149/
#                   • Supplement associated with preprint.

## System definitions and extension sign-outs - DO NOT CHANGE
print("Start initializing")
import sys, os, subprocess, arcpy, glob, os.path, time, math, shutil
print("Import OS SYS ARCPY good")
from subprocess import Popen, PIPE
print("Import SubProcess good")
from arcpy import env
print("Import ENV good")
from arcpy.sa import *
print("Import Spatial Analyst good")
from time import localtime, strftime
print("Import Time good")
print("")

try:
    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        # raise a custom exception
        raise LicenseError
except LicenseError:
    print("Spatial Analyst license is unavailable")

try:
    if arcpy.CheckExtension("3D") == "Available":
        arcpy.CheckOutExtension("3D")
    else:
        # raise a custom exception
        raise LicenseError
except LicenseError:
    print("3D Analyst license is unavailable")

env.overwriteOutput = True
totalstart = time.clock()
print("Initialization Complete. Start time:",time.strftime("%c", time.localtime()))
print("")

# Set Variables - Change as Needed
env.workspace  =    r"C:\Workspace"             # Full path to input shapefile deliniations of slumps
inputDEM = r"C:\Workspace\DEMs\DEM.tif"         # Full path to input high resolution DEM (Geotif), which can represents topography in disturbed or undisturbed state (see Purpose)
DEMres   =    1.0    # in meters                # In Step 5 adjust string "CELLSIZE 1.0" appropriately if the DEM does not have a spatial resolution of 1 m.
IDattribute = 3                                 # The integer location of the column 'UniqueID' in the attribute table, following a n-1 naming convention with Python lists. (E.g., 11th column is the 10th attribute in the list including the standard columns such as FID and Shape. Enter "10" in this case)

# Initiate loop
inputVector = arcpy.ListFeatureClasses()
wspace = env.workspace
print("The following vector files will be considered: " + str(inputVector))

for slumpset in inputVector:

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 1: Buffer slump digitizations by 50 m
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    ## Shorten and reset workspace string for every iteration
    arcpy.env.workspace = wspace

    print(time.strftime("%X", time.localtime()), "  Start slump buffering")
    buf_start = time.clock()

    ### string for the output buffered digitizations folder
    bufFolder_name =          "00_SlumpBuffers"
    newrawpath = os.path.join(env.workspace, bufFolder_name)
    if not os.path.exists(newrawpath): os.makedirs(newrawpath)
    bufFolder_joined = os.path.join(env.workspace, bufFolder_name)

    ### Prepare output file name
    desc = arcpy.Describe(slumpset)
    outputVector = bufFolder_joined + "\\" + desc.baseName + "_50m_buf.shp"

    arcpy.Buffer_analysis(slumpset, outputVector , "50 Meters", "FULL", "ROUND", "NONE")

    ### Report happy ending for buffer
    buf_end = time.clock()
    buf_elapsed = buf_end - buf_start
    buf_elapsedmin = buf_elapsed / 60
    buf_elapsedhr = buf_elapsedmin / 60
    print(time.strftime("%X", time.localtime()), "  buffer ran successfully")
    print("Buffer processing time: " + str(round(buf_elapsed,0)) + " seconds " + " or " + str(round(buf_elapsedmin,1)) + " minutes " + " or " + str(round(buf_elapsedhr,4)) + " hours ")
    print(" ")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Step 2: Clip Input DEM to minimum bounding rectangle of buffered slump and convert to point shapefile
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    print(time.strftime("%X", time.localtime()), "  Start DEM Clip")
    clip_start = time.clock()

    ### string for the output clipped DEM folder
    clipFolder_name =          "01_ClippedDEMs"
    newrawpath = os.path.join(env.workspace, clipFolder_name)
    if not os.path.exists(newrawpath): os.makedirs(newrawpath)
    clipFolder_joined = os.path.join(env.workspace, clipFolder_name)

    ## Prepare slump naming convention by iteratively pulling values from attribute table, then clipping individual DEMs based on FOR-loop
    shapeName = arcpy.Describe(outputVector).shapeFieldName
    fieldnames = [field.name for field in arcpy.ListFields(outputVector)]
    cursor = arcpy.SearchCursor(outputVector,['UniqueID'])
    for row in cursor:

        ## Parse UniqueID for further file naming
        rowclean = row.getValue(fieldnames[IDattribute])       ## Here UniqueID is the n-th field in the attribute table; n-1 naming convention with Python lists. (E.g., 11th column is the 10th attribute in the list)
        slumpname = "_SlumpID_" + str(rowclean)

        ## Retrieve inputDEM raster name from which to concatenate with slump name and output path for output files
        descDEM = arcpy.Describe(inputDEM)
        clipDEMname = desc.baseName + slumpname + ".tif"
        clipDEMoutput = clipFolder_joined + "\\" + clipDEMname

        ## Obtain the extent of each feature
        feat = row.getValue(shapeName)
        extent = feat.extent
        extentstr = str('{} {} {} {}'.format(extent.XMin, extent.YMin, extent.XMax, extent.YMax))

        print(slumpname)
        print(extentstr)

        ## Clip the DEM based on invidual extent
        arcpy.Clip_management(inputDEM, extentstr,clipDEMoutput, "", "-9999", "NONE", "NO_MAINTAIN_EXTENT")
        print(clipDEMname + " successfully clipped")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 3: Convert clipped DEMs to points
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        ## Prepare slump naming convention
        postpointsDEMname = desc.baseName + slumpname + "_points_post"
        postpointsDEMoutput = r"memory" + "\X" + postpointsDEMname

        ## Convert clupped DEMs to points
        arcpy.RasterToPoint_conversion(clipDEMoutput, postpointsDEMoutput, "VALUE")
        print(postpointsDEMname + " successfully processed")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 4: Select by location and remove points
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        ## Prepare slump naming convention for interpreter
        prepointsDEMname = desc.baseName + slumpname + "_points_pre"

        ## Create a temporary layer for selection of points and the respective slump area (cannot rely solely on a Select by Location/intersect in case of overlapping randomized slump areas)
        tempPoints = "pointsLayer"
        tempArea = "areaLayer"
        arcpy.MakeFeatureLayer_management(postpointsDEMoutput, tempPoints)
        arcpy.MakeFeatureLayer_management(slumpset, tempArea)

        ## Select by attribute for correct slump area
        whereClause = '"UniqueID" = ' + str(rowclean)
        print(whereClause)
        arcpy.SelectLayerByAttribute_management(tempArea,'NEW_SELECTION',whereClause)

        ## Select by location and delete points
        arcpy.SelectLayerByLocation_management(tempPoints, "Intersect",tempArea)
        arcpy.DeleteFeatures_management(tempPoints)

        ## Save left-over points to new feature class
        print(prepointsDEMname + " successfully processed")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 5: Reinterpolate clipped points to create pre-disturbance model
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        ### string for the output predisturbance DEM folder
        predisFolder_name =          "02_PredisturbDEMs"
        newrawpath = os.path.join(env.workspace, predisFolder_name)
        if not os.path.exists(newrawpath): os.makedirs(newrawpath)
        predisFolder_joined = os.path.join(env.workspace, predisFolder_name)

        ## Retrieve inputDEM raster name from which to concatenate with slump name and predisturbance DEM
        predisTINname = desc.baseName + slumpname + "_predisturbanceTIN"
        predisTINoutput = predisFolder_joined + "\\" + predisTINname

        predisDEMname = desc.baseName + slumpname + "_predisturbance.tif"
        predisDEMoutput = predisFolder_joined + "\\" + predisDEMname

        ## Concatename in_features string for TIN parameters and describe spatial reference for TIN output
        featuresTIN = postpointsDEMoutput + " GRID_CODE Mass_Points <None>"
        spref = arcpy.Describe(clipDEMoutput).spatialReference

        ## Execute NN interpolation to obtain predisturbance DEM
        arcpy.CreateTin_3d(predisTINoutput,spref, featuresTIN, "DELAUNAY")
        arcpy.TinRaster_3d(predisTINoutput, predisDEMoutput,"FLOAT","NATURAL_NEIGHBORS","CELLSIZE 1.0","1")
        arcpy.Delete_management(postpointsDEMoutput)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 6: Create DOD through Raster Algebra
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        ### string for the output DOD folder
        dodFolder_name =          "03_DODs"
        newrawpath = os.path.join(env.workspace, dodFolder_name)
        if not os.path.exists(newrawpath): os.makedirs(newrawpath)
        dodFolder_joined = os.path.join(env.workspace, dodFolder_name)

        ## Retrieve inputDEM raster name from which to concatenate with slump name and predisturbance DEM to name DOD raster
        dodname = desc.baseName + slumpname + "_dod.tif"
        dodoutput = dodFolder_joined + "\\" + dodname

        ## Execute Raster Algebra for DOD
        outMinus = Minus(clipDEMoutput, predisDEMoutput)
        outMinus.save(dodoutput)
        arcpy.CalculateStatistics_management(dodoutput, "", "", "")

        print(dodname + " successfully processed")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 7: Create Summary Statistics of area and volume through Zonal Statistics
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        ### string for the output Zonal Statistics folder
        zsFolder_name =          "04_IndZonalStats_" + desc.baseName
        newrawpath = os.path.join(env.workspace, zsFolder_name)
        if not os.path.exists(newrawpath): os.makedirs(newrawpath)
        zsFolder_joined = os.path.join(env.workspace, zsFolder_name)

        ## Retrieve inputDEM raster name from which to concatenate with slump name and DOD raster to name Zonal Statistics table
        iTablename = desc.baseName + slumpname + "_zs.dbf"
        iTableoutput = zsFolder_joined + "\\" + iTablename

        ## Select by Attribute on slump name to ensure that only the feature itself is used (closely located slumps cause multiple records for each when merged; next step). Start with temporary layer creation, then concatenation of SQL query string before running zonal statistics tool.

        ## Create a temporary layer for selection
        zoneLayer = "zone"
        arcpy.MakeFeatureLayer_management(slumpset, zoneLayer)

        ## Concatenate query
        qry = ""'UniqueID'" = " + str(rowclean) + " "

        ## Select by Attribute
        arcpy.SelectLayerByAttribute_management(zoneLayer, "NEW_SELECTION",qry)

        ## Run Zonal Statistics as a Table for each individual unbuffered feature
        arcpy.gp.ZonalStatisticsAsTable_sa(zoneLayer, "UniqueID", dodoutput, iTableoutput, "DATA", "ALL")

        print(iTablename + " successfully processed")
        print("")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 8: Calculate RMSE and save as table
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        ### string for the output Zonal Statistics folder
        rmseFolder_name =          "05_IndRMSEStats_" + desc.baseName
        newrawpath = os.path.join(env.workspace, rmseFolder_name)
        if not os.path.exists(newrawpath): os.makedirs(newrawpath)
        rmseFolder_joined = os.path.join(env.workspace, rmseFolder_name)

        ## Retrieve inputDEM raster name from which to concatenate with slump name and DOD raster to name Zonal Statistics table
        rmseTablename = desc.baseName + slumpname + "_rmse.dbf"
        rmseTableoutput = rmseFolder_joined + "\\" + rmseTablename

        ## Retrieve inputDEM raster name from which to concatenate with slump name to name squared DOD raster
        dodsqname = desc.baseName + slumpname + "_dodsq.tif"
        dodsqoutput = dodFolder_joined + "\\" + dodsqname

        ## Execute Raster Algebra for RMSE
        outSquare = Square(dodoutput)
        outSquare.save(dodsqoutput)
        arcpy.CalculateStatistics_management(dodsqoutput, "", "", "")

        print(dodsqname + " successfully processed")

        ## Select by Attribute on slump name to ensure that only the feature itself is used (closely located slumps cause multiple records for each when merged; next step). Start with temporary layer creation, then concatenation of SQL query string before running zonal statistics tool.

        ## Create a temporary layer for selection
        zoneLayer = "zone"
        arcpy.MakeFeatureLayer_management(slumpset, zoneLayer)

        ## Concatenate query
        qry = ""'UniqueID'" = " + str(rowclean) + " "

        ## Select by Attribute
        arcpy.SelectLayerByAttribute_management(zoneLayer, "NEW_SELECTION",qry)

        ## Run Zonal Statistics as a Table for each individual unbuffered feature
        arcpy.gp.ZonalStatisticsAsTable_sa(zoneLayer, "UniqueID", dodsqoutput, rmseTableoutput, "DATA", "MEAN")

        print(rmseTablename + " successfully processed")
        print("")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 9: Create Summary Statistics through Zonal Statistics and merge per iteration
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    ### string for the output Zonal Statistics folder
    fzsFolder_name =          "06_FinalZonalStats"
    newrawpath = os.path.join(env.workspace, fzsFolder_name)
    if not os.path.exists(newrawpath): os.makedirs(newrawpath)
    fzsFolder_joined = os.path.join(env.workspace, fzsFolder_name)

    ## Retrieve a list of tables with statistics of each individual slump
    arcpy.env.workspace = zsFolder_joined
    tableList = arcpy.ListTables()

    print("The following tables will be merged:")
    for table in tableList:
        print(table)

    ## Retrieve input vector name from which to concatenate with final statistics name for merged DBF table output
    fTablename = desc.baseName +  "_FinalStatistics" + ".dbf"
    fTableoutput = fzsFolder_joined + "\\" + fTablename

    ## Merge list of tables in final table output
    arcpy.Merge_management(tableList,fTableoutput)

    print(fTablename + " successfully processed")
    print("")

    ### Report happy ending for iteration
    clip_end = time.clock()
    clip_elapsed = clip_end - clip_start
    clip_elapsedmin = clip_elapsed / 60
    clip_elapsedhr = clip_elapsedmin / 60
    print(time.strftime("%X", time.localtime()), "  iteration ran successfully")
    print("Area/volume processing time for iteration: " + str(round(clip_elapsed,0)) + " seconds " + " or " + str(round(clip_elapsedmin,1)) + " minutes " + " or " + str(round(clip_elapsedhr,4)) + " hours ")
    print(" ")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 10: Create Summary Statistics through Zonal Statistics and merge per iteration
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    ## Shorten and reset workspace string for every iteration
    arcpy.env.workspace = wspace

    ### string for the output Zonal Statistics folder
    frmseFolder_name =          "07_FinalRMSEStats"
    newrawpath = os.path.join(env.workspace, frmseFolder_name)
    if not os.path.exists(newrawpath): os.makedirs(newrawpath)
    frmseFolder_joined = os.path.join(env.workspace, frmseFolder_name)

    ## Retrieve a list of tables with statistics of each individual slump
    arcpy.env.workspace = rmseFolder_joined
    rmsetableList = arcpy.ListTables()

    print("The following tables will be merged:")
    for rmsetable in rmsetableList:
        print(rmsetable)

    ## Retrieve input vector name from which to concatenate with final statistics name for merged DBF table output
    frmseTablename = desc.baseName +  "_FinalRMSE" + ".dbf"
    frmseTableoutput = frmseFolder_joined + "\\" + frmseTablename

    ## Merge list of tables in final table output
    arcpy.Merge_management(rmsetableList,frmseTableoutput)

    print(frmseTablename + " successfully processed")
    print("")

    ### Report happy ending for iteration
    clip_end = time.clock()
    clip_elapsed = clip_end - clip_start
    clip_elapsedmin = clip_elapsed / 60
    clip_elapsedhr = clip_elapsedmin / 60
    print(time.strftime("%X", time.localtime()), "  iteration ran successfully")
    print("Area/volume processing time for iteration: " + str(round(clip_elapsed,0)) + " seconds " + " or " + str(round(clip_elapsedmin,1)) + " minutes " + " or " + str(round(clip_elapsedhr,4)) + " hours ")
    print(" ")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 11: Merge zonal statistics of all iterations together
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Retrieve a list of tables with statistics of each individual slump
arcpy.env.workspace = fzsFolder_joined
fsTableList = arcpy.ListTables()

print("The following tables will be merged:")
print(fsTableList)

## Retrieve input vector name from which to concatenate with final statistics name for merged DBF table output
fsTablename = desc.baseName +  "_FinalStatistics_merged" + ".dbf"
fsTableoutput = fzsFolder_joined + "\\" + fsTablename

## Merge list of tables in final table output
arcpy.Merge_management(fsTableList,fsTableoutput)

print(fsTablename + " successfully processed")
print("")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 12: Merge RMSE statistics of all iterations together
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Retrieve a list of tables with statistics of each individual slump
arcpy.env.workspace = frmseFolder_joined
fmrmseTableList = arcpy.ListTables()

print("The following tables will be merged:")
print(fmrmseTableList)

## Retrieve input vector name from which to concatenate with final statistics name for merged DBF table output
fmrmseTablename = desc.baseName +  "_FinalRMSE_merged" + ".dbf"
fmrmseTableoutput = frmseFolder_joined + "\\" + fmrmseTablename

## Merge list of tables in final table output
arcpy.Merge_management(fmrmseTableList,fmrmseTableoutput)

print(fmrmseTablename + " successfully processed")
print("")

### Report happy ending for iteration
totalend = time.clock()
total_elapsed = totalend - totalstart
total_elapsedmin = total_elapsed / 60
total_elapsedhr = total_elapsedmin / 60
print(time.strftime("%X", time.localtime()), "  script ran successfully")
print("Area/volume processing time for script: " + str(round(total_elapsed,0)) + " seconds " + " or " + str(round(total_elapsedmin,1)) + " minutes " + " or " + str(round(total_elapsedhr,4)) + " hours ")
print(" ")