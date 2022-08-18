#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('getwd()')


# In[2]:


get_ipython().system('ls')


# In[3]:


import zipfile
with zipfile.ZipFile('Updated_version.zip', 'r') as zip_ref:
    zip_ref.extractall()


# In[6]:


import math, time, geopandas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from dateutil import parser


# In[7]:


#Input start and end dates for evaluation
start_date = "2022-03-15 00:00"
end_date = "2022-07-27 00:00"


# In[12]:


# # Reading data
EQuIS_Analytical_Results = "Input_EQuIS_Analytical_Results.csv"
EQuIS_Location_Parameters = "Input_EQuIS_Location_Parameters.csv"
location_information = "Input_location_information.csv"


# In[13]:


analytical_results_indexes = ['user_report_id','facility_id','facility_code','sys_loc_code','loc_name','sample_id','sys_sample_code','sample_name','sample_date','sampledate','sampletime','sample_type_code','start_depth','end_depth','depth_unit','matrix_code','task_code','task_code_2','parent_sample_code','field_sdg','analysis_location','lab_sample_id','lab_matrix_code','lab_name_code','analytic_method','analysis_date','column_number','fraction','test_type','prep_method','leachate_method','leachate_date','lab_sdg','percent_moisture','dilution_factor','test_id','cas_rn','chemical_name','organic_yn','report_result_text','report_result_value','report_result_unit','report_result_limit','report_method_detection_limit','report_reporting_limit','report_quantitation_limit','reportable_result','detect_flag','interpreted_qualifiers','validator_qualifiers','lab_qualifiers','quantitation_limit','method_detection_limit','reporting_detection_limit','detection_limit_unit','approval_code','result_text','result_numeric','result_unit','result_type_code','x_coord','y_coord','z_coord_avg','zfrom','zto','longitude','latitude','method_analyte_group','mag_report_order','locparam_submergestatus','remark','locparam_remark','locparam_param_value','locparam_param_code','locparam_systemstatus','locparam_param_unit','loc1','loc2','coord_type_code','coord_unit','well_id','well_status','dt_result_approval_code','well_purpose','subfacility_code','validated_yn','coord_desc','srid','esri_spatial_ref','analysis_batch','prep_batch','batch_validator_name','geologic_unit_code','custom_field_2','custom_field_3','facility_full_name','well_owner','reporting_units','unit_conversion_note','ID','loc_group_code','loc_report_order']
lp_cols_headers = ['ID','user_report_id','facility_id','facility_code','sys_loc_code','loc_name','loc_group','loc_report_order','measurement_date','param_code','param_value','param_text','param_unit','measurement_method','task_code','x_coord','y_coord','longitude','latitude']


# In[14]:


#Read EQuIS inputs
df_EQuIS_Analytical_Results = pd.read_csv(EQuIS_Analytical_Results , names= analytical_results_indexes,low_memory = False) #after resetting header
df_EQuIS_Location_Parameters = pd.read_csv(EQuIS_Location_Parameters,names = lp_cols_headers,low_memory=False)


# In[15]:


# # Start of script

# ## Set variables

# In[ ]:


#Set values for step, number of iterations on each point, and distance tolerance to biovent points
step_default = 2
iterations = 100
sub_iterations = 2
tolerance = 7

#Vadose zone attributes
n = 0.2 #porosity (no units)
H = 1 #thickness (feet)

#Text string IDs for each system area
system_areas = ["A03", "A3", "A5"]


# In[16]:



#Spatial reference as well-known text (WKT)
sr = str('''PROJCS["NAD_1983_StatePlane_Illinois_West_FIPS_1202_Feet",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",2296583.333333333],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",-90.16666666666667],PARAMETER["Scale_Factor",0.9999411764705882],PARAMETER["Latitude_Of_Origin",36.66666666666666],UNIT["Foot_US",0.30480060960121924]]''')


# ## Define Functions


# In[17]:



#Function used to calculate psi, Q, W, and coordinate complexes for points
def stream(point_eval, df_biovents):
   
    #Set initial variables from agruments
    point_x = point_eval.x
    point_y = point_eval.y
    point_complex = complex(point_x,point_y)
    stream_psi = 0
    stream_Qx = 0
    stream_Qy = 0
   
    #Calculate psi and W
    for index_bv, row_bv in df_biovents.iterrows():
        x_bv = df_biovents.loc[index_bv, "x_coord"]
        y_bv = df_biovents.loc[index_bv, "y_coord"]
        af_bv = df_biovents.loc[index_bv, "param_value"]

        #Calculate psi
        theta = (point_y-y_bv)/(point_x-x_bv)
        stream_psi = stream_psi + (af_bv/(2*math.pi)*(math.atan(theta)))

        #Calculate Qx and Qy
        Qx_new = (af_bv/(2*math.pi))*(point_x - x_bv)/((point_x-x_bv)**2+(point_y - y_bv)**2)
        Qy_new = (af_bv/(2*math.pi))*(point_y - y_bv)/((point_x-x_bv)**2+(point_y - y_bv)**2)
        stream_Qx = stream_Qx + Qx_new
        stream_Qy = stream_Qy + Qy_new
       
    #Calculate W
    stream_W = complex(stream_Qx, stream_Qy)
   
    return stream_psi, stream_Qx, stream_Qy, stream_W, point_complex


# In[18]:



#Function to advance step for points during the iteration
def step_advance(ad_Qx, ad_Qy, ad_W, ad_complex):
    ad_x = ((step*ad_Qx)/abs(ad_W) + ad_complex.real)
    ad_y = ((step*ad_Qy)/abs(ad_W) + ad_complex.imag)
    ad_point = Point(ad_x,ad_y)
    return ad_point

#Function to calculate W using Newton-Rhapsom
def W_NewtonRhapson(wnr_Qx, wnr_Qy, wnr_psi, wnr_psi_vmp):
    wnr_W = complex(
        ((wnr_Qx)*(wnr_psi-(wnr_psi_vmp))/((wnr_Qx**2)+(wnr_Qy**2))),
        ((wnr_Qy)*(wnr_psi-(wnr_psi_vmp))/((wnr_Qx**2)+(wnr_Qy**2)))
    )
    return wnr_W


# In[19]:



#Function to check if points cross a discontinuity and create coordinate modifier as needed.
def discontinuity_check(point_eval_prev, point_eval_post, df_biovent):
    #Set variables from agruments
    dc_x_prev = point_eval_prev.x
    dc_y_prev = point_eval_prev.y
    dc_x_post = point_eval_post.x
    dc_y_post = point_eval_post.y
    modifier_x = 0
    modifier_y = 0
   
    #Filter BV dataframe for any records that have an x coordinate between the previous and post evaluation points
    df_discon = df_biovent.copy(deep = True)
    df_discon = df_discon[((df_biovent["x_coord"] <= dc_x_post) & (df_discon["x_coord"] > dc_x_prev)) | ((df_biovent["x_coord"] > dc_x_post) & (df_discon["x_coord"] <= dc_x_prev))]
    df_discon.reset_index(inplace=True)
   
    #Check if the new evaluation point crosses a biovent point.
    if df_discon.empty:
        pass
    else:
        if df_discon.shape[0] > 1:
            pass #Need to determine how to address if there are multiple records that meet the criteria of a discontinuity
        else:
            #Set biovent information to variables
            biovent_x = df_discon.loc[0,"x_coord"]
            biovent_af = df_discon.loc[0,"param_value"]
           
            #Determine x modifier
            if dc_x_post >= biovent_x:
                modifier_x = -(biovent_af/2)
            else:
                modifier_x = (biovent_af/2)
    #Filter BV dataframe for any records that have an x coordinate between the previous and post evaluation points
    df_discon = df_biovent.copy(deep = True)
    df_discon = df_discon[((df_discon["y_coord"] <= dc_y_post) & (df_discon["y_coord"] > dc_y_prev)) | ((df_discon["y_coord"] > dc_y_post) & (df_discon["y_coord"] <= dc_y_prev))]
    df_discon.reset_index(inplace=True)
    
        #Check if the new evaluation point crosses a biovent point.
    if df_discon.empty:
        pass
    else:
        if df_discon.shape[0] > 1:
            pass #Need to determine how to address if there are multiple records that meet the criteria of a discontinuity
        else:
            #Set biovent information to variables
            biovent_y = df_discon.loc[0,"y_coord"]
            biovent_af = df_discon.loc[0,"param_value"]
           
            #Determine y modifier
            if dc_y_post >= biovent_y:
                modifier_y = -(biovent_af/2)
            else:
                modifier_y = (biovent_af/2)
           
    return modifier_x, modifier_y


# In[20]:



#Function to calculate flowpath distance and travel time for specific points on the flow path.
def travel(point_eval, point_prior, zone_porosity, zone_thickness, df_biovents):
   
    #Calculate k_phi
    k_phi = 0
    distances_point_eval = df_biovents.distance(point_eval).tolist()
    distances_point_prior = df_biovents.distance(point_prior).tolist()
   
    for i in range(0,len(distances_point_eval),1):
        airflow_bv = df_biovents.loc[i, "param_value"]
        k_phi = k_phi + (airflow_bv/(2*math.pi*zone_thickness))*math.log(distances_point_eval[i]/distances_point_prior[i])
   
    #Calculate travel time
    old_point = geopandas.GeoSeries(point_prior)
    distance_segment = old_point.distance(point_eval).tolist()[0]
    time_segment = ((distance_segment**2)*zone_porosity)/k_phi
   
    return distance_segment, time_segment


# In[21]:



# ## Create dataframes of VMP, biovents, and system readings

# In[ ]:


#Create dataframe of vapor monitoring points (VMP)
#Read input file
df_results = df_EQuIS_Analytical_Results

#Filter results dataframe for VMP results only.
df_vmp_data = df_results[df_results.report_result_unit.eq("%V/V")]

#Reset index
df_vmp_data.reset_index(
    drop = True,
    inplace = True
)

#Cast sample_date to datetime
df_vmp_data = df_vmp_data.rename(
    mapper = {"sample_date": "sample_date_text"},
    axis = 1
)
df_vmp_data.loc[:,"sample_date"] = df_vmp_data["sample_date_text"].apply(lambda x: parser.parse(x))

#Reduce to relevant columns
vmp_col01 = ["sys_loc_code", "sample_date", "chemical_name", "report_result_value", "x_coord", "y_coord"]
df_vmp_data = df_vmp_data[vmp_col01]

#Create a new dataframe of each VMP reading by
#reducing number of columns, dropping duplicate records, and resetting the index
vmp_col02 = ["sys_loc_code", "sample_date", "x_coord", "y_coord"]
df_vmp_readings = df_vmp_data[vmp_col02]
df_vmp_readings = df_vmp_readings.drop_duplicates().reset_index(drop = True)

#Drop records without coordinates
df_vmp_readings.dropna(
    axis = 0,
    subset = ["x_coord"],
    inplace = True
)

#Convert VMP dataframe to geodataframe
df_vmp_readings = geopandas.GeoDataFrame(
    df_vmp_readings,
    geometry=geopandas.points_from_xy(df_vmp_readings.x_coord, df_vmp_readings.y_coord),
    crs = sr
)

#Rename column for consistency with calculation loops
df_vmp_readings = df_vmp_readings.rename({"geometry": "geometry_0"},axis=1)

#Create new dataframe for model iterations by filtering VMP readings for the selected dates.
vmp_event = df_vmp_readings.copy(deep = True)
vmp_event = vmp_event[
    (vmp_event["sample_date"] >= parser.parse(start_date)) & (vmp_event["sample_date"] < parser.parse(end_date))
].reset_index(drop = True)

#Create dataframe of biovent (BV) points
#Read input file
df_bv_all = pd.read_csv(
    filepath_or_buffer = location_information , low_memory=False
)

#Reduce number of columns
bv_col = ["sys_loc_code", "x_coord", "y_coord"]
df_bv_all = df_bv_all[bv_col]

#Filter for BV locations only
df_bv_all = df_bv_all[df_bv_all.sys_loc_code.str.contains("BV")]

#Drop any records without coordinates
df_bv_all.dropna(
    axis = 0,
    subset = ["x_coord"],
    inplace = True
)

#Convert BV dataframe to geodataframe
df_bv_all = geopandas.GeoDataFrame(
    df_bv_all,
    geometry=geopandas.points_from_xy(df_bv_all.x_coord, df_bv_all.y_coord),
    crs = sr
)

#Create dataframe of system airflow readings
#Read input file
df_system = df_EQuIS_Location_Parameters

#Filter for non-zero VMP readings only
df_system = df_system[df_system.sys_loc_code.str.contains("VMP") & ~df_system.param_value.eq(0)]

#Reduce number of columns
system_col = ["sys_loc_code", "measurement_date", "param_value", "x_coord", "y_coord"]
df_system = df_system[system_col]

#Cast measurement_date to datetime
df_system = df_system.rename(
    mapper = {"measurement_date": "measurement_date_text"},
    axis = 1
)
df_system.loc[:,"measurement_date"] = df_system["measurement_date_text"].apply(lambda x: parser.parse(x))

#Create empty dataframe to append any VMP that do not have the required system information for the reading date.
df_vmp_lack = pd.DataFrame(
    columns = ["sys_loc_code", "sample_date", "reason"]
)


# ## Iterate for each VMP record

# In[ ]:


# In[22]:



#For each record in the VMP dataframe, iterate from that point until within the tolerance distance of a BV point
for index, row in vmp_event.iterrows():
    loc = row["sys_loc_code"]
    date = row["sample_date"]
    print(loc + "\t" + str(date))
    df_vmp = vmp_event.copy(deep = True)
    df_vmp = vmp_event[vmp_event.sys_loc_code.eq(loc) & vmp_event.sample_date.eq(date)]
    step = step_default
    flowpath_distance = 0
    flowpath_time = 0
    bv_distances = []
    bv_ends = []
    path_distances = []
    path_times = []
   
    #Check if any of the system designations are in the sys_loc_code value for the VMP location.
    if any(area in loc for area in system_areas):
        system_temp = df_system[df_system.sys_loc_code.eq(loc)]

        #Check if there are system readings associated with the VMP location.
        if system_temp.empty:
            print("\tThere is no recorded system airflow for this sample date. Appending VMP info to separate output for review.")
            df_vmp_lack.loc[len(df_vmp_lack)] = [loc, date, "No recorded system airflow for this sample date."]
        else:
            #Determine in which system the VMP is located.
            for area in system_areas:
                if area in str(row["sys_loc_code"]):
                    #Filter for only the relevant BV locations in the system
                    df_bv = df_bv_all[df_bv_all.sys_loc_code.str.contains(area)]
                    bv_number = df_bv.shape[0]
                else:
                    pass

            #Determine most recent system reading
            print("\tDetermining most recent system reading.")

            #point of error in determining date_system

            date_system = min(system_temp["measurement_date"], key=lambda x: abs(x - date))


            #print(system_temp["measurement_date"])

            #Calculate average airflow for BV locations in the system at the time of the most recent system reading
            print("\tCalculating average airflow for BV locations.")
            system_flow = system_temp.loc[system_temp["measurement_date"].eq(date_system),"param_value"]

            #Check if system_flow != 0. If it is zero, then script cannot model the flow path (because there's no system information)
            #Append VMP name and reading date to separate dataframe for output
            if system_flow.values[0] == '0.0':
                print("\tRecorded system airflow is zero for this sample date. Appending VMP info to separate output for review.")
                df_vmp_lack.loc[len(df_vmp_lack)] = [loc, date, "Recorded system flow is zero for this date."]
            else:
                airflow = float( float(system_flow.values[0]) / bv_number) #Divide system airflow by the number of BV locations

                #Set air flow value to BV locations
                print("\tSetting air flow value to BV locations.")
                df_bv = df_bv.assign(
                    param_value = -airflow
                )

                #Reset index
                df_bv.reset_index(
                    drop = True,
                    inplace = True
                )

                print("\tBegin iterations.")
                for i in range(0,iterations,1):
                    print("\t\t"+str(i))
                    #Set point variable from the dataframe
                    point = df_vmp.loc[index, "geometry_"+str(i)]

                    #Check distance to BV points
                    measurements = df_bv.distance(point).tolist()

                    #If point is within the specified tolerance distance to a BV point, then append the point to the main dataframe
                    print("\t\tDistance to nearest BV: " + str(round(min(measurements),3)))
                    if min(measurements) <= tolerance:
                        print("\t\t\tWithin set tolerance of BV!")
                        print("\t\t\t\tDistance: " + str(round(flowpath_distance,2)))
                        print("\t\t\t\tTime: " + str(round(flowpath_time,2)))

                        bv_distances.append(min(measurements))
                        bv_ends.append(df_bv.loc[measurements.index(min(measurements)),"sys_loc_code"])
                        path_distances.append(flowpath_distance)
                        path_times.append(flowpath_time)

                        break
                    else:
                        #If this is the original VMP, then set the initial psi value to be used in all iterations.
                        if i == 0:

                            #Set initial VMP values
                            psi_vmp, Qx_vmp, Qy_vmp, W_vmp, complex_vmp = stream(point,df_bv)

                            print("\t\tSetting initial VMP values.")
                            print("\t\t\tpsi_vmp: " + str(psi_vmp))
                            print("\t\t\tQx_vmp: " + str(Qx_vmp))
                            print("\t\t\tQy_vmp: " + str(Qy_vmp))
                            print("\t\t\tW_vmp: " + str(W_vmp))

                            #Step advance
                            psi, Qx, Qy, W, complex_point = stream(point,df_bv)
                            point = step_advance(Qx, Qy, W, complex_point)

                            print("\t\tStep advance.")
                            print("\t\t\tpsi: " + str(psi))
                            print("\t\t\tQx: " + str(Qx))
                            print("\t\t\tQy: " + str(Qy))
                            print("\t\t\tW: " + str(W))

                            #Check if point has crossed a discontinuity, then modify as needed.
                            print("\t\tDiscontinuity check.")
                            point_previous = df_vmp.loc[index, "geometry_"+str(i)]
                            mod_x, mod_y = discontinuity_check(point_previous, point, df_bv)

                            #Modify psi_vmp value
                            if mod_x != 0:
                                psi_vmp = psi_vmp + mod_x
                                print("\t\t\tModify psi value because a discontinuity has been crossed.")
                                print("\t\t\tpsi_vmp: " + str(psi_vmp))
                            else:
                                pass

                        #Otherwise, advance step the point
                        else:
                            #Step advance
                            psi, Qx, Qy, W, complex_point = stream(point,df_bv)
                            point = step_advance(Qx, Qy, W, complex_point)

                            print("\t\tStep advance.")
                            print("\t\t\tpsi: " + str(psi))
                            print("\t\t\tQx: " + str(Qx))
                            print("\t\t\tQy: " + str(Qy))
                            print("\t\t\tW: " + str(W))

                            #Check if point has crossed a discontinuity, then modify as needed.
                            print("\t\tDiscontinuity check.")
                            point_previous = df_vmp.loc[index, "geometry_"+str(i)]
                            mod_x, mod_y = discontinuity_check(point_previous, point, df_bv)

                            #Modify psi_vmp value
                            if mod_x != 0:
                                psi_vmp = psi_vmp + mod_x
                                print("\t\t\tModify psi value because a discontinuity has been crossed.")
                                print("\t\t\tpsi_vmp: " + str(psi_vmp))
                            else:
                                pass

                        print("\t\t"+str(point))
                        print("\t\tBegin subiterations")
                        for j in range (0,sub_iterations,1):
                            print("\t\t\tsub " +str(j))

                            #Set previous point to new variable for discontinuity check
                            point_sub = point

                            #Calculate iteration point
                            psi, Qx, Qy, W, complex_point = stream(point,df_bv)
                            print("\t\t\t\tCalculating iteration point.")
                            print("\t\t\t\t\tpsi: " + str(psi))
                            print("\t\t\t\t\tQx: " + str(Qx))
                            print("\t\t\t\t\tQy: " + str(Qy))
                            print("\t\t\t\t\tW: " + str(W))

                            #Determine complex of next point
                            z = W + complex_point

                            #Create geometry object of the point
                            point = Point(z.real,z.imag)
                            print("\t\t\t\t\t"+str(point))

                            #Check if point has crossed a discontinuity, then modify as needed.
                            print("\t\t\t\tDiscontinuity check.")
                            mod_x, mod_y = discontinuity_check(point_sub, point, df_bv)

                            #Modify psi_vmp value
                            if mod_x != 0:
                                psi_vmp = psi_vmp + mod_x
                                print("\t\t\t\tModify psi value because a discontinuity has been crossed.")
                                print("\t\t\t\tpsi_vmp: " + str(psi_vmp))
                            else:
                                pass

                            #Check distance to BV points
                            measurements = df_bv.distance(point).tolist()

                            #If point is within the specified tolerance distance to a BV point, then append the point to the main dataframe
                            print("\t\t\t\tDistance to nearest BV: " + str(round(min(measurements),2)))
                            if (min(measurements) <= tolerance) or (j == sub_iterations-1):

                                df_vmp.loc[index, "geometry_"+str(i+1)] = 0 #Need to create a "blank" record first since can not directly append geometry objects to the dataframe as a new column
                                df_vmp.loc[index, "geometry_"+str(i+1)] = point

                                #Calculate distance and travel time for segment, then add to total distance and time for the flowpath.
                                fp_dist_seg, fp_time_seg = travel(point, point_previous, n, H, df_bv)
                                flowpath_distance = flowpath_distance + fp_dist_seg
                                flowpath_time = flowpath_time + fp_time_seg

                                #Append iteration values to list for retrieval later
                                bv_distances.append(min(measurements))
                                bv_ends.append(df_bv.loc[measurements.index(min(measurements)),"sys_loc_code"])
                                path_distances.append(flowpath_distance)
                                path_times.append(flowpath_time)

                                if j == sub_iterations-1: #Append the final iteration on a point to the main dataframe, regardless whether it is close to a BV point.
                                    print("\t\t\t\tAppend final iteration to main dataframe.")
                                    print("\t\t\t\t\tDistance: " + str(round(flowpath_distance,2))+"(+" + str(round(fp_dist_seg,2)) + ")")
                                    print("\t\t\t\t\tTime: " + str(round(flowpath_time,2))+"(+" + str(round(fp_time_seg,2)) + ")")
                                else:
                                    pass

                                break
                            else:
                                pass

                    iter = bv_distances.index(min(bv_distances))
                    flowpath_distance = path_distances[iter]
                    flowpath_time = path_times[iter]
                    bv_end = bv_ends[iter]
                    bv_distance = bv_distances[iter]

                    #Append distance and travel time values to associated VMP records.
                    print("\tAppend distance and travel time values to associated VMP records.")
                    mask = (vmp_event.sys_loc_code.eq(loc) & vmp_event.sample_date.eq(date))

                    vmp_event.loc[mask, "travel_distance"] = flowpath_distance #Travel distance in feet.
                    vmp_event.loc[mask, "travel_time"] = flowpath_time / 60 #Convert travel time from minutes to hours.
                    vmp_event.loc[mask, "bv_end"] = bv_end #Final BV location
                    vmp_event.loc[mask, "bv_distance"] = bv_distance #Distance to BV in feet.


# In[23]:



#Create dataframe of VMP reading, travel distance, travel time, and final BV information
df_final = vmp_event.copy(deep = True)
df_final = df_final[df_final.travel_time.notnull()].reset_index(drop = True)

#Remove geometry column
df_final.drop(
    labels = "geometry_0",
    axis = 1,
    inplace = True
)

#Export to csv
df_final.to_csv(
    path_or_buf = "results.csv",
    index = False
)

#Export dataframe of VMP that were missing the required system info to be modeled
df_vmp_lack.to_csv(
    path_or_buf = "VMP_results.csv",
    index = False
)


# In[ ]:




