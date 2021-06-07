import cygnsslib
import numpy as np 
import datetime as dt 

# cygnss_l1_path = 'C:/Cygnss/drive/files/allData/cygnss/L1/v3.0'

month_to_day_range = {
	'Jan':(1,32),
	'Feb':(32,60),
	'March':(60, 91),
	'April':(91, 121),
	'May':(121,152),
	'June':(152,182),
	'July':(182, 213),
	'August':(213, 244),
	'Sept':(244,274),
	'Oct': (274, 305),
	'Nov': (305, 335),
	'Dec':(335, 366)

}

month = 'Dec'
cygnss_l1_path = 'C:/Cygnss/drive/files/allData/cygnss/L1/'+month+'/v3.0' # Have to change the month here
year = 2019

# days_list = np.arange(1,32)
# days_list = np.arange(32,60) # Have to change day numbers accordingly
month_lower, month_upper = month_to_day_range[month]
days_list = np.arange(month_lower, month_upper)

ref_pos = [30.31183, -98.7759]

radius = 18e3
thresh_ddm_snr = 2
thresh_noise = 1
month_tag = month.lower()
# kml_out_tag = f'feb_quality_rough' # Have to change month name here
kml_out_tag = month_tag+'_quality_rough'
save_podaac_pass = False
quality_check = True

options = {
	'save_cvs': True, # Check this for true or false
	'save_ddm_img': False,
	'plt_tag': '',
	'title_img_imc_ddm_time': False,
	'img_save_type': ['png'],
	'sheet_type': 'xls' 
}

cygnsslib.write_sp_within_radius(cygnss_l1_path, year=2019, daylist = days_list, ref_pos= ref_pos,
	radius=radius, out_root=kml_out_tag, thresh_ddm_snr = thresh_ddm_snr, thresh_noise=thresh_noise, download_cygnss_data = False,
	out_options = options, save_podaac_pass= save_podaac_pass, quality_check = quality_check)