#!bin/bash

filename="dates_folder_all.txt"


#change USERID to your earthdata login (which has podaac permissons enabled)
#change PASSWORD to your earthdata login
# /cygdrive/i/cygnss/data/raw/${yyyy}/ is the data directory for me.  Change as nessary 
# yyyy is the year 
#dates_folder_all.txt has the Julian day is the format of %03d

yyyy=$(printf "%04d" 2019)
while read line ;do

echo $line
wget --http-user=USERID --http-password=PASSWORD --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off -nH --cut-dirs 8 -l1 --no-parent -A "*.nc" https://podaac-opendap.jpl.nasa.gov/opendap/allData/cygnss/L1/v3.0/${yyyy}/$line/ -P ../drive/files/allData/cygnss/L1/Nov/v3.0/${yyyy}/$line/

done < "$filename"

