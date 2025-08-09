mkdir LESnew
cd LESnew
ln -s ../run/* .
# Want a hard copy to edit and a backup
cp namelist.input namelist.input0
rm namelist.input
cp namelist.input0 namelist.input
#cp input_sounding input_sounding0
cd ../