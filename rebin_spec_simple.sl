
variable filn = __argv[1];
Remove_Spectrum_Gaps=1;
variable id=load_data(filn);
%group(id;min_sn=5,min_chan=10,bounds=1.0,unit="keV");
group(id;min_sn=5,min_chan=2,bounds=1.0,unit="keV");
regroup_file(id);

exit;
