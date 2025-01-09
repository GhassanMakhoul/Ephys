import mne


edf_file = mne.io.read_raw_edf("Z:/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/Epat21/Epat21_11112019_20211303_clabel.EDF")
edf_file.crop(tmin=0.0, tmax=200.0)

out_fpath = "Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/nights_of_sleep/Epat21/test.edf"

edf_file.export(out_fpath)