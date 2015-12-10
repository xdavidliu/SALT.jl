#!/usr/bin/python
Nxy = 120*120
Nc = 2  # run3d always had Nc = 3
Nznew = 3
infile = 'above_file.m' # 'aboveA_file.m'
outfile = 'above3d_file.m'
initial_useless_lines = 3    # 3 for mode files
Nr = 2 # for clarity
fin = open(infile)
fout = open(outfile, 'w')
for i in range(initial_useless_lines): # copy first few lines in .m file
	s = fin.readline()
	fout.write(s)
for i in range(Nr*Nxy*Nc):
	s = fin.readline()
	for j in range(Nznew):
		fout.write(s)
	# no need to do below unless Nc == 2
	
	if Nc == 2 and (i+1) % (Nxy*Nc) == 0: 
		for k in range(Nxy*Nznew):
			fout.write('0.0\n') # empty padding for Ez
		fout.flush() # just to be safe. also future writes append correctly
	
s = fin.read() # just read entire rest of file
fout.write(s)
print('must manually put 1e-15 for kz in new .m file! putting it in run script does not work because that is only for passive. Also, remember to rerun this script with Nr = 1, Nc =1 and initial_useless_lines = 0 for the eps and fprof text files!')
