import glob
output_file = open("output.csv", "w")
output_file.write('filename'+', '+'V(G)'+', '+'V(e)'+', '+'Vp'+', '+'V(G)/Vp'+', '+'logL'+', '+ 'logL0'+', '+ 'LRT'+', '+ 'Pval'+', '+ 'n' +'\n')
for file in glob.glob('*.hsq'):
    input_file    = file
    with open(input_file) as f:
        variables = []
        content = f.readlines()
        for line in content:
            line = line.split()
            variables.append(line[1])
        VG, Ve, Vp, VG_Vp, logL, logL0, LRT, Pval, n = variables[1:]  
        output_file.write(input_file +', '+ VG+', '+ Ve +', '+ Vp +', '+ VG_Vp +', '+ logL +', '+ logL0 +', '+ LRT +', '+ Pval +', '+ n +'\n')     
output_file.close()